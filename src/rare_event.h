#ifndef QUDEC_RARE_EVENT_H
#define QUDEC_RARE_EVENT_H

#include <stim.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

std::vector<double> compute_probability_polynomial(const stim::DetectorErrorModel&, size_t);

// Independent DEM Bernoulli channels conditioned on their total active count.
// At each tree split, P(left=j, right=k-j) determines the conditional choice.
// Normalization of truncated node polynomials cancels in these choices.
class ConditionalErrorSampler {
public:
    ConditionalErrorSampler(const stim::DetectorErrorModel& dem, size_t max_count)
        : max_count_(max_count),
          num_detectors_(dem.count_detectors()), num_observables_(dem.count_observables()),
          count_probabilities_(compute_probability_polynomial(dem, max_count)) {
        if (!max_count) throw std::invalid_argument("ConditionalErrorSampler: max_count must be positive");
        probabilities_.reserve(dem.count_errors());
        dem.iter_flatten_error_instructions([&](const auto& inst) {
            const double p = inst.arg_data[0];
            if (!std::isfinite(p) || p < 0 || p > 1)
                throw std::invalid_argument("Invalid DEM error probability");
            probabilities_.push_back(p);
            patterns_.emplace_back(inst.target_data.begin(), inst.target_data.end());
        });
        nodes_.reserve(probabilities_.empty() ? 1 : 2 * probabilities_.size() - 1);
        root_ = build_tree(0, probabilities_.size());
        long double retained = 0;
        for (double p : count_probabilities_) retained += p;
        if (!(retained > 0))
            throw std::domain_error("No probability mass within the retained error counts");
        const auto& tree = nodes_[root_].probabilities;
        for (size_t k = 0; k < max_count_; ++k) {
            const double expected = count_probabilities_[k] / retained;
            const double actual = k < tree.size() ? tree[k] : 0;
            if (std::abs(expected - actual) > 2e-10)
                throw std::runtime_error("Conditional count distributions disagree");
        }
    }

    const std::vector<double>& count_probabilities() const noexcept { return count_probabilities_; }

    std::pair<stim::simd_bit_table<64>, stim::simd_bit_table<64>>
    sample_syndromes(size_t count, size_t shots, std::mt19937_64& rng) const {
        stim::simd_bit_table<64> detectors(shots, num_detectors_);
        stim::simd_bit_table<64> observables(shots, num_observables_);
        for (size_t shot = 0; shot < shots; ++shot)
            for (size_t id : sample_errors(count, rng))
                for (auto t : patterns_[id]) {
                    if (t.is_relative_detector_id()) detectors[shot][t.val()] ^= 1;
                    else if (t.is_observable_id()) observables[shot][t.val()] ^= 1;
                }
        return {std::move(detectors), std::move(observables)};
    }

    // Sorted flattened DEM instruction IDs; impossible counts are rejected.
    std::vector<size_t> sample_errors(size_t count, std::mt19937_64& rng) const {
        if (count >= max_count_) throw std::out_of_range("Error count exceeds retained range");
        const auto& p = nodes_[root_].probabilities;
        if (count >= p.size() || !(p[count] > 0) || !(count_probabilities_[count] > 0))
            throw std::domain_error("Requested error count has zero probability");
        std::vector<size_t> selected;
        selected.reserve(count);
        sample_node(root_, count, rng, selected);
        return selected;
    }

    // Chernoff upper bound on P(K >= first_omitted). Independent of decoder
    // failures, hence also bounds the missing contribution to logical error.
    double tail_probability_bound(size_t first_omitted) const {
        if (!first_omitted) return 1;
        size_t possible = 0;
        long double mean = 0;
        for (double p : probabilities_) { possible += p > 0; mean += p; }
        if (first_omitted > possible) return 0;
        if (first_omitted <= mean) return 1;
        auto tilted_mean = [&](long double theta) {
            long double sum = 0;
            for (double p : probabilities_) {
                if (p == 0) continue;
                if (p == 1) { sum += 1; continue; }
                const long double x = std::log((long double)p) - std::log1p(-(long double)p) + theta;
                sum += x >= 0 ? 1 / (1 + std::exp(-x)) : std::exp(x) / (1 + std::exp(x));
            }
            return sum;
        };
        long double log_bound = 0;
        if (first_omitted == possible) {
            for (double p : probabilities_) if (p > 0) log_bound += std::log((long double)p);
        } else {
            long double low = 0, high = 1;
            while (tilted_mean(high) < first_omitted) high *= 2;
            for (unsigned i = 0; i < 80; ++i) {
                const long double mid = (low + high) / 2;
                if (tilted_mean(mid) < first_omitted) low = mid; else high = mid;
            }
            const long double theta = (low + high) / 2;
            log_bound = -theta * first_omitted;
            for (double p : probabilities_) {
                if (p == 0) continue;
                if (p == 1) { log_bound += theta; continue; }
                const long double a = std::log1p(-(long double)p);
                const long double x = std::log((long double)p) + theta - a;
                log_bound += a + (x > 0 ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x)));
            }
        }
        log_bound += 64 * std::numeric_limits<long double>::epsilon() *
                     (1 + std::abs(log_bound) + probabilities_.size());
        return std::min(1.0, std::nextafter(double(std::exp(log_bound)), 1.0));
    }

private:
    struct Node {
        size_t begin = 0, end = 0, left = 0, right = 0;
        std::vector<double> probabilities;
    };
    size_t max_count_, num_detectors_, num_observables_, root_ = 0;
    std::vector<double> probabilities_, count_probabilities_;
    std::vector<std::vector<stim::DemTarget>> patterns_;
    std::vector<Node> nodes_;

    static void normalize(std::vector<double>& values) {
        long double sum = 0;
        for (double p : values) sum += p;
        if (sum > 0) for (double& p : values) p /= sum;
    }

    size_t build_tree(size_t begin, size_t end) {
        Node node{begin, end, 0, 0, {}};
        if (begin == end) node.probabilities = {1};
        else if (end - begin == 1) {
            node.probabilities = {1 - probabilities_[begin]};
            if (max_count_ > 1) node.probabilities.push_back(probabilities_[begin]);
            normalize(node.probabilities);
        } else {
            const size_t mid = begin + (end - begin) / 2;
            node.left = build_tree(begin, mid);
            node.right = build_tree(mid, end);
            const auto& a = nodes_[node.left].probabilities;
            const auto& b = nodes_[node.right].probabilities;
            node.probabilities.assign(std::min(max_count_, a.size() + b.size() - 1), 0);
            for (size_t i = 0; i < a.size(); ++i)
                for (size_t j = 0; j < b.size() && i + j < node.probabilities.size(); ++j)
                    node.probabilities[i + j] += a[i] * b[j];
            normalize(node.probabilities);
        }
        nodes_.push_back(std::move(node));
        return nodes_.size() - 1;
    }

    void sample_node(size_t id, size_t count, std::mt19937_64& rng, std::vector<size_t>& selected) const {
        if (!count) return;
        const auto& node = nodes_[id];
        if (node.end - node.begin == 1) {
            if (count != 1) throw std::logic_error("Invalid leaf count");
            selected.push_back(node.begin);
            return;
        }
        const auto& left = nodes_[node.left];
        const auto& right = nodes_[node.right];
        const size_t low = count > right.end - right.begin ? count - (right.end - right.begin) : 0;
        const size_t high = std::min({count, left.end - left.begin, left.probabilities.size() - 1});
        auto weight = [&](size_t j) {
            return count - j < right.probabilities.size()
                ? left.probabilities[j] * right.probabilities[count - j] : 0;
        };
        double total = 0;
        for (size_t j = low; j <= high; ++j) total += weight(j);
        if (!(total > 0)) throw std::runtime_error("Conditional split has zero mass");
        double choice = std::uniform_real_distribution<double>(0, total)(rng);
        size_t chosen = low;
        for (size_t j = low; j <= high; ++j) {
            const double p = weight(j);
            if (p <= 0) continue;
            chosen = j; // Last positive bin also handles rounding at the upper endpoint.
            if (choice < p) break;
            choice -= p;
        }
        sample_node(node.left, chosen, rng, selected);
        sample_node(node.right, count - chosen, rng, selected);
    }
};

#endif
