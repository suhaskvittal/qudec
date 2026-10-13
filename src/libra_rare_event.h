/*
 * Exact fixed-count sampling for independent DEM error instructions.
 *
 * The tree stores truncated, normalized Poisson-binomial count distributions.
 * For a requested total k, each split samples j with weight
 * P_left(j) * P_right(k-j), yielding the exact DEM Bernoulli distribution
 * conditioned on exactly k instructions firing.
 */
#ifndef QUDEC_LIBRA_RARE_EVENT_H
#define QUDEC_LIBRA_RARE_EVENT_H

#include "sampler.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

class LibraRareEventSampler
{
public:
    LibraRareEventSampler(const stim::DetectorErrorModel& dem, size_t max_count)
        : max_count_(validate_max_count(max_count)),
          count_probabilities_(compute_probability_polynomial(dem, max_count_))
    {
        std::vector<double> probabilities;
        probabilities.reserve(dem.count_errors());
        dem.iter_flatten_error_instructions([&probabilities](const auto& inst) {
            probabilities.push_back(inst.arg_data[0]);
        });
        nodes_.reserve(probabilities.empty() ? 1 : 2 * probabilities.size() - 1);
        root_ = build_tree(probabilities, 0, probabilities.size());
        check_count_distribution();
    }

    const std::vector<double>& count_probabilities() const noexcept
    {
        return count_probabilities_;
    }

    // Returns sorted, flattened DEM error-instruction IDs.
    std::vector<size_t> sample_errors(size_t count, std::mt19937_64& rng) const
    {
        if (count >= max_count_)
            throw std::out_of_range("LibraRareEventSampler: count exceeds truncated range");
        if (count >= nodes_[root_].probabilities.size() ||
            nodes_[root_].probabilities[count] == 0 ||
            count_probabilities_[count] == 0)
            throw std::domain_error("LibraRareEventSampler: requested count has zero probability");

        std::vector<size_t> selected;
        selected.reserve(count);
        sample_node(root_, count, rng, selected);
        return selected;
    }

private:
    struct Node
    {
        size_t begin{};
        size_t end{};
        size_t left{std::numeric_limits<size_t>::max()};
        size_t right{std::numeric_limits<size_t>::max()};
        std::vector<double> probabilities;
        bool leaf() const noexcept { return end - begin == 1; }
    };

    size_t max_count_;
    std::vector<double> count_probabilities_;
    std::vector<Node> nodes_;
    size_t root_{};

    static size_t validate_max_count(size_t max_count)
    {
        if (max_count == 0)
            throw std::invalid_argument("LibraRareEventSampler: max_count must be positive");
        return max_count;
    }

    size_t build_tree(const std::vector<double>& p, size_t begin, size_t end)
    {
        Node node;
        node.begin = begin;
        node.end = end;
        if (end - begin == 1)
        {
            const double q = p[begin];
            if (!(q >= 0 && q <= 1))
                throw std::invalid_argument("LibraRareEventSampler: invalid DEM error probability");
            node.probabilities.assign(std::min<size_t>(2, max_count_), 0.0);
            node.probabilities[0] = 1.0 - q;
            if (max_count_ > 1)
                node.probabilities[1] = q;
            normalize(node.probabilities);
        }
        else if (begin == end)
        {
            node.probabilities = {1.0};
        }
        else
        {
            const size_t mid = begin + (end - begin) / 2;
            node.left = build_tree(p, begin, mid);
            node.right = build_tree(p, mid, end);
            const auto& a = nodes_[node.left].probabilities;
            const auto& b = nodes_[node.right].probabilities;
            const size_t length = std::min(max_count_, a.size() + b.size() - 1);
            node.probabilities.assign(length, 0.0);
            for (size_t i = 0; i < a.size(); ++i)
                for (size_t j = 0; j < b.size() && i + j < length; ++j)
                    node.probabilities[i + j] += a[i] * b[j];
            normalize(node.probabilities);
        }
        nodes_.push_back(std::move(node));
        return nodes_.size() - 1;
    }

    static void normalize(std::vector<double>& values)
    {
        double sum = 0;
        for (double v : values)
            sum += v;
        if (sum > 0)
            for (double& v : values)
                v /= sum;
    }

    void check_count_distribution() const
    {
        const auto& tree = nodes_[root_].probabilities;
        double reference_sum = 0;
        for (double p : count_probabilities_)
            reference_sum += p;
        if (!(reference_sum > 0))
            throw std::runtime_error("LibraRareEventSampler: truncated count distribution has zero mass");
        for (size_t k = 0; k < max_count_; ++k)
        {
            const double expected = count_probabilities_[k] / reference_sum;
            const double actual = k < tree.size() ? tree[k] : 0;
            if (std::abs(expected - actual) > 2e-10 * std::max(1.0, std::abs(expected)))
                throw std::runtime_error("LibraRareEventSampler: tree and sampler count distributions disagree");
        }
    }

    void sample_node(size_t node_id, size_t count, std::mt19937_64& rng,
                     std::vector<size_t>& selected) const
    {
        const Node& node = nodes_[node_id];
        if (count == 0)
            return;
        if (node.leaf())
        {
            if (count != 1)
                throw std::logic_error("LibraRareEventSampler: invalid leaf count");
            selected.push_back(node.begin);
            return;
        }

        const Node& left = nodes_[node.left];
        const Node& right = nodes_[node.right];
        const size_t low = count > right.end - right.begin ? count - (right.end - right.begin) : 0;
        const size_t high = std::min({count, left.end - left.begin, left.probabilities.size() - 1});
        double total = 0;
        for (size_t left_count = low; left_count <= high; ++left_count)
        {
            const size_t right_count = count - left_count;
            if (left_count < left.probabilities.size() &&
                right_count < right.probabilities.size())
                total += left.probabilities[left_count] * right.probabilities[right_count];
        }
        if (!(total > 0))
            throw std::runtime_error("LibraRareEventSampler: conditional split has zero mass");
        std::uniform_real_distribution<double> uniform(0.0, total);
        double choice = uniform(rng);
        size_t left_count = low;
        for (size_t i = low; i <= high; ++i)
        {
            const size_t right_count = count - i;
            const double weight = i < left.probabilities.size() &&
                                  right_count < right.probabilities.size()
                ? left.probabilities[i] * right.probabilities[right_count]
                : 0.0;
            if (choice < weight || i == high)
            {
                left_count = i;
                break;
            }
            choice -= weight;
        }
        sample_node(node.left, left_count, rng, selected);
        sample_node(node.right, count - left_count, rng, selected);
    }
};

#endif // QUDEC_LIBRA_RARE_EVENT_H
