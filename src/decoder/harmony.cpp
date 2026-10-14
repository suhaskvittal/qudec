/* Contributor: OpenAI Codex. Harmony prototype, arXiv:2401.12434. */
#include "harmony.h"
#include <pymatching/sparse_blossom/driver/mwpm_decoding.h>
#include <pymatching/sparse_blossom/driver/user_graph.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>

#if defined(ENABLE_MPI)
#include <mpi.h>
#endif

namespace dec {
namespace {
double sample_probability(double p, double alpha, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> uniform(-alpha, alpha);
    // PyMatching's implied-weight representation assumes nonnegative weights.
    return std::clamp(p * (1.0 + uniform(rng)), 1e-15, 0.499999999999);
}
double probability_from_weight(double weight) {
    return 1.0 / (1.0 + std::exp(weight));
}
void perturb_graph(pm::UserGraph& graph, double alpha_edge, double alpha_rule,
                   std::mt19937_64& rng) {
    for (auto& edge : graph.edges) {
        const double p = sample_probability(edge.error_probability, alpha_edge, rng);
        edge.error_probability = p;
        edge.weight = std::log((1 - p) / p);
        for (auto& rule : edge.implied_weights_for_other_edges) {
            const double q = sample_probability(probability_from_weight(rule.implied_weight),
                                                alpha_rule, rng);
            rule.implied_weight = std::log((1 - q) / q);
        }
    }
}
}

struct Harmony::Impl {
    struct Member {
        std::unique_ptr<pm::Mwpm> first;
        std::unique_ptr<pm::Mwpm> second;
    };
    stim::DetectorErrorModel dem;
    Config config;
    Statistics stats;
    std::vector<Member> members;

    Impl(const stim::DetectorErrorModel& input, Config cfg)
        : dem(input), config(cfg), members(cfg.ensemble_size) {
        if (!cfg.ensemble_size || !cfg.sentinel_size || cfg.sentinel_size > cfg.ensemble_size)
            throw std::invalid_argument("Harmony needs 1 <= sentinel_size <= ensemble_size");
        for (double a : {cfg.alpha_first, cfg.alpha_second, cfg.alpha_conditional})
            if (!(a >= 0 && a <= 1))
                throw std::invalid_argument("Harmony perturbations must lie in [0, 1]");
    }
    void prepare(size_t index) {
        auto& member = members[index];
        if (member.first) return;
        std::mt19937_64 rng(config.seed + 0x9e3779b97f4a7c15ULL * (index + 1));
        auto first = pm::detector_error_model_to_user_graph(dem, true, pm::NUM_DISTINCT_WEIGHTS);
        auto second = pm::detector_error_model_to_user_graph(dem, true, pm::NUM_DISTINCT_WEIGHTS);
        perturb_graph(first, config.alpha_first, 0, rng);
        perturb_graph(second, config.alpha_second, config.alpha_conditional, rng);
        member.first = std::make_unique<pm::Mwpm>(first.to_mwpm(pm::NUM_DISTINCT_WEIGHTS, true));
        member.second = std::make_unique<pm::Mwpm>(second.to_mwpm(pm::NUM_DISTINCT_WEIGHTS, true));
    }
    ObsType run(size_t index, const std::vector<uint64_t>& events, size_t nobs) {
        prepare(index);
        auto& member = members[index];
        std::vector<int64_t> initial_edges;
        pm::decode_detection_events_to_edges(*member.first, events, initial_edges);
        auto& final_graph = member.second->flooder.graph;
        auto& final_search_graph = member.second->search_flooder.graph;
        final_graph.reweight_for_edges(initial_edges);
        final_search_graph.reweight_for_edges(initial_edges);
        std::vector<uint8_t> bytes(nobs);
        pm::total_weight_int ignored_weight = 0;
        try {
            pm::decode_detection_events(*member.second, events, bytes.data(), ignored_weight, false);
        } catch (...) {
            final_graph.undo_reweights();
            final_search_graph.undo_reweights();
            throw;
        }
        final_graph.undo_reweights();
        final_search_graph.undo_reweights();
        ObsType obs(nobs);
        for (size_t k = 0; k < nobs; ++k) obs[k] = bytes[k] & 1;
        return obs;
    }
};

Harmony::Harmony(const stim::DetectorErrorModel& dem) : Harmony(dem, Config{}) {}
Harmony::Harmony(const stim::DetectorErrorModel& dem, Config cfg)
    : num_detectors(dem.count_detectors()), num_observables(dem.count_observables()),
      impl_(std::make_unique<Impl>(dem, cfg)) {
    if (num_observables != 1)
        throw std::invalid_argument("Harmony prototype supports one observable");
}
Harmony::~Harmony() = default;

void Harmony::prepare_ensemble() {
    for (size_t i = 0; i < impl_->config.ensemble_size; ++i)
        impl_->prepare(i);
}

Result Harmony::decode(SyndromeRef syndrome, ObsRef) {
    const auto start = std::chrono::steady_clock::now();
    std::vector<uint64_t> events;
    for (size_t d = 0; d < num_detectors; ++d)
        if (syndrome[d]) events.push_back(d);
    const auto count = [&](size_t end, size_t& ones) {
        for (size_t i = 0; i < end; ++i)
            ones += impl_->run(i, events, num_observables)[0] ? 1 : 0;
        impl_->stats.members_run += end;
    };
    size_t ones = 0;
    const size_t sentinels = impl_->config.layered
        ? impl_->config.sentinel_size : impl_->config.ensemble_size;
    count(sentinels, ones);
    size_t used = sentinels;
    if (impl_->config.layered && ones != 0 && ones != sentinels) {
        for (size_t i = sentinels; i < impl_->config.ensemble_size; ++i)
            ones += impl_->run(i, events, num_observables)[0] ? 1 : 0;
        impl_->stats.members_run += impl_->config.ensemble_size - sentinels;
        ++impl_->stats.expanded_shots;
        used = impl_->config.ensemble_size;
    }
    Result result;
    result.flipped_obs = ObsType(num_observables);
    result.flipped_obs[0] = 2 * ones > used;
    result.matching_data.frame_flips = result.flipped_obs;
    result.matching_data.probability = double(std::max(ones, used - ones)) / used;
    const auto ns = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now() - start).count());
    result.decode_latency_ns = ns;
    auto& stats = impl_->stats;
    ++stats.shots;
    stats.decode_total_ns += ns;
    stats.decode_min_ns = std::min(stats.decode_min_ns, ns);
    stats.decode_max_ns = std::max(stats.decode_max_ns, ns);
    return result;
}

Harmony::Statistics Harmony::statistics() const {
    auto result = impl_->stats;
    if (!result.shots) result.decode_min_ns = 0;
    return result;
}
void Harmony::print_stats(std::ostream& stream) const {
    const auto s = statistics();
    stream << "Harmony shots=" << s.shots << " expanded=" << s.expanded_shots
           << " members_run=" << s.members_run
           << " mean_members_per_shot=" << (s.shots ? double(s.members_run) / s.shots : 0)
           << '\n';
    stream << "Harmony decode latency (us): timed_shots=" << s.shots << " average="
           << (s.shots ? double(s.decode_total_ns) / s.shots / 1000 : 0)
           << " min=" << double(s.decode_min_ns) / 1000
           << " max=" << double(s.decode_max_ns) / 1000
           << " includes_empty_syndromes=true scope=full_history\n";
}
void Harmony::mpi_accumulate() {
#if defined(ENABLE_MPI)
    auto& s = impl_->stats;
    uint64_t sums[]{s.shots, s.expanded_shots, s.members_run, s.decode_total_ns};
    MPI_Allreduce(MPI_IN_PLACE, sums, 4, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    s.shots = sums[0]; s.expanded_shots = sums[1]; s.members_run = sums[2]; s.decode_total_ns = sums[3];
    MPI_Allreduce(MPI_IN_PLACE, &s.decode_min_ns, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &s.decode_max_ns, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
#endif
}
} // namespace dec
