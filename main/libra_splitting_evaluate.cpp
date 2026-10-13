/*
 *  contributor: OpenAI Codex
 * */

/*
 * SI1000 failure-conditioned likelihood-ratio splitting evaluator.
 *
 * Usage:
 *   libra_splitting_evaluate D WORKERS BASE_MAX_SHOTS SEED OUTPUT
 *       [BASE_THETA=0.8] [LEVELS=8] [BURNIN=1000] [RECORDS=1000] [STRIDE=4]
 */
#include "circuit_generator.h"
#include "decoder/libra.h"
#include "decoder/surface_code.h"
#include "libra_splitting.h"
#if defined(LIBRA_SPLITTING_LOCAL)
#include "libra_local_splitting.h"
#endif

#include <stim.h>
#include <stim/util_top/circuit_to_dem.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;
#if defined(LIBRA_SPLITTING_LOCAL)
using SplittingModel = LibraLocalSplittingModel;
constexpr const char* kProposalMode = "half_global_birth_death_half_nearby_error_swap";
#else
using SplittingModel = LibraSplittingModel;
constexpr const char* kProposalMode = "global_birth_death";
#endif
#if defined(LIBRA_SPLITTING_SYSTEMATIC)
constexpr const char* kResamplingMode = "systematic_over_worker_weight_divided_across_four_reservoir_slots_then_shuffled";
#else
constexpr const char* kResamplingMode = "multinomial_worker_then_uniform_reservoir_slot";
#endif
constexpr size_t kBatch = 1024;
constexpr uint64_t kErrorTarget = 25;
constexpr uint64_t kSamplingSeedDomain = 0xd6e8feb86659fd93ULL;

uint64_t splitmix64_permute(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

uint64_t sampling_seed(uint64_t root_seed, uint32_t phase, uint32_t stream) {
    const uint64_t key = (uint64_t(phase) << 32) | stream;
    return splitmix64_permute(root_seed ^ kSamplingSeedDomain ^ key);
}

#if defined(LIBRA_SPLITTING_FAST)
constexpr const char* kCallbackMode = "pm_first_union_shortcut";
#else
constexpr const char* kCallbackMode = "full_pm_and_libra_proposals";
#endif

struct Options {
    int distance = 0, workers = 0, levels = 8;
    uint64_t base_max_shots = 0, seed = 0, burnin = 1000, records = 1000, stride = 4;
    double base_theta = 0.8;
    std::string output;
};

uint64_t u64(const char* text, const char* label) {
    size_t used = 0;
    std::string s(text);
    uint64_t v = std::stoull(s, &used);
    if (used != s.size()) throw std::invalid_argument(std::string("invalid ") + label);
    return v;
}

Options parse(int argc, char** argv) {
    if (argc < 6) throw std::invalid_argument(
        "usage: libra_splitting_evaluate D WORKERS BASE_MAX_SHOTS SEED OUTPUT "
        "[BASE_THETA=0.8] [LEVELS=8] [BURNIN=1000] [RECORDS=1000] [STRIDE=4]");
    Options o;
    o.distance = static_cast<int>(u64(argv[1], "distance"));
    o.workers = static_cast<int>(u64(argv[2], "workers"));
    o.base_max_shots = u64(argv[3], "base shot cap");
    o.seed = u64(argv[4], "seed");
    o.output = argv[5];
    if (argc >= 7) o.base_theta = std::stod(argv[6]);
    if (argc >= 8) o.levels = static_cast<int>(u64(argv[7], "levels"));
    if (argc >= 9) o.burnin = u64(argv[8], "burn-in");
    if (argc >= 10) o.records = u64(argv[9], "records");
    if (argc >= 11) o.stride = u64(argv[10], "stride");
    if (argc > 11) throw std::invalid_argument("too many arguments");
    if (o.distance < 2 || o.workers < 1 || o.workers > 64 || !o.base_max_shots ||
        !std::isfinite(o.base_theta) || o.base_theta < 0 || o.levels < 1 ||
        !o.records || !o.stride)
        throw std::invalid_argument("invalid distance, workers, theta, levels, records, or stride");
    return o;
}

struct WorkerDecoder {
    std::unique_ptr<dec::Libra> libra;
    std::unique_ptr<dec::PyMatching> pm;
};

struct BaseStats {
    uint64_t shots = 0, pm_errors = 0, libra_errors = 0, union_errors = 0;
    std::vector<std::vector<size_t>> seeds;
};

struct StageStats {
    long double ratio_sum = 0, ratio_sum2 = 0;
    long double weighted_pm_fail = 0, weighted_libra_fail = 0;
    long double weighted_both_fail = 0, weighted_union_fail = 0;
    long double weight_sum = 0;
    uint64_t ratio_count = 0;
    std::array<std::vector<size_t>, 4> reservoir;
    std::array<bool, 4> reservoir_filled{};
    uint64_t endpoint_records = 0, endpoint_both = 0, endpoint_pm_only = 0, endpoint_libra_only = 0;
    uint64_t attempts = 0, physical_accepts = 0, event_accepts = 0;
    uint64_t proposal_libra_calls = 0, proposal_libra_skips = 0, record_libra_calls = 0;
};

dec::Libra::Dimensions build_dimensions(const stim::DetectorErrorModel& dem) {
    dec::Libra::Config cfg;
    cfg.ensemble_size = 0;
    return dec::Libra(dem, cfg).dimensions();
}

LibraSplittingModel::Labels decide(dec::PyMatching& pm, dec::Libra& libra,
                                   dec::SyndromeRef syndrome, dec::ObsRef truth) {
    const auto p = pm.decode(syndrome, truth);
    const auto l = libra.decode(syndrome, truth);
    return {p.flipped_obs[0] != truth[0], l.flipped_obs[0] != truth[0]};
}

LibraSplittingModel::Labels decide_union(dec::PyMatching& pm, dec::Libra& libra,
        dec::SyndromeRef syndrome, dec::ObsRef truth, StageStats& stats) {
    const auto p = pm.decode(syndrome, truth);
    if (p.flipped_obs[0] != truth[0]) {
        ++stats.proposal_libra_skips;
        return {true, false};
    }
    ++stats.proposal_libra_calls;
    const auto l = libra.decode(syndrome, truth);
    return {false, l.flipped_obs[0] != truth[0]};
}

std::vector<double> theta_schedule(const Options& o) {
    std::vector<double> values(o.levels + 1);
    for (int i = 0; i <= o.levels; ++i)
        values[i] = o.base_theta * (o.levels - i) / o.levels;
    return values;
}

stim::DetectorErrorModel make_tilted_dem(const stim::DetectorErrorModel& dem, double theta) {
    stim::DetectorErrorModel tilted;
    dem.iter_flatten_error_instructions([&](const auto& inst) {
        const double p = inst.arg_data[0];
        const double odds = p == 0 ? 0 : p / (1 - p) * std::exp(theta);
        const double p_tilted = odds / (1 + odds);
        tilted.append_error_instruction(p_tilted, inst.target_data, "");
    });
    if (tilted.count_errors() != dem.count_errors() ||
        tilted.count_detectors() != dem.count_detectors() ||
        tilted.count_observables() != dem.count_observables())
        throw std::runtime_error("tilted DEM changed flattened channel/detector dimensions");
    return tilted;
}

void run(const Options& o, const stim::DetectorErrorModel& dem,
         const dec::Libra::Dimensions& dims, std::ostream& out) {
    SplittingModel model(dem);
    const stim::DetectorErrorModel tilted_dem = make_tilted_dem(dem, o.base_theta);
    const auto thetas = theta_schedule(o);
    std::vector<double> log_z(thetas.size());
    for (size_t i = 0; i < thetas.size(); ++i)
        log_z[i] = model.log_normaliser(thetas[i]);

    std::vector<WorkerDecoder> decoders(o.workers);
    std::vector<BaseStats> base(o.workers);
    std::atomic<uint64_t> base_ticket{0}, base_shots{0}, base_pm_errors{0},
        base_libra_errors{0}, base_union_errors{0};
    std::atomic<int> base_workers{o.workers};
    std::atomic<bool> failed{false};
    std::mutex error_mutex;
    std::string error_message;
    std::vector<std::thread> threads;
    const auto started = Clock::now();

    // Baseline tilted independent Bernoulli sampling uses Stim's vectorized
    // sampler on a DEM whose instruction probabilities are odds-tilted.
    for (int w = 0; w < o.workers; ++w) {
        threads.emplace_back([&, w] {
            try {
                dec::Libra::Config lc;
                lc.seed = o.seed;
                decoders[w].libra = std::make_unique<dec::Libra>(dem, lc);
                decoders[w].pm = std::make_unique<dec::PyMatching>(dem, false);
                decoders[w].libra->prepare_ensemble();

                const uint64_t worker_seed = sampling_seed(o.seed, 0, static_cast<uint32_t>(w));
                stim::DemSampler<64> sampler(
                    tilted_dem, std::mt19937_64(worker_seed), kBatch);
                for (;;) {
                    if (base_pm_errors.load() >= kErrorTarget &&
                        base_libra_errors.load() >= kErrorTarget) break;
                    const uint64_t start = base_ticket.fetch_add(kBatch);
                    if (start >= o.base_max_shots) break;
                    const size_t count = std::min<uint64_t>(kBatch, o.base_max_shots - start);
                    sampler.resample(false);
                    const auto dets = sampler.det_buffer.transposed();
                    const auto obs = sampler.obs_buffer.transposed();
                    const auto errors = sampler.err_buffer.transposed();
                    for (size_t s = 0; s < count; ++s) {
                        if (base_pm_errors.load() >= kErrorTarget &&
                            base_libra_errors.load() >= kErrorTarget) break;
                        const auto p = decoders[w].pm->decode(dets[s], obs[s]);
                        const auto l = decoders[w].libra->decode(dets[s], obs[s]);
                        const bool pm_error = p.flipped_obs[0] != obs[s][0];
                        const bool libra_error = l.flipped_obs[0] != obs[s][0];
                        const bool union_error = pm_error || libra_error;
                        ++base[w].shots;
                        ++base_shots;
                        if (pm_error) { ++base[w].pm_errors; ++base_pm_errors; }
                        if (libra_error) { ++base[w].libra_errors; ++base_libra_errors; }
                        if (union_error) {
                            ++base[w].union_errors;
                            ++base_union_errors;
                            std::vector<size_t> active;
                            for (size_t e = 0; e < model.probabilities.size(); ++e)
                                if (errors[s][e]) active.push_back(e);
                            base[w].seeds.push_back(std::move(active));
                        }
                    }
                }
            } catch (const std::exception& e) {
                failed.store(true);
                std::lock_guard<std::mutex> lock(error_mutex);
                if (error_message.empty()) error_message = e.what();
            }
            base_workers.fetch_sub(1);
        });
    }
    auto last_progress = started;
    while (base_workers.load() > 0) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        const auto now = Clock::now();
        if (std::chrono::duration_cast<std::chrono::seconds>(now-last_progress).count() >= 10) {
            std::cerr << "progress phase=base shots=" << base_shots.load()
                      << " cap=" << o.base_max_shots
                      << " pm_errors=" << base_pm_errors.load()
                      << " libra_errors=" << base_libra_errors.load()
                      << " union_errors=" << base_union_errors.load()
                      << " elapsed_s=" << std::chrono::duration<double>(now-started).count() << "\n";
            last_progress = now;
        }
    }
    for (auto& t : threads) t.join();
    if (failed) throw std::runtime_error(error_message);
    if (base_union_errors.load() == 0)
        throw std::runtime_error("base tilted Monte Carlo found no union failures; cannot initialize splitting");

    std::vector<std::vector<size_t>> initial_states;
    for (auto& worker : base)
        for (auto& state : worker.seeds)
            initial_states.push_back(std::move(state));
    if (initial_states.empty()) throw std::runtime_error("no captured union-failure states");

    const double base_union_rate = double(base_union_errors.load()) / base_shots.load();
    std::vector<std::vector<size_t>> current_seeds(o.workers);
    for (int w = 0; w < o.workers; ++w)
        current_seeds[w] = initial_states[(size_t(w) * initial_states.size()) / o.workers];
    std::vector<std::vector<StageStats>> stages(
        o.workers, std::vector<StageStats>(o.levels));
    std::vector<double> level_means(o.levels), stage_ess(o.levels);
    std::vector<size_t> stage_unique(o.levels), stage_reservoir_unique(o.levels);
    std::vector<uint64_t> stage_attempts(o.levels), stage_physical(o.levels), stage_event(o.levels);
    std::vector<uint64_t> stage_proposal_libra_calls(o.levels), stage_proposal_libra_skips(o.levels);
    std::vector<uint64_t> stage_record_libra_calls(o.levels);
    std::atomic<uint64_t> progress_records{0};
    std::atomic<int> live_chains{0};
    long double likelihood_ratio = 1;

    auto run_parallel_phase = [&](const std::string& phase, auto&& body) {
        threads.clear();
        live_chains.store(o.workers);
        for (int w = 0; w < o.workers; ++w) {
            threads.emplace_back([&, w] {
                try {
                    body(w);
                } catch (const std::exception& e) {
                    failed.store(true);
                    std::lock_guard<std::mutex> lock(error_mutex);
                    if (error_message.empty()) error_message = e.what();
                }
                live_chains.fetch_sub(1);
            });
        }
        auto progress_at = Clock::now();
        while (live_chains.load() > 0) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            const auto now = Clock::now();
            if (std::chrono::duration_cast<std::chrono::seconds>(now-progress_at).count() >= 10) {
                std::cerr << "progress phase=" << phase << " records=" << progress_records.load()
                          << " workers_active=" << live_chains.load()
                          << " elapsed_s=" << std::chrono::duration<double>(now-started).count() << "\n";
                progress_at = now;
            }
        }
        for (auto& t : threads) t.join();
        if (failed) throw std::runtime_error(error_message);
    };

    for (int stage = 0; stage < o.levels; ++stage) {
        const double from = thetas[stage], to = thetas[stage + 1];
        const double log_constant = log_z[stage] - log_z[stage + 1];
        run_parallel_phase("stage_" + std::to_string(stage), [&](int w) {
            auto& pm = *decoders[w].pm;
            auto& libra = *decoders[w].libra;
            std::mt19937_64 rng(sampling_seed(o.seed, static_cast<uint32_t>(stage + 1),
                                               static_cast<uint32_t>(w)));
            std::uniform_real_distribution<double> unit(0, 1);
            auto full_callback = [&](dec::SyndromeRef syn, dec::ObsRef obs) {
                return decide(pm, libra, syn, obs);
            };
            auto& stats = stages[w][stage];
            auto state = model.make_state(current_seeds[w], full_callback);
#if defined(LIBRA_SPLITTING_FAST)
            auto callback = [&](dec::SyndromeRef syn, dec::ObsRef obs) {
                return decide_union(pm, libra, syn, obs, stats);
            };
#else
            auto& callback = full_callback;
#endif
            for (uint64_t j = 0; j < o.burnin; ++j)
                model.step(state, from, rng, callback);
            for (uint64_t r = 0; r < o.records; ++r) {
                for (uint64_t s = 0; s < o.stride; ++s)
                    model.step(state, from, rng, callback);
#if defined(LIBRA_SPLITTING_FAST)
                state.labels = full_callback(state.syndrome, state.truth);
                ++stats.record_libra_calls;
#endif
                const double ratio = std::exp(log_constant +
                    (to - from) * state.active.size());
                if (!std::isfinite(ratio) || !(ratio > 0))
                    throw std::runtime_error("non-finite stage likelihood ratio");
                stats.ratio_sum += ratio;
                stats.ratio_sum2 += ratio * ratio;
                stats.weighted_pm_fail += ratio * state.labels.pymatching;
                stats.weighted_libra_fail += ratio * state.labels.libra;
                stats.weighted_both_fail += ratio *
                    (state.labels.pymatching && state.labels.libra);
                stats.weighted_union_fail += ratio * state.labels.any();
                stats.weight_sum += ratio;
                ++stats.ratio_count;
                for (size_t slot = 0; slot < stats.reservoir.size(); ++slot) {
                    if (!stats.reservoir_filled[slot] ||
                        unit(rng) < ratio / stats.weight_sum) {
                        stats.reservoir[slot] = state.active;
                        stats.reservoir_filled[slot] = true;
                    }
                }
                ++progress_records;
            }
            stats.attempts = state.attempts;
            stats.physical_accepts = state.physical_accepts;
            stats.event_accepts = state.event_accepts;
        });

        long double sum_weight = 0, sum_weight2 = 0, sum_ratio = 0;
        uint64_t ratio_count = 0;
        std::set<std::vector<size_t>> unique_reservoir_states;
        std::vector<double> parent_weights(o.workers);
        for (int w = 0; w < o.workers; ++w) {
            const auto& stats = stages[w][stage];
            sum_weight += stats.weight_sum;
            sum_weight2 += stats.weight_sum * stats.weight_sum;
            sum_ratio += stats.ratio_sum;
            ratio_count += stats.ratio_count;
            parent_weights[w] = static_cast<double>(stats.weight_sum);
            for (size_t slot = 0; slot < stats.reservoir.size(); ++slot)
                if (stats.reservoir_filled[slot])
                    unique_reservoir_states.insert(stats.reservoir[slot]);
            stage_attempts[stage] += stats.attempts;
            stage_physical[stage] += stats.physical_accepts;
            stage_event[stage] += stats.event_accepts;
            stage_proposal_libra_calls[stage] += stats.proposal_libra_calls;
            stage_proposal_libra_skips[stage] += stats.proposal_libra_skips;
            stage_record_libra_calls[stage] += stats.record_libra_calls;
        }
        if (!(sum_weight > 0) || !(sum_weight2 > 0))
            throw std::runtime_error("stage population has zero importance weight");
        level_means[stage] = static_cast<double>(sum_ratio / ratio_count);
        stage_ess[stage] = static_cast<double>(sum_weight * sum_weight / sum_weight2);
        stage_reservoir_unique[stage] = unique_reservoir_states.size();

        std::mt19937_64 resample_rng(sampling_seed(o.seed, static_cast<uint32_t>(stage + 1),
                                                   UINT32_MAX));
        std::set<std::vector<size_t>> resampled_unique;
        std::vector<std::vector<size_t>> next_seeds(o.workers);
#if defined(LIBRA_SPLITTING_SYSTEMATIC)
        std::vector<const std::vector<size_t>*> candidates;
        std::vector<double> cumulative_weights;
        candidates.reserve(size_t(o.workers) * 4);
        cumulative_weights.reserve(size_t(o.workers) * 4);
        double total_candidate_weight = 0;
        for (int parent = 0; parent < o.workers; ++parent) {
            const auto& parent_stats = stages[parent][stage];
            for (size_t slot = 0; slot < parent_stats.reservoir.size(); ++slot) {
                if (!parent_stats.reservoir_filled[slot])
                    throw std::runtime_error("weighted reservoir unexpectedly empty");
                total_candidate_weight += static_cast<double>(parent_stats.weight_sum / 4);
                candidates.push_back(&parent_stats.reservoir[slot]);
                cumulative_weights.push_back(total_candidate_weight);
            }
        }
        const double spacing = total_candidate_weight / o.workers;
        std::uniform_real_distribution<double> choose_offset(0, spacing);
        const double offset = choose_offset(resample_rng);
        std::vector<std::vector<size_t>> selected;
        selected.reserve(o.workers);
        size_t candidate = 0;
        for (int j = 0; j < o.workers; ++j) {
            const double position = offset + j * spacing;
            while (candidate + 1 < cumulative_weights.size() &&
                   position >= cumulative_weights[candidate])
                ++candidate;
            selected.push_back(*candidates[candidate]);
        }
        std::shuffle(selected.begin(), selected.end(), resample_rng);
        for (int w = 0; w < o.workers; ++w) {
            next_seeds[w] = std::move(selected[w]);
            resampled_unique.insert(next_seeds[w]);
        }
#else
        // Multinomial resampling: select a worker by aggregate importance,
        // then one of its four independent weighted-reservoir representatives.
        std::discrete_distribution<int> choose_worker(parent_weights.begin(), parent_weights.end());
        std::uniform_int_distribution<size_t> choose_reservoir(0, 3);
        for (int w = 0; w < o.workers; ++w) {
            const int parent = choose_worker(resample_rng);
            const auto& parent_stats = stages[parent][stage];
            if (!parent_stats.reservoir_filled[0] ||
                !parent_stats.reservoir_filled[1] ||
                !parent_stats.reservoir_filled[2] ||
                !parent_stats.reservoir_filled[3])
                throw std::runtime_error("weighted reservoir unexpectedly empty");
            const auto& seed = parent_stats.reservoir[choose_reservoir(resample_rng)];
            next_seeds[w] = seed;
            resampled_unique.insert(seed);
        }
#endif
        stage_unique[stage] = resampled_unique.size();
        current_seeds = std::move(next_seeds);
        likelihood_ratio *= level_means[stage];
    }

    std::vector<StageStats> endpoints(o.workers);
    run_parallel_phase("endpoint_theta_0", [&](int w) {
        auto& pm = *decoders[w].pm;
        auto& libra = *decoders[w].libra;
        std::mt19937_64 rng(sampling_seed(o.seed, static_cast<uint32_t>(o.levels + 1),
                                           static_cast<uint32_t>(w)));
        auto full_callback = [&](dec::SyndromeRef syn, dec::ObsRef obs) {
            return decide(pm, libra, syn, obs);
        };
        auto& stats = endpoints[w];
        auto state = model.make_state(current_seeds[w], full_callback);
#if defined(LIBRA_SPLITTING_FAST)
        auto callback = [&](dec::SyndromeRef syn, dec::ObsRef obs) {
            return decide_union(pm, libra, syn, obs, stats);
        };
#else
        auto& callback = full_callback;
#endif
        for (uint64_t j = 0; j < o.burnin; ++j)
            model.step(state, 0, rng, callback);
        for (uint64_t r = 0; r < o.records; ++r) {
            for (uint64_t s = 0; s < o.stride; ++s)
                model.step(state, 0, rng, callback);
#if defined(LIBRA_SPLITTING_FAST)
            state.labels = full_callback(state.syndrome, state.truth);
            ++stats.record_libra_calls;
#endif
            ++stats.endpoint_records;
            stats.endpoint_both += state.labels.pymatching && state.labels.libra;
            stats.endpoint_pm_only += state.labels.pymatching && !state.labels.libra;
            stats.endpoint_libra_only += !state.labels.pymatching && state.labels.libra;
            ++progress_records;
        }
        stats.attempts = state.attempts;
        stats.physical_accepts = state.physical_accepts;
        stats.event_accepts = state.event_accepts;
    });
    long double ratio_before_last = 1;
    for (int stage = 0; stage + 1 < o.levels; ++stage)
        ratio_before_last *= level_means[stage];
    long double last_pm_weighted = 0, last_libra_weighted = 0;
    long double last_both_weighted = 0, last_union_weighted = 0;
    uint64_t last_bridge_records = 0;
    const int last_stage = o.levels - 1;
    for (int w = 0; w < o.workers; ++w) {
        const auto& s = stages[w][last_stage];
        last_pm_weighted += s.weighted_pm_fail;
        last_libra_weighted += s.weighted_libra_fail;
        last_both_weighted += s.weighted_both_fail;
        last_union_weighted += s.weighted_union_fail;
        last_bridge_records += s.ratio_count;
    }
    const double before_last_rate = base_union_rate * static_cast<double>(ratio_before_last);
    const double union_rate_theta0 = before_last_rate *
        static_cast<double>(last_union_weighted / last_bridge_records);
    const double pm_ler = before_last_rate *
        static_cast<double>(last_pm_weighted / last_bridge_records);
    const double libra_ler = before_last_rate *
        static_cast<double>(last_libra_weighted / last_bridge_records);
    uint64_t endpoint_n = 0, endpoint_both = 0, endpoint_pm_only = 0, endpoint_libra_only = 0;
    for (const auto& stats : endpoints) {
        endpoint_n += stats.endpoint_records;
        endpoint_both += stats.endpoint_both;
        endpoint_pm_only += stats.endpoint_pm_only;
        endpoint_libra_only += stats.endpoint_libra_only;
    }
    const double pm_fraction = double(endpoint_both + endpoint_pm_only) / endpoint_n;
    const double libra_fraction = double(endpoint_both + endpoint_libra_only) / endpoint_n;
    const double pm_ler_endpoint = union_rate_theta0 * pm_fraction;
    const double libra_ler_endpoint = union_rate_theta0 * libra_fraction;
    uint64_t endpoint_proposal_calls = 0, endpoint_proposal_skips = 0,
        endpoint_record_calls = 0;
    for (const auto& stats : endpoints) {
        endpoint_proposal_calls += stats.proposal_libra_calls;
        endpoint_proposal_skips += stats.proposal_libra_skips;
        endpoint_record_calls += stats.record_libra_calls;
    }

    out << std::setprecision(14)
        << "method=failure_conditioned_likelihood_ratio_splitting\n"
        << "algorithm=sequential_importance_resampling_with_conditional_MH\n"
        << "proposal_kernel=" << kProposalMode << "\n"
        << "proposal_callback_mode=" << kCallbackMode << "\n"
        << "random_stream_scheme=bijective_splitmix64_per_phase_and_worker_with_separate_sampling_domain\n"
        << "resampling_scheme=" << kResamplingMode << "\n"
        << "weighted_reservoirs_per_worker=4\n"
        << "distance=" << o.distance << "\nrounds=" << o.distance
        << "\np=0.001\nmemory_basis=Z\ninclude_opposite_basis_detectors=true\n"
        << "dem_error_instructions=" << model.probabilities.size() << "\n"
        << "detectors=" << dims.detectors << "\nhyperedges=" << dims.hyperedges
        << "\ngraph_edges=" << dims.graph_edges << "\ncorrelation_rules=" << dims.correlation_rules << "\n"
        << "workers=" << o.workers << "\nseed=" << o.seed << "\n"
        << "base_theta=" << o.base_theta << "\ntheta_schedule=";
    for (double theta : thetas) out << theta << ",";
    out << "\nlikelihood_ratio_formula=exp((theta_next-theta)*K + logM(theta)-logM(theta_next))\n"
        << "M_theta=product_i(1-p_i+p_i*exp(theta))\n"
        << "base_max_shots=" << o.base_max_shots << "\n"
        << "base_shots=" << base_shots.load() << "\n"
        << "base_pm_failures=" << base_pm_errors.load() << "\n"
        << "base_libra_failures=" << base_libra_errors.load() << "\n"
        << "base_union_failures=" << base_union_errors.load() << "\n"
        << "base_pm_rate_tilted=" << double(base_pm_errors.load()) / base_shots.load() << "\n"
        << "base_libra_rate_tilted=" << double(base_libra_errors.load()) / base_shots.load() << "\n"
        << "base_union_rate_tilted=" << base_union_rate << "\n"
        << "baseline_sampler=Stim_DemSampler_on_odds_tilted_DEM\n"
        << "mc_error_target_each_decoder=" << kErrorTarget << "\n"
        << "base_stopping_rule=stop_when_both_decoders_reach_target_or_base_shot_cap\n"
        << "burnin=" << o.burnin << "\nrecords_per_chain_per_stage=" << o.records
        << "\nstride=" << o.stride << "\n";
    for (int i = 0; i < o.levels; ++i) {
        long double stage_pm = 0, stage_libra = 0, stage_both = 0, stage_union = 0;
        for (int w = 0; w < o.workers; ++w) {
            stage_pm += stages[w][i].weighted_pm_fail;
            stage_libra += stages[w][i].weighted_libra_fail;
            stage_both += stages[w][i].weighted_both_fail;
            stage_union += stages[w][i].weighted_union_fail;
        }
        const uint64_t stage_n = uint64_t(o.workers) * o.records;
        out << "ratio_stage_" << i << "_mean=" << level_means[i] << "\n"
            << "ratio_stage_" << i << "_mean_R_times_pm_fail=" << static_cast<double>(stage_pm / stage_n) << "\n"
            << "ratio_stage_" << i << "_mean_R_times_libra_fail=" << static_cast<double>(stage_libra / stage_n) << "\n"
            << "ratio_stage_" << i << "_mean_R_times_both_fail=" << static_cast<double>(stage_both / stage_n) << "\n"
            << "ratio_stage_" << i << "_mean_R_times_union_fail=" << static_cast<double>(stage_union / stage_n) << "\n"
            << "ratio_stage_" << i << "_worker_weight_ess=" << stage_ess[i] << "\n"
            << "ratio_stage_" << i << "_unique_reservoir_configurations="
            << stage_reservoir_unique[i] << "\n"
            << "ratio_stage_" << i << "_unique_resampled_configurations=" << stage_unique[i] << "\n"
            << "ratio_stage_" << i << "_attempts=" << stage_attempts[i] << "\n"
            << "ratio_stage_" << i << "_physical_accepts=" << stage_physical[i] << "\n"
            << "ratio_stage_" << i << "_event_accepts=" << stage_event[i] << "\n"
            << "ratio_stage_" << i << "_proposal_libra_calls=" << stage_proposal_libra_calls[i] << "\n"
            << "ratio_stage_" << i << "_proposal_libra_skips=" << stage_proposal_libra_skips[i] << "\n"
            << "ratio_stage_" << i << "_record_libra_calls=" << stage_record_libra_calls[i] << "\n";
    }
    out << "likelihood_ratio_product=" << static_cast<double>(likelihood_ratio) << "\n"
        << "ratio_product_before_last_bridge=" << static_cast<double>(ratio_before_last) << "\n"
        << "last_bridge_mean_R=" << static_cast<double>(last_union_weighted / last_bridge_records) << "\n"
        << "last_bridge_mean_R_pm_fail=" << static_cast<double>(last_pm_weighted / last_bridge_records) << "\n"
        << "last_bridge_mean_R_libra_fail=" << static_cast<double>(last_libra_weighted / last_bridge_records) << "\n"
        << "last_bridge_mean_R_both_fail=" << static_cast<double>(last_both_weighted / last_bridge_records) << "\n"
        << "union_failure_rate_theta0_estimate=" << union_rate_theta0 << "\n"
        << "endpoint_records=" << endpoint_n << "\n"
        << "endpoint_both_fail=" << endpoint_both << "\n"
        << "endpoint_pymatching_only_fail=" << endpoint_pm_only << "\n"
        << "endpoint_libra_only_fail=" << endpoint_libra_only << "\n"
        << "endpoint_proposal_libra_calls=" << endpoint_proposal_calls << "\n"
        << "endpoint_proposal_libra_skips=" << endpoint_proposal_skips << "\n"
        << "endpoint_record_libra_calls=" << endpoint_record_calls << "\n"
        << "endpoint_pymatching_fraction_given_union=" << static_cast<double>(pm_fraction) << "\n"
        << "endpoint_libra_fraction_given_union=" << static_cast<double>(libra_fraction) << "\n"
        << "endpoint_pymatching_ler_mixing_check=" << pm_ler_endpoint << "\n"
        << "endpoint_libra_ler_mixing_check=" << libra_ler_endpoint << "\n"
        << "pymatching_ler_estimate=" << pm_ler << "\n"
        << "libra_ler_estimate=" << libra_ler << "\n"
        << "ler_estimator=pre_last_bridge_union_rate_times_mean_R_times_decoder_failure_indicator\n"
        << "uncertainty_note=finite-particle resampling and correlated chains; no exact confidence interval claimed\n"
        << "resampling_note=each worker retains four independent ratio-weighted reservoirs; reservoir slot weight is worker importance sum divided by four\n"
        << "endpoint_note=theta0 raw-label estimates are mixing cross-checks only; zero endpoint class counts do not establish zero failure probability\n"
        << "base_pool_note=base failure pool is stopped using per-decoder error targets and may have stopping-composition bias\n"
        << "sampling_note=baseline errors are drawn from the odds-tilted DEM; decoders always use original p=1e-3 DEM\n"
        << "elapsed_wall_seconds=" << std::chrono::duration<double>(Clock::now()-started).count() << "\n";

    for (int w = 0; w < o.workers; ++w) {
        out << "worker_" << w << "_endpoint_pm_fraction="
            << double(endpoints[w].endpoint_both + endpoints[w].endpoint_pm_only) /
                endpoints[w].endpoint_records << "\n"
            << "worker_" << w << "_endpoint_libra_fraction="
            << double(endpoints[w].endpoint_both + endpoints[w].endpoint_libra_only) /
                endpoints[w].endpoint_records << "\n"
            << "worker_" << w << "_endpoint_both_fail=" << endpoints[w].endpoint_both << "\n"
            << "worker_" << w << "_endpoint_pm_only_fail=" << endpoints[w].endpoint_pm_only << "\n"
            << "worker_" << w << "_endpoint_libra_only_fail=" << endpoints[w].endpoint_libra_only << "\n"
            << "worker_" << w << "_endpoint_attempts=" << endpoints[w].attempts << "\n"
            << "worker_" << w << "_endpoint_physical_accepts=" << endpoints[w].physical_accepts << "\n"
            << "worker_" << w << "_endpoint_event_accepts=" << endpoints[w].event_accepts << "\n";
        out << "worker_" << w << "_endpoint_proposal_libra_calls="
            << endpoints[w].proposal_libra_calls << "\n"
            << "worker_" << w << "_endpoint_proposal_libra_skips="
            << endpoints[w].proposal_libra_skips << "\n"
            << "worker_" << w << "_endpoint_record_libra_calls="
            << endpoints[w].record_libra_calls << "\n";
        for (int stage = 0; stage < o.levels; ++stage)
            out << "worker_" << w << "_ratio_mean_" << stage << "="
                << static_cast<double>(stages[w][stage].ratio_sum / stages[w][stage].ratio_count) << "\n"
                << "worker_" << w << "_ratio_weight_sum_" << stage << "="
                << static_cast<double>(stages[w][stage].weight_sum) << "\n"
                << "worker_" << w << "_attempts_" << stage << "=" << stages[w][stage].attempts << "\n"
                << "worker_" << w << "_physical_accepts_" << stage << "="
                << stages[w][stage].physical_accepts << "\n"
                << "worker_" << w << "_event_accepts_" << stage << "="
                << stages[w][stage].event_accepts << "\n"
                << "worker_" << w << "_proposal_libra_calls_" << stage << "="
                << stages[w][stage].proposal_libra_calls << "\n"
                << "worker_" << w << "_proposal_libra_skips_" << stage << "="
                << stages[w][stage].proposal_libra_skips << "\n"
                << "worker_" << w << "_record_libra_calls_" << stage << "="
                << stages[w][stage].record_libra_calls << "\n";
        out << "worker_" << w << "_last_bridge_sum_R_pm_fail="
            << static_cast<double>(stages[w][last_stage].weighted_pm_fail) << "\n"
            << "worker_" << w << "_last_bridge_sum_R_libra_fail="
            << static_cast<double>(stages[w][last_stage].weighted_libra_fail) << "\n";
    }
}
}  // namespace

int main(int argc, char** argv) {
    try {
        const Options o = parse(argc, argv);
        const auto circuit = sc_si1000(o.distance, o.distance, 1e-3, false, true);
        const auto dem = stim::circuit_to_dem(circuit, {.decompose_errors = true});
        const auto dims = build_dimensions(dem);
        std::ofstream out(o.output);
        if (!out) throw std::runtime_error("cannot create output");
        run(o, dem, dims, out);
    } catch (const std::exception& e) {
        std::cerr << "libra_splitting_evaluate: " << e.what() << "\n";
        return 2;
    }
    return 0;
}
