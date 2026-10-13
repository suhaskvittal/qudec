/*
 *  contributor: OpenAI Codex
 * */

/*
 * Standalone Libra / PyMatching SI1000 evaluator.
 *
 * Usage:
 *   libra_evaluate mc DISTANCE WORKERS MAX_SHOTS SEED OUTPUT [ERROR_TARGET] [ENSEMBLE]
 *   libra_evaluate latency DISTANCE WORKERS SHOTS SEED OUTPUT [ENSEMBLE]
 *   libra_evaluate rare DISTANCE WORKERS SHOTS SEED OUTPUT LEVEL [ENSEMBLE]
 */
#include "circuit_generator.h"
#include "decoder/libra.h"
#include "decoder/surface_code.h"
#include "libra_rare_event.h"
#include "sampler.h"

#include <stim.h>
#include <stim/util_top/circuit_to_dem.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <sstream>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;
constexpr size_t kBatchSize = 1024;

struct Counts {
    uint64_t shots = 0, pymatching_errors = 0, libra_errors = 0;
    uint64_t both_errors = 0, pymatching_only = 0, libra_only = 0;
    long double pymatching_ns = 0, libra_ns = 0;
    uint64_t ensemble_shots = 0, fallback_shots = 0, improving_cycles = 0;
    uint64_t complementary_ns = 0, ensemble_ns = 0, synthesis_ns = 0, degeneracy_ns = 0;
    uint64_t fallback_ns = 0;
};

struct Options {
    std::string mode, output;
    int distance = 0, workers = 0, level = -1;
    uint64_t shots = 0, seed = 0, error_target = 25;
    size_t ensemble_size = 60;
    size_t max_count = 256;
    double theta = 0;
    uint64_t progress_seconds = 10;
};

uint64_t parse_u64(const char* s, const char* what) {
    size_t used = 0;
    const std::string text(s);
    const uint64_t value = std::stoull(text, &used);
    if (used != text.size()) throw std::invalid_argument(std::string("invalid ") + what);
    return value;
}

double parse_double(const char* s, const char* what) {
    size_t used = 0;
    const std::string text(s);
    const double value = std::stod(text, &used);
    if (used != text.size() || !std::isfinite(value))
        throw std::invalid_argument(std::string("invalid ") + what);
    return value;
}

Options parse_args(int argc, char** argv) {
    if (argc < 7) throw std::invalid_argument(
        "usage: libra_evaluate mc|latency DISTANCE WORKERS SHOTS SEED OUTPUT "
        "[--max-errors N] [--level K] [--ensemble-size N] [--progress-seconds N]");
    Options o;
    o.mode = argv[1];
    o.distance = static_cast<int>(parse_u64(argv[2], "distance"));
    o.workers = static_cast<int>(parse_u64(argv[3], "worker count"));
    o.shots = parse_u64(argv[4], "shot count");
    o.seed = parse_u64(argv[5], "seed");
    o.output = argv[6];
    if (o.mode != "mc" && o.mode != "latency" && o.mode != "rare" && o.mode != "rare_mix")
        throw std::invalid_argument("mode must be mc, latency, rare, or rare_mix");
    bool target_seen = false, level_seen = false;
    for (int i = 7; i < argc; ++i) {
        const std::string arg(argv[i]);
        auto value = [&](const std::string& name) -> std::string {
            const std::string prefix = name + "=";
            if (arg.rfind(prefix, 0) == 0) return arg.substr(prefix.size());
            if (arg == name && i + 1 < argc) return argv[++i];
            return {};
        };
        if (auto v = value("--max-errors"); !v.empty())
            o.error_target = parse_u64(v.c_str(), "maximum errors");
        else if (auto v = value("--ensemble-size"); !v.empty())
            o.ensemble_size = parse_u64(v.c_str(), "ensemble size");
        else if (auto v = value("--progress-seconds"); !v.empty())
            o.progress_seconds = parse_u64(v.c_str(), "progress interval");
        else if (auto v = value("--max-count"); !v.empty())
            o.max_count = parse_u64(v.c_str(), "maximum error count");
        else if (auto v = value("--theta"); !v.empty())
            o.theta = parse_double(v.c_str(), "proposal tilt");
        else if (auto v = value("--level"); !v.empty()) {
            o.level = static_cast<int>(parse_u64(v.c_str(), "error level"));
            level_seen = true;
        } else if (arg.rfind("--", 0) == 0)
            throw std::invalid_argument("unknown option: " + arg);
        else if (o.mode == "mc" && !target_seen) {
            o.error_target = parse_u64(argv[i], "error target");
            target_seen = true;
        } else if (o.mode == "rare" && !level_seen) {
            o.level = static_cast<int>(parse_u64(argv[i], "error level"));
            level_seen = true;
        } else
            o.ensemble_size = parse_u64(argv[i], "ensemble size");
    }
    if (o.mode == "rare" && o.level < 0)
        throw std::invalid_argument("rare mode requires --level K");
    if (o.distance < 2 || o.workers < 1 || o.workers > 64 || !o.shots ||
        (o.mode != "latency" && o.error_target == 0))
        throw std::invalid_argument("distance >= 2, workers in [1,64], shots > 0 required");
    return o;
}

stim::Circuit make_circuit(const Options& o) {
    return sc_si1000(o.distance, o.distance, 1e-3, false, true);
}

dec::Libra make_libra(const stim::DetectorErrorModel& dem, uint64_t seed,
                      size_t ensemble_size, bool profile) {
    dec::Libra::Config config;
    config.seed = seed;
    config.ensemble_size = ensemble_size;
    config.profile = profile;
    return dec::Libra(dem, config);
}

void write_metadata(std::ostream& out, const Options& o, const stim::DetectorErrorModel& dem,
                    const dec::Libra::Dimensions& dims) {
    out << "mode=" << o.mode << "\n"
        << "distance=" << o.distance << "\n"
        << "rounds=" << o.distance << "\n"
        << "physical_error_rate=0.001\n"
        << "memory_basis=Z\n"
        << "include_opposite_basis_detectors=true\n"
        << "workers=" << o.workers << "\n"
        << "requested_shots=" << o.shots << "\n"
        << "seed=" << o.seed << "\n"
        << "detectors=" << dims.detectors << "\n"
        << "dem_error_instructions=" << dem.count_errors() << "\n"
        << "libra_hyperedges=" << dims.hyperedges << "\n"
        << "libra_graph_edges=" << dims.graph_edges << "\n"
        << "libra_correlation_rules=" << dims.correlation_rules << "\n"
        << "libra_ensemble_size=" << o.ensemble_size << "\n"
        << (o.mode == "rare_mix" ? "max_count=" + std::to_string(o.max_count) + "\ntheta=" + std::to_string(o.theta) + "\n" : "")
        << (o.mode != "latency" ? "max_errors=" + std::to_string(o.error_target) + "\n" : "")
        << (o.mode == "rare" ? "level=" + std::to_string(o.level) + "\n" : "");
}

void add_stats(Counts& c, const dec::Libra::Statistics& s) {
    c.ensemble_shots += s.ensemble_shots;
    c.fallback_shots += s.fallback_shots;
    c.improving_cycles += s.improving_cycles;
    c.complementary_ns += s.complementary_ns;
    c.fallback_ns += s.fallback_ns;
    c.ensemble_ns += s.ensemble_ns;
    c.synthesis_ns += s.synthesis_ns;
    c.degeneracy_ns += s.degeneracy_ns;
}

void run_monte_carlo(const Options& o, const stim::DetectorErrorModel& dem,
                     const dec::Libra::Dimensions& dims, std::ostream& out, bool latency) {
    std::atomic<uint64_t> ticket{0};
    std::atomic<uint64_t> total_pm_errors{0}, total_libra_errors{0};
    std::atomic<uint64_t> processed{0};
    std::atomic<int> live_workers{o.workers};
    std::atomic<bool> failed{false};
    std::mutex error_mutex;
    std::string error_message;
    std::vector<Counts> results(o.workers);
    std::vector<std::thread> threads;
    const auto run_started = Clock::now();
    for (int worker = 0; worker < o.workers; ++worker) {
        threads.emplace_back([&, worker] {
            try {
                const uint64_t worker_seed = o.seed + 0x9e3779b97f4a7c15ULL * uint64_t(worker + 1);
                auto libra = make_libra(dem, o.seed, o.ensemble_size, latency);
                dec::PyMatching pm(dem, false);
                libra.prepare_ensemble();
                std::mt19937_64 rng(worker_seed);
                stim::DemSampler<64> sampler(dem, std::move(rng), kBatchSize);
                Counts& c = results[worker];
                for (;;) {
                    if (!latency && total_pm_errors.load() >= o.error_target &&
                        total_libra_errors.load() >= o.error_target) break;
                    const uint64_t start = ticket.fetch_add(kBatchSize);
                    if (start >= o.shots) break;
                    const size_t take = std::min<uint64_t>(kBatchSize, o.shots - start);
                    sampler.resample(false);
                    const auto dets = sampler.det_buffer.transposed();
                    const auto obs = sampler.obs_buffer.transposed();
                    for (size_t i = 0; i < take; ++i) {
                        Clock::time_point t0, t1, t2;
                        if (latency) t0 = Clock::now();
                        const auto pm_result = pm.decode(dets[i], obs[i]);
                        if (latency) t1 = Clock::now();
                        const auto libra_result = libra.decode(dets[i], obs[i]);
                        if (latency) t2 = Clock::now();
                        const bool pm_error = pm_result.flipped_obs[0] != obs[i][0];
                        const bool libra_error = libra_result.flipped_obs[0] != obs[i][0];
                        ++c.shots;
                        processed.fetch_add(1);
                        c.pymatching_errors += pm_error;
                        c.libra_errors += libra_error;
                        if (pm_error) total_pm_errors.fetch_add(1);
                        if (libra_error) total_libra_errors.fetch_add(1);
                        c.both_errors += pm_error && libra_error;
                        c.pymatching_only += pm_error && !libra_error;
                        c.libra_only += libra_error && !pm_error;
                        if (latency) {
                            c.pymatching_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
                            c.libra_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
                        }
                        if (!latency && total_pm_errors.load() >= o.error_target &&
                            total_libra_errors.load() >= o.error_target) break;
                    }
                }
                add_stats(c, libra.statistics());
            } catch (const std::exception& e) {
                failed.store(true);
                std::lock_guard<std::mutex> lock(error_mutex);
                if (error_message.empty()) error_message = e.what();
            }
            live_workers.fetch_sub(1);
        });
    }
    auto last_progress = run_started;
    while (live_workers.load() > 0) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        const auto now = Clock::now();
        if (o.progress_seconds &&
            std::chrono::duration_cast<std::chrono::seconds>(now - last_progress).count() >=
                static_cast<int64_t>(o.progress_seconds)) {
            std::cerr << "progress mode=" << o.mode << " completed=" << processed.load()
                      << " requested=" << o.shots << " pymatching_errors="
                      << total_pm_errors.load() << " libra_errors=" << total_libra_errors.load()
                      << " workers_active=" << live_workers.load()
                      << " elapsed_s="
                      << std::chrono::duration<double>(now - run_started).count() << "\n";
            last_progress = now;
        }
    }
    for (auto& t : threads) t.join();
    if (failed) throw std::runtime_error(error_message);

    Counts total;
    for (const auto& c : results) {
        total.shots += c.shots;
        total.pymatching_errors += c.pymatching_errors;
        total.libra_errors += c.libra_errors;
        total.both_errors += c.both_errors;
        total.pymatching_only += c.pymatching_only;
        total.libra_only += c.libra_only;
        total.pymatching_ns += c.pymatching_ns;
        total.libra_ns += c.libra_ns;
        total.ensemble_shots += c.ensemble_shots;
        total.fallback_shots += c.fallback_shots;
        total.improving_cycles += c.improving_cycles;
        total.complementary_ns += c.complementary_ns;
        total.fallback_ns += c.fallback_ns;
        total.ensemble_ns += c.ensemble_ns;
        total.synthesis_ns += c.synthesis_ns;
        total.degeneracy_ns += c.degeneracy_ns;
    }
    write_metadata(out, o, dem, dims);
    out << std::setprecision(12)
        << "wall_time_seconds=" << std::chrono::duration<double>(Clock::now() - run_started).count() << "\n"
        << "shots=" << total.shots << "\n"
        << "pymatching_errors=" << total.pymatching_errors << "\n"
        << "libra_errors=" << total.libra_errors << "\n"
        << "both_fail=" << total.both_errors << "\n"
        << "pymatching_only_fail=" << total.pymatching_only << "\n"
        << "libra_only_fail=" << total.libra_only << "\n"
        << "pymatching_ler=" << double(total.pymatching_errors) / total.shots << "\n"
        << "libra_ler=" << double(total.libra_errors) / total.shots << "\n";
    if (latency) {
        out << "pymatching_mean_decode_ns=" << double(total.pymatching_ns / total.shots) << "\n"
            << "libra_mean_decode_ns=" << double(total.libra_ns / total.shots) << "\n"
            << "libra_ensemble_trigger_fraction=" << double(total.ensemble_shots) / total.shots << "\n"
            << "libra_fallback_trigger_fraction=" << double(total.fallback_shots) / total.shots << "\n"
            << "libra_profile_complementary_ns=" << total.complementary_ns << "\n"
            << "libra_profile_fallback_ns=" << total.fallback_ns << "\n"
            << "libra_profile_ensemble_ns=" << total.ensemble_ns << "\n"
            << "libra_profile_synthesis_ns=" << total.synthesis_ns << "\n"
            << "libra_profile_degeneracy_ns=" << total.degeneracy_ns << "\n"
            << "libra_improving_cycles=" << total.improving_cycles << "\n";
    }
}

struct ErrorPattern { std::vector<uint32_t> detectors; bool observable = false; };

std::vector<ErrorPattern> make_error_patterns(const stim::DetectorErrorModel& dem) {
    std::vector<ErrorPattern> patterns;
    patterns.reserve(dem.count_errors());
    dem.iter_flatten_error_instructions([&](const auto& inst) {
        ErrorPattern p;
        inst.for_separated_targets([&](const auto& group) {
            for (auto t : group) {
                if (t.is_relative_detector_id()) p.detectors.push_back(t.val());
                else if (t.is_observable_id()) p.observable = !p.observable;
            }
        });
        std::sort(p.detectors.begin(), p.detectors.end());
        std::vector<uint32_t> parity;
        for (size_t i = 0; i < p.detectors.size();) {
            size_t j = i + 1;
            while (j < p.detectors.size() && p.detectors[j] == p.detectors[i]) ++j;
            if ((j - i) & 1) parity.push_back(p.detectors[i]);
            i = j;
        }
        p.detectors = std::move(parity);
        patterns.push_back(std::move(p));
    });
    return patterns;
}

void run_rare(const Options& o, const stim::DetectorErrorModel& dem,
              const dec::Libra::Dimensions& dims, std::ostream& out) {
    if (o.level < 0 || static_cast<size_t>(o.level) >= dem.count_errors())
        throw std::invalid_argument("rare LEVEL must be in [0, dem_error_instructions-1]");
    LibraRareEventSampler conditional(dem, size_t(o.level) + 1);
    const double pk = conditional.count_probabilities()[o.level];
    if (pk <= 0) throw std::invalid_argument("requested rare-event level has zero probability");
    const auto patterns = make_error_patterns(dem);
    const uint64_t base = o.shots / o.workers, rem = o.shots % o.workers;
    std::vector<uint64_t> pm_fail(o.workers), libra_fail(o.workers), samples(o.workers);
    std::vector<std::thread> threads;
    std::atomic<uint64_t> processed{0};
    std::atomic<uint64_t> progress_pm_fail{0}, progress_libra_fail{0};
    std::atomic<int> live_workers{o.workers};
    std::mutex error_mutex;
    std::string error_message;
    std::atomic<bool> failed{false};
    const auto run_started = Clock::now();
    for (int w = 0; w < o.workers; ++w) {
        const uint64_t n = base + (uint64_t(w) < rem);
        threads.emplace_back([&, w, n] {
            try {
                auto libra = make_libra(dem, o.seed, o.ensemble_size, false);
                dec::PyMatching pm(dem, false);
                libra.prepare_ensemble();
                std::mt19937_64 rng(o.seed + 0xd1b54a32d192ed03ULL * uint64_t(w + 1));
                for (uint64_t i = 0; i < n; ++i) {
                    if (progress_pm_fail.load() >= o.error_target &&
                        progress_libra_fail.load() >= o.error_target) break;
                    const auto selected = conditional.sample_errors(o.level, rng);
                    dec::SyndromeType syn(dem.count_detectors());
                    dec::ObsType truth(1);
                    for (size_t id : selected) {
                        const auto& p = patterns[id];
                        for (uint32_t d : p.detectors) syn[d] ^= 1;
                        truth[0] ^= p.observable;
                    }
                    const auto pm_result = pm.decode(syn, truth);
                    const auto libra_result = libra.decode(syn, truth);
                    const bool pm_error = pm_result.flipped_obs[0] != truth[0];
                    const bool libra_error = libra_result.flipped_obs[0] != truth[0];
                    pm_fail[w] += pm_error;
                    libra_fail[w] += libra_error;
                    if (pm_error) progress_pm_fail.fetch_add(1);
                    if (libra_error) progress_libra_fail.fetch_add(1);
                    ++samples[w];
                    processed.fetch_add(1);
                }
            } catch (const std::exception& e) {
                failed.store(true);
                std::lock_guard<std::mutex> lock(error_mutex);
                if (error_message.empty()) error_message = e.what();
            }
            live_workers.fetch_sub(1);
        });
    }
    auto last_progress = run_started;
    while (live_workers.load() > 0) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        const auto now = Clock::now();
        if (o.progress_seconds &&
            std::chrono::duration_cast<std::chrono::seconds>(now - last_progress).count() >=
                static_cast<int64_t>(o.progress_seconds)) {
            std::cerr << "progress mode=rare level=" << o.level << " completed=" << processed.load()
                      << " requested=" << o.shots << " pymatching_errors=" << progress_pm_fail.load()
                      << " libra_errors=" << progress_libra_fail.load()
                      << " workers_active=" << live_workers.load()
                      << " elapsed_s="
                      << std::chrono::duration<double>(now - run_started).count()
                      << "\n";
            last_progress = now;
        }
    }
    for (auto& t : threads) t.join();
    if (failed) throw std::runtime_error(error_message);
    uint64_t n = 0, pm_errors = 0, libra_errors = 0;
    for (int w = 0; w < o.workers; ++w) {
        n += samples[w]; pm_errors += pm_fail[w]; libra_errors += libra_fail[w];
    }
    write_metadata(out, o, dem, dims);
    auto report = [&](const char* name, uint64_t errors) {
        const double phat = double(errors) / n;
        const double variance = phat * (1 - phat) / n;
        out << name << "_conditional_failures=" << errors << "\n"
            << name << "_conditional_rate=" << phat << "\n"
            << name << "_conditional_variance=" << variance << "\n"
            << name << "_conditional_standard_error=" << std::sqrt(variance) << "\n"
            << name << "_level_contribution=" << pk * phat << "\n";
    };
    out << std::setprecision(12) << "level=" << o.level << "\n"
        << "wall_time_seconds=" << std::chrono::duration<double>(Clock::now() - run_started).count() << "\n"
        << "level_probability=" << pk << "\n"
        << "samples=" << n << "\n";
    report("pymatching", pm_errors);
    report("libra", libra_errors);
    out << "level_scope=only_explicitly_requested_level; submit all K levels separately for a full sum\n";
}

void run_rare_mix(const Options& o, const stim::DetectorErrorModel& dem,
                  const dec::Libra::Dimensions& dims, std::ostream& out) {
    if (o.max_count == 0 || o.max_count > dem.count_errors() + 1)
        throw std::invalid_argument("rare_mix --max-count must be in [1, DEM instruction count + 1]");
    LibraRareEventSampler conditional(dem, o.max_count);
    const auto& pk = conditional.count_probabilities();
    std::vector<double> log_q(pk.size(), -std::numeric_limits<double>::infinity());
    double max_log_q = -std::numeric_limits<double>::infinity();
    for (size_t k = 0; k < pk.size(); ++k) {
        if (pk[k] > 0) {
            log_q[k] = std::log(pk[k]) + o.theta * static_cast<double>(k);
            max_log_q = std::max(max_log_q, log_q[k]);
        }
    }
    std::vector<double> q(pk.size(), 0);
    double q_sum = 0, covered_mass = 0;
    for (size_t k = 0; k < pk.size(); ++k) {
        if (pk[k] > 0)
            q[k] = std::exp(std::max(log_q[k] - max_log_q,
                                     std::log(std::numeric_limits<double>::min())));
        q_sum += q[k];
        covered_mass += pk[k];
    }
    if (!(q_sum > 0)) throw std::runtime_error("rare_mix proposal has no support");
    for (double& probability : q) probability /= q_sum;
    const double tail_bound = conditional.tail_probability_bound(o.max_count);
    const auto patterns = make_error_patterns(dem);
    const uint64_t base = o.shots / o.workers, rem = o.shots % o.workers;
    struct WorkerStats {
        uint64_t n = 0;
        uint64_t pm_errors = 0, libra_errors = 0;
        long double pm_sum = 0, pm_sum2 = 0, libra_sum = 0, libra_sum2 = 0;
        std::vector<uint64_t> k_n, k_pm_errors, k_libra_errors;
        explicit WorkerStats(size_t levels)
            : k_n(levels), k_pm_errors(levels), k_libra_errors(levels) {}
    };
    std::vector<WorkerStats> results;
    results.reserve(o.workers);
    for (int w = 0; w < o.workers; ++w) results.emplace_back(pk.size());
    std::vector<std::thread> threads;
    std::mutex error_mutex;
    std::string error_message;
    std::atomic<bool> failed{false};
    std::atomic<uint64_t> processed{0}, pm_error_progress{0}, libra_error_progress{0};
    std::atomic<int> live_workers{o.workers};
    const auto run_started = Clock::now();
    for (int w = 0; w < o.workers; ++w) {
        const uint64_t n = base + (uint64_t(w) < rem);
        threads.emplace_back([&, w, n] {
            try {
                auto libra = make_libra(dem, o.seed, o.ensemble_size, false);
                dec::PyMatching pm(dem, false);
                libra.prepare_ensemble();
                std::mt19937_64 rng(o.seed + 0xd1b54a32d192ed03ULL * uint64_t(w + 1));
                std::discrete_distribution<size_t> choose_k(q.begin(), q.end());
                WorkerStats& stats = results[w];
                for (uint64_t i = 0; i < n; ++i) {
                    if (pm_error_progress.load() >= o.error_target &&
                        libra_error_progress.load() >= o.error_target) break;
                    const size_t k = choose_k(rng);
                    const auto selected = conditional.sample_errors(k, rng);
                    dec::SyndromeType syn(dem.count_detectors());
                    dec::ObsType truth(1);
                    for (size_t id : selected) {
                        const auto& p = patterns[id];
                        for (uint32_t d : p.detectors) syn[d] ^= 1;
                        truth[0] ^= p.observable;
                    }
                    const auto pm_result = pm.decode(syn, truth);
                    const auto libra_result = libra.decode(syn, truth);
                    const bool pm_error = pm_result.flipped_obs[0] != truth[0];
                    const bool libra_error = libra_result.flipped_obs[0] != truth[0];
                    const double weight = pk[k] / q[k];
                    const long double pm_x = pm_error ? weight : 0;
                    const long double libra_x = libra_error ? weight : 0;
                    stats.pm_sum += pm_x;
                    stats.pm_sum2 += pm_x * pm_x;
                    stats.libra_sum += libra_x;
                    stats.libra_sum2 += libra_x * libra_x;
                    stats.pm_errors += pm_error;
                    stats.libra_errors += libra_error;
                    ++stats.k_n[k];
                    stats.k_pm_errors[k] += pm_error;
                    stats.k_libra_errors[k] += libra_error;
                    ++stats.n;
                    processed.fetch_add(1);
                    if (pm_error) pm_error_progress.fetch_add(1);
                    if (libra_error) libra_error_progress.fetch_add(1);
                }
            } catch (const std::exception& e) {
                failed.store(true);
                std::lock_guard<std::mutex> lock(error_mutex);
                if (error_message.empty()) error_message = e.what();
            }
            live_workers.fetch_sub(1);
        });
    }
    auto last_progress = run_started;
    while (live_workers.load() > 0) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        const auto now = Clock::now();
        if (o.progress_seconds &&
            std::chrono::duration_cast<std::chrono::seconds>(now - last_progress).count() >=
                static_cast<int64_t>(o.progress_seconds)) {
            std::cerr << "progress mode=rare_mix completed=" << processed.load()
                      << " requested=" << o.shots << " pymatching_errors=" << pm_error_progress.load()
                      << " libra_errors=" << libra_error_progress.load()
                      << " workers_active=" << live_workers.load()
                      << " elapsed_s=" << std::chrono::duration<double>(now - run_started).count()
                      << "\n";
            last_progress = now;
        }
    }
    for (auto& t : threads) t.join();
    if (failed) throw std::runtime_error(error_message);

    WorkerStats total(pk.size());
    for (const auto& s : results) {
        total.n += s.n;
        total.pm_errors += s.pm_errors;
        total.libra_errors += s.libra_errors;
        total.pm_sum += s.pm_sum;
        total.pm_sum2 += s.pm_sum2;
        total.libra_sum += s.libra_sum;
        total.libra_sum2 += s.libra_sum2;
        for (size_t k = 0; k < pk.size(); ++k) {
            total.k_n[k] += s.k_n[k];
            total.k_pm_errors[k] += s.k_pm_errors[k];
            total.k_libra_errors[k] += s.k_libra_errors[k];
        }
    }
    auto write_estimate = [&](const char* name, long double sum, long double sum2) {
        const long double mean = sum / total.n;
        const long double sample_variance = total.n > 1
            ? std::max<long double>(0, (sum2 - sum * sum / total.n) / (total.n - 1))
            : 0;
        const double variance_of_mean = static_cast<double>(sample_variance / total.n);
        out << name << "_partial_ler=" << static_cast<double>(mean) << "\n"
            << name << "_partial_variance=" << variance_of_mean << "\n"
            << name << "_partial_standard_error=" << std::sqrt(variance_of_mean) << "\n"
            << name << "_truncated_support_estimate=" << static_cast<double>(mean) << "\n"
            << name << "_estimate_plus_tail_bound="
            << std::min(1.0, static_cast<double>(mean + tail_bound)) << "\n";
    };
    write_metadata(out, o, dem, dims);
    out << std::setprecision(12) << "samples=" << total.n << "\n"
        << "wall_time_seconds=" << std::chrono::duration<double>(Clock::now() - run_started).count() << "\n"
        << "proposal=proportional_to_max_P_K_times_exp_theta_K_floor_double_min\n"
        << "proposal_theta=" << o.theta << "\n"
        << "stopping_rule=stop_after_both_decoders_reach_max_errors_or_requested_shots\n"
        << "covered_count_mass=" << covered_mass << "\n"
        << "omitted_tail_probability_bound=" << tail_bound << "\n"
        << "tail_bound_scope=LER contribution omitted above max_count; statistical interval from sampling error reported separately\n"
        << "pymatching_failures=" << total.pm_errors << "\n"
        << "libra_failures=" << total.libra_errors << "\n"
        << "pymatching_importance_sum=" << static_cast<double>(total.pm_sum) << "\n"
        << "pymatching_importance_sum_squares=" << static_cast<double>(total.pm_sum2) << "\n"
        << "libra_importance_sum=" << static_cast<double>(total.libra_sum) << "\n"
        << "libra_importance_sum_squares=" << static_cast<double>(total.libra_sum2) << "\n";
    write_estimate("pymatching", total.pm_sum, total.pm_sum2);
    write_estimate("libra", total.libra_sum, total.libra_sum2);
    for (size_t k = 0; k < pk.size(); ++k) {
        if (!total.k_n[k]) continue;
        const double pm_rate = double(total.k_pm_errors[k]) / total.k_n[k];
        const double libra_rate = double(total.k_libra_errors[k]) / total.k_n[k];
        out << "K." << k << ".P=" << pk[k] << "\n"
            << "K." << k << ".q=" << q[k] << "\n"
            << "K." << k << ".samples=" << total.k_n[k] << "\n"
            << "K." << k << ".pymatching_failures=" << total.k_pm_errors[k] << "\n"
            << "K." << k << ".libra_failures=" << total.k_libra_errors[k] << "\n"
            << "K." << k << ".pymatching_conditional_rate=" << pm_rate << "\n"
            << "K." << k << ".pymatching_conditional_variance="
            << pm_rate * (1 - pm_rate) / total.k_n[k] << "\n"
            << "K." << k << ".libra_conditional_rate=" << libra_rate << "\n"
            << "K." << k << ".libra_conditional_variance="
            << libra_rate * (1 - libra_rate) / total.k_n[k] << "\n";
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_args(argc, argv);
        const auto circuit = make_circuit(options);
        const auto dem = stim::circuit_to_dem(circuit, {.decompose_errors = true});
        dec::Libra probe(dem, dec::Libra::Config{.ensemble_size=0});
        const auto dims = probe.dimensions();
        std::ofstream output(options.output);
        if (!output) throw std::runtime_error("cannot open output file: " + options.output);
        if (options.mode == "rare") run_rare(options, dem, dims, output);
        else if (options.mode == "rare_mix") run_rare_mix(options, dem, dims, output);
        else run_monte_carlo(options, dem, dims, output, options.mode == "latency");
    } catch (const std::exception& e) {
        std::cerr << "libra_evaluate: " << e.what() << "\n";
        return 2;
    }
    return 0;
}
