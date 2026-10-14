/*
 *  contributor: OpenAI Codex
 * */

#include "decoder/pdz_libra.h"

#include <ostream>
#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include <fcntl.h>
#include <unistd.h>
#include <zlib.h>

#if defined(ENABLE_MPI)
#include <mpi.h>
#endif

namespace dec
{

// Versioned little-endian wire format; never serialize compiler-dependent
// structs. See docs/pdz-syndrome-recording.md for the header and record layout.
struct PDZLibra::SyndromeRecorder
{
    gzFile file{nullptr};
    const size_t detector_count;
    const size_t syndrome_bytes;
    std::vector<uint8_t> buffer;

    static void put_integer(uint8_t* destination, uint64_t value, size_t bytes)
    {
        for (size_t i = 0; i < bytes; i++)
            destination[i] = uint8_t(value >> (8 * i));
    }

    SyndromeRecorder(const std::string& path, size_t detectors, size_t observables, uint32_t rank)
        :detector_count(detectors), syndrome_bytes((detectors + 7) / 8),
        buffer(24 + 2 * syndrome_bytes, 0)
    {
        if (detectors > std::numeric_limits<uint32_t>::max() ||
            buffer.size() > std::numeric_limits<int>::max())
            throw std::invalid_argument("PDZ syndrome recording exceeds the binary format limits");
        const int fd = open(path.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0644);
        if (fd < 0)
            throw std::runtime_error("Cannot create PDZ syndrome file (must not exist): " + path);
        file = gzdopen(fd, "wb1");
        if (!file)
        {
            close(fd);
            throw std::runtime_error("Cannot initialize gzip for PDZ syndrome file: " + path);
        }
        // gzbuffer must precede the first write.
        if (gzbuffer(file, 64 * 1024) != 0)
        {
            gzclose(file);
            file = nullptr;
            throw std::runtime_error("Cannot allocate PDZ syndrome gzip buffer");
        }
        uint8_t header[40]{};
        const char magic[] = "PDZLIBR1";
        std::copy_n(magic, 8, header);
        put_integer(header + 8, 2, 4);
        put_integer(header + 12, buffer.size(), 4);
        put_integer(header + 16, detectors, 8);
        put_integer(header + 24, observables, 8);
        put_integer(header + 32, rank, 4);
        put_integer(header + 36, sizeof(header), 4);
        if (gzwrite(file, header, sizeof(header)) != sizeof(header))
        {
            gzclose(file);
            file = nullptr;
            throw std::runtime_error("Cannot write PDZ syndrome file header");
        }
    }

    ~SyndromeRecorder()
    {
        if (file)
            gzclose(file);
    }

    void pack_syndrome(SyndromeRef syndrome, uint8_t* destination) const
    {
        for (size_t i = 0; i < syndrome_bytes; i++)
            destination[i] = uint8_t(syndrome.u64[i / 8] >> (8 * (i % 8)));
        if (detector_count % 8)
            destination[syndrome_bytes - 1] &= uint8_t((1u << (detector_count % 8)) - 1);
    }

    void add(SyndromeRef input, SyndromeRef remaining, size_t input_hw, size_t post_hw,
             double gap_db, bool ensemble, bool changed, bool initial_frame,
             bool final_frame, bool true_frame, bool bypassed)
    {
        put_integer(buffer.data(), input_hw, 4);
        put_integer(buffer.data() + 4, post_hw, 4);
        put_integer(buffer.data() + 8, std::bit_cast<uint64_t>(gap_db), 8);
        buffer[16] = ensemble;
        buffer[17] = changed;
        buffer[18] = initial_frame;
        buffer[19] = final_frame;
        buffer[20] = true_frame;
        buffer[21] = bypassed;
        pack_syndrome(input, buffer.data() + 24);
        pack_syndrome(remaining, buffer.data() + 24 + syndrome_bytes);
        if (gzwrite(file, buffer.data(), unsigned(buffer.size())) != int(buffer.size()))
            throw std::runtime_error("Cannot write PDZ syndrome record");
    }

    void finish()
    {
        const int status = gzclose(file);
        file = nullptr;
        if (status != Z_OK)
            throw std::runtime_error("Cannot finish PDZ syndrome gzip file");
    }
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

PDZLibra::PDZLibra(const stim::DetectorErrorModel& dem, pp::PDZ::config_type pdz_config,
                  Libra::Config libra_config, double removal_threshold)
    :s_pdz_removed_pairs("PDZ_LIBRA_REMOVED_PAIRS", 0,
                         double(dem.count_detectors() / 2 + 1),
                         dem.count_detectors() / 2 + 1),
    s_ensemble_by_pdz_removed_pairs("PDZ_LIBRA_ENSEMBLE_BY_REMOVED_PAIRS", 0,
                                    double(dem.count_detectors() / 2 + 1),
                                    dem.count_detectors() / 2 + 1),
    s_ensemble_frame_changes_by_pdz_removed_pairs("PDZ_LIBRA_ENSEMBLE_FRAME_CHANGES_BY_REMOVED_PAIRS", 0,
                                                 double(dem.count_detectors() / 2 + 1),
                                                 dem.count_detectors() / 2 + 1),
    s_pdz_same_basis_removed_pairs("PDZ_LIBRA_SAME_BASIS_REMOVED_PAIRS", 0,
                                    double(dem.count_detectors() / 2 + 1),
                                    dem.count_detectors() / 2 + 1),
    s_ensemble_by_pdz_same_basis_removed_pairs("PDZ_LIBRA_ENSEMBLE_BY_SAME_BASIS_REMOVED_PAIRS", 0,
                                               double(dem.count_detectors() / 2 + 1),
                                               dem.count_detectors() / 2 + 1),
    s_ensemble_frame_changes_by_pdz_same_basis_removed_pairs(
            "PDZ_LIBRA_ENSEMBLE_FRAME_CHANGES_BY_SAME_BASIS_REMOVED_PAIRS", 0,
            double(dem.count_detectors() / 2 + 1), dem.count_detectors() / 2 + 1),
    s_post_pdz_hamming_weight("PDZ_LIBRA_POST_PDZ_HAMMING_WEIGHT", 0,
                              double(dem.count_detectors() + 1), dem.count_detectors() + 1),
    s_ensemble_by_post_pdz_hamming_weight("PDZ_LIBRA_ENSEMBLE_BY_POST_PDZ_HAMMING_WEIGHT", 0,
                                          double(dem.count_detectors() + 1), dem.count_detectors() + 1),
    s_ensemble_frame_changes_by_post_pdz_hamming_weight(
            "PDZ_LIBRA_ENSEMBLE_FRAME_CHANGES_BY_POST_PDZ_HAMMING_WEIGHT", 0,
            double(dem.count_detectors() + 1), dem.count_detectors() + 1),
    initial_(dem, {.enable_gap_estimation=false, .report_matching_data=true}),
    residual_(dem, libra_config),
    pipeline_(dem, &initial_, &residual_, pdz_config, removal_threshold, "pymatching", "libra"),
    removal_threshold_(removal_threshold)
{
    if (!std::isfinite(removal_threshold_) || removal_threshold_ < 0 || removal_threshold_ > 1)
        throw std::invalid_argument("PDZ removal threshold must be in [0, 1]");
    residual_.prepare_ensemble();
}

PDZLibra::~PDZLibra() = default;

void
PDZLibra::record_syndromes(const std::string& path, uint32_t rank)
{
    if (recorder_ || counts_[0])
        throw std::logic_error("PDZ syndrome recording must start once, before decoding");
    recorder_ = std::make_unique<SyndromeRecorder>(path, pipeline_.predecoder().num_detectors,
                                                  pipeline_.predecoder().num_observables, rank);
}

void
PDZLibra::finish_recording()
{
    if (recorder_)
    {
        recorder_->finish();
        recorder_.reset();
    }
}

Result
PDZLibra::decode(SyndromeRef syndrome, ObsRef truth)
{
    // Preserve the original sample for the sampler and its error callbacks.
    auto decoded = pipeline_.decode_detailed(syndrome, truth);
    const auto& initial = decoded.initial;
    const auto& remaining = decoded.remaining;
    const auto& partial = decoded.partial;
    const auto input_events = syndrome.popcnt();
    const auto residual_events = remaining.popcnt();
    const double removed_pairs = double(input_events - residual_events) / 2.0;
    size_t same_basis_removed_events = 0;
    for (const auto& assignment : partial.matching_data.assignments)
        for (auto d : {assignment.d1, assignment.d2})
            if (d >= 0 && !pipeline_.predecoder().is_opposing_basis_detector(d))
                same_basis_removed_events++;
    const double same_basis_removed_pairs = double(same_basis_removed_events) / 2.0;

    const bool bypassed = decoded.bypassed;
    auto result = std::move(decoded.result);
    const bool ensembled = decoded.ensembled;
    if (bypassed)
        counts_[4]++;
    s_pdz_removed_pairs.add(removed_pairs);
    s_post_pdz_hamming_weight.add(residual_events);
    s_pdz_same_basis_removed_pairs.add(same_basis_removed_pairs);
    const bool frame_changed = result.flipped_obs != initial.flipped_obs;
    if (ensembled)
    {
        s_ensemble_by_post_pdz_hamming_weight.add(residual_events);
        s_ensemble_by_pdz_removed_pairs.add(removed_pairs);
        s_ensemble_by_pdz_same_basis_removed_pairs.add(same_basis_removed_pairs);
        if (frame_changed)
        {
            s_ensemble_frame_changes_by_post_pdz_hamming_weight.add(residual_events);
            s_ensemble_frame_changes_by_pdz_removed_pairs.add(removed_pairs);
            s_ensemble_frame_changes_by_pdz_same_basis_removed_pairs.add(same_basis_removed_pairs);
        }
    }
    // Record an unavailable gap for bypassed shots without changing the
    // returned original PyMatching Result or running confidence estimation.
    if (recorder_)
        recorder_->add(syndrome, remaining, input_events, residual_events,
                       bypassed ? std::numeric_limits<double>::quiet_NaN() : result.matching_data.gap,
                       ensembled, ensembled && frame_changed, initial.flipped_obs[0],
                       result.flipped_obs[0], truth[0], bypassed);

    bool failed = false;
    for (size_t i = 0; i < pipeline_.predecoder().num_observables; i++)
        failed |= result.flipped_obs[i] != truth[i];
    counts_[0]++;
    counts_[1] += failed;
    counts_[2] += input_events;
    counts_[3] += residual_events;
    return result;
}

PDZLibra::Statistics
PDZLibra::statistics() const
{
    // Derive calls from the MPI-reduced shot and bypass counts, keeping
    // sampled_shots == bypassed_shots + libra_called_shots on every rank.
    const auto calls = counts_[0] - counts_[4];
    const double fraction = counts_[0] ? double(calls) / double(counts_[0]) : 0.0;
    const auto timing = residual_.statistics();
    const double mean_us = timing.timed_shots
        ? double(timing.decode_total_ns) / double(timing.timed_shots) / 1000.0 : 0.0;
    return {counts_[0], counts_[1], counts_[2], counts_[3], counts_[4], calls,
            fraction, mean_us, fraction * mean_us};
}

void
PDZLibra::print_stats(std::ostream& out) const
{
    const auto stats = statistics();
    out << "PDZLibra sampled_shots=" << stats.sampled_shots
        << " sampled_failures=" << stats.sampled_failures
        << " input_detection_events=" << stats.input_events
        << " residual_detection_events=" << stats.residual_events
        << " libra_bypassed_shots=" << stats.bypassed_shots
        << " libra_called_shots=" << stats.libra_called_shots
        << " pdz_removal_threshold=" << removal_threshold_ << '\n';
    out << "Libra weighted decode latency (us): called_fraction=" << stats.libra_called_fraction
        << " mean_when_called=" << stats.libra_mean_latency_us
        << " weighted=" << stats.libra_weighted_latency_us
        << " scope=full_history excludes_initial_pymatching=true\n";
    /*
    out << "PDZ histograms: X = removed syndrome bits / 2; bin k covers "
           "2k or 2k+1 removed bits (boundary matches can remove one bit)\n";
    out << s_pdz_removed_pairs.to_string_full();
    out << s_ensemble_by_pdz_removed_pairs.to_string_full();
    out << "PDZ ensemble frame-change histogram: final PDZ + Libra frame differs "
           "from initial uncorrelated PyMatching frame\n";
    out << s_ensemble_frame_changes_by_pdz_removed_pairs.to_string_full();
    out << "PDZ same-basis histograms: X = removed memory-basis syndrome bits / 2; "
           "opposing-basis pairs excluded\n";
    out << s_pdz_same_basis_removed_pairs.to_string_full();
    out << s_ensemble_by_pdz_same_basis_removed_pairs.to_string_full();
    out << s_ensemble_frame_changes_by_pdz_same_basis_removed_pairs.to_string_full();
    out << "PDZ Hamming-weight histograms: X = post-PDZ syndrome popcount "
           "including opposing-basis detectors\n";
    out << s_post_pdz_hamming_weight.to_string_full();
    out << s_ensemble_by_post_pdz_hamming_weight.to_string_full();
    out << s_ensemble_frame_changes_by_post_pdz_hamming_weight.to_string_full();
    initial_.print_stats(out);
    residual_.print_stats(out);
    */
}

void
PDZLibra::mpi_accumulate()
{
#if defined(ENABLE_MPI)
    MPI_Allreduce(MPI_IN_PLACE, counts_.data(), counts_.size(), MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
#endif
    s_pdz_removed_pairs.mpi_accumulate();
    s_ensemble_by_pdz_removed_pairs.mpi_accumulate();
    s_ensemble_frame_changes_by_pdz_removed_pairs.mpi_accumulate();
    s_pdz_same_basis_removed_pairs.mpi_accumulate();
    s_ensemble_by_pdz_same_basis_removed_pairs.mpi_accumulate();
    s_ensemble_frame_changes_by_pdz_same_basis_removed_pairs.mpi_accumulate();
    s_post_pdz_hamming_weight.mpi_accumulate();
    s_ensemble_by_post_pdz_hamming_weight.mpi_accumulate();
    s_ensemble_frame_changes_by_post_pdz_hamming_weight.mpi_accumulate();
    pipeline_.mpi_accumulate();
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace dec
