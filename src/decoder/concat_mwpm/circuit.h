/*
 *  author: OpenAI Codex
 *  date:   22 October 2026
 *
 *  Bell-flagged SI1000 adapter for concatenated MWPM.
 * */
#ifndef DECODER_CONCAT_MWPM_CIRCUIT_h
#define DECODER_CONCAT_MWPM_CIRCUIT_h

#include "decoder/concat_mwpm/matching.h"
#include <string>
#include <tuple>

namespace dec::concat
{

struct CircuitFault
{
    double probability;
    std::vector<size_t> detectors;
    uint64_t observable;
};

struct FlagWeight
{
    std::vector<size_t> checks; // Local memory-sector detector IDs.
    uint64_t observable;
    double unflagged_probability;
    double effective_probability;
    size_t relevant_flags;
    size_t fired_flags;
};

struct CircuitResult
{
    Result correction;
    std::array<Result, 3> colors;
    std::vector<FlagWeight> weights;
    Model model;
    size_t audited_graph_edges;
};

class CircuitDecoder
{
public:
    // Uses the full original DEM and infers each flag's isolated-flip prior.
    CircuitDecoder(const stim::DetectorErrorModel&, bool memory_x,
                   bool use_flags = true, size_t cache_capacity = 16);
    uint64_t decode(const std::vector<uint8_t>& full_syndrome);
    CircuitResult inspect(const std::vector<uint8_t>& full_syndrome);
    std::vector<uint8_t> fault_syndrome(const std::vector<size_t>& ids) const;
    uint64_t fault_observable(const std::vector<size_t>& ids) const;

    const size_t num_detectors;
    const std::vector<CircuitFault> faults;
    const std::map<uint64_t, std::vector<double>> coordinates;
    std::vector<size_t> check_globals, flag_globals;
    size_t class_count() const { return classes_.size(); }
    const std::vector<double>& flag_measurement_probabilities() const { return measurement_p_; }

private:
    struct Variant
    {
        std::vector<size_t> flags;
        double probability;
    };
    struct Class
    {
        std::vector<size_t> checks;
        uint64_t observable;
        std::vector<size_t> flags;
        std::vector<Variant> variants;
    };
    struct Context
    {
        std::unique_ptr<Decoder> decoder;
        std::vector<FlagWeight> weights;
    };
    std::vector<double> measurement_p_; // Indexed by original global flag ID.
    const bool use_flags_;
    const size_t cache_capacity_;
    std::vector<int> colors_;
    std::vector<Class> classes_;
    std::map<std::string, Context> contexts_;

    Context make_context(const std::vector<uint8_t>&) const;
    Context& context(const std::vector<uint8_t>&);
    std::vector<uint8_t> memory_syndrome(const std::vector<uint8_t>&) const;
};
} // namespace dec::concat
#endif
