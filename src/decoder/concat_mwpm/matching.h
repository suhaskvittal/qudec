/*
 *  author: OpenAI Codex
 *  date:   22 October 2026
 *
 *  Lee, Li, Bartlett, Quantum 9, 1609 (2025), Secs. 3 and 4.2.
 *  C++ concatenated MWPM decoder; PyMatching performs MWPM.
 * */
#ifndef DECODER_CONCAT_MWPM_MATCHING_h
#define DECODER_CONCAT_MWPM_MATCHING_h

#include <stim.h>
#include <pymatching/sparse_blossom/driver/mwpm_decoding.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <vector>

namespace dec::concat
{

struct Event
{
    double probability;
    std::vector<size_t> detectors;
    uint64_t observable;
};

struct Model
{
    std::vector<int> colors;
    std::vector<Event> events;
};

struct Result
{
    uint64_t observable;
    double weight;
    int color;
    std::vector<size_t> events;
    std::vector<std::vector<size_t>> restricted_edges;
    // Ordered original-model event IDs on each second-stage matched path.
    std::vector<std::vector<size_t>> second_stage_paths;
};

class Decoder
{
public:
    // Section 3 uses uniform restricted-lattice weights for code capacity.
    // Algorithm 1 uses compressed projected probabilities for noisy readout.
    explicit Decoder(Model, bool uniform_restricted_weights = false);
    Result decode(const std::vector<uint8_t>& syndrome);
    std::array<Result, 3> decode_colors(const std::vector<uint8_t>& syndrome, bool trace = true);
    // Verify the actual integer matching/search weights against log odds.
    size_t audit_weights() const;
    const Model model;

private:
    struct Graph;
    struct SubDecoder;
    std::array<std::shared_ptr<SubDecoder>, 3> subdecoders_;
};

} // namespace dec::concat
#endif
