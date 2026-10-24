/* Author: OpenAI Codex. Connected, fixed-weight DEM stress tests. */
#ifndef ERROR_CHAIN_H
#define ERROR_CHAIN_H

#include "decoder/types.h"
#include <functional>
#include <iosfwd>
#include <string>

struct ErrorChainConfig
{
    size_t length{2};
    uint64_t samples{10000};
    uint64_t seed{0};
    bool exhaustive{false};
    uint64_t print_failures{1};
    bool print_progress{false};
    uint64_t print_chains{0};
    std::string code{};
    uint32_t distance{0};
    bool superdense{false}; // Author: OpenAI Codex. Select matching display geometry.
};
struct ErrorChainSummary
{
    uint64_t tested{0};
    uint64_t failures{0};
    uint64_t combinations_visited{0};
};

// Exact decimal binomial count, including counts exceeding uint64_t.
std::string error_chain_combination_count(size_t events, size_t length);
using ErrorChainDecode = std::function<dec::Result(dec::SyndromeRef, dec::ObsRef)>;
ErrorChainSummary evaluate_error_chains(const stim::DetectorErrorModel&, const ErrorChainConfig&,
    const ErrorChainDecode&, std::ostream&);

#endif
