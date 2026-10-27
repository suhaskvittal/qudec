/*
 *  author: OpenAI Codex
 *  date:   27 October 2026
 *  Auxiliary-detector sum-product preprocessing; memory checks remain with MWPM.
 * */
#ifndef DECODER_CONCAT_MWPM_FLAG_BP_h
#define DECODER_CONCAT_MWPM_FLAG_BP_h

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace dec::concat
{
struct FlagBPVariable
{
    double probability;
    std::vector<size_t> flags; // Original auxiliary detector IDs: flags or opposite-basis checks.
    size_t class_id{SIZE_MAX}; // Pure auxiliary faults have no matching event.
    uint64_t observable{0};
};

class FlagBP
{
public:
    FlagBP(size_t num_detectors,std::vector<FlagBPVariable>);
    std::vector<double> posterior(const std::vector<uint8_t>& syndrome) const;
    const std::vector<FlagBPVariable>& variables() const { return variables_; }
    bool independent_stars() const { return stars_; }
private:
    std::vector<FlagBPVariable> variables_;
    std::vector<std::vector<size_t>> checks_,variable_edges_;
    std::vector<size_t> owners_;
    std::vector<std::array<double,2>> star_posterior_;
    std::vector<double> prior_llr_;
    bool stars_{true};
};
} // namespace dec::concat
#endif
