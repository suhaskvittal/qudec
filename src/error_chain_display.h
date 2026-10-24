/* Author: OpenAI Codex. Color-code detector identities and ASCII syndromes. */
#ifndef ERROR_CHAIN_DISPLAY_H
#define ERROR_CHAIN_DISPLAY_H

#include "decoder/types.h"
#include <cstdint>
#include <iosfwd>
#include <map>
#include <string>
#include <vector>

class ErrorChainDisplay
{
public:
    ErrorChainDisplay(const stim::DetectorErrorModel&, const std::string& code, uint32_t distance,
        bool superdense = false);
    std::string detector_label(size_t id) const;
    void draw(dec::SyndromeRef syndrome, std::ostream&) const;
private:
    struct Check
    {
        double x, y;
        int color;
        std::vector<uint32_t> qubits;
    };
    std::string code_;
    std::map<uint64_t, std::vector<double>> coordinates_;
    std::map<uint64_t, std::vector<double>> qubits_;
    std::vector<Check> checks_;
    std::map<size_t,size_t> detector_checks_;
};

#endif
