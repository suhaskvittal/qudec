/*
 *  contributor: OpenAI Codex (decoder types migration).
 *  author: OpenAI GPT-6-Luna
 *  date:   23 Sep 2026
 *  purpose: Implement the Tesseract wrapper in its dedicated decoder translation unit.
 * */

#include "decoder/tesseract.h"

#include <vector>
#include <iomanip>
#include <sstream>
#include <iostream>

namespace dec
{

Tesseract::Tesseract(const stim::DetectorErrorModel& dem)
    :num_detectors(dem.count_detectors()),
    num_observables(dem.count_observables()),
    decoder_(tesseract_decoder::TesseractConfig{.dem=dem})
{}

// OpenAI GPT-6-Sol: Preserve the existing default constructor while accepting
// explicit search settings for reproducible benchmark comparisons.
Tesseract::Tesseract(tesseract_decoder::TesseractConfig config)
    :num_detectors(config.dem.count_detectors()),
    num_observables(config.dem.count_observables()),
    decoder_(std::move(config))
{}

Result
Tesseract::decode(SyndromeRef syn, ObsRef)
{
    // Collect indices of fired detectors -- `TesseractDecoder::decode` wants
    // the sparse list of detection events rather than a dense bitvector.
    std::vector<uint64_t> det_events;
    for (size_t i = 0; i < num_detectors; i++)
        if (syn[i])
            det_events.push_back(i);

    // `decode()` returns the sorted list of flipped observable ids.
    auto flipped = decoder_.decode(det_events);

    Result out{.flipped_obs=ObsType(num_observables)};
    for (int k : flipped)
        out.flipped_obs[k] ^= 1;
    // OpenAI Codex: Report retained hyperedges using original flattened DEM IDs.
    if (debug_ || decoder_.config.verbose)
    {
        std::ostringstream text;
        text << std::setprecision(17)
             << "Tesseract final solution (original flattened DEM error IDs)\n"
             << "Total cost=" << decoder_.cost_from_errors(decoder_.predicted_errors_buffer)
             << " low_confidence=" << decoder_.low_confidence_flag << '\n';
        for (auto id:decoder_.predicted_errors_buffer)
        {
            const auto& error=decoder_.errors.at(decoder_.dem_error_to_error.at(id));
            text << "DEM_E" << id << " p=" << error.get_probability()
                 << " cost=" << error.likelihood_cost << " detectors={";
            for (auto d:error.symptom.detectors) text << " D" << d;
            text << " } observables={";
            for (auto l:error.symptom.observables) text << " L" << l;
            text << " }\n";
        }
        out.debug_info=text.str();
        if (decoder_.config.verbose) std::cout << out.debug_info;
    }
    return out;
}

} // namespace dec
