/*
 *  author: OpenAI GPT-6-Luna
 *  date:   23 Sep 2026
 *  purpose: Implement the Tesseract wrapper in its dedicated decoder translation unit.
 * */

#include "decoder/tesseract.h"

#include <vector>

namespace decoder
{

Tesseract::Tesseract(const stim::DetectorErrorModel& dem)
    :num_detectors(dem.count_detectors()),
    num_observables(dem.count_observables()),
    decoder_(tesseract_decoder::TesseractConfig{.dem=dem})
{}

result_type
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

    result_type out{.flipped_obs=ObsType(num_observables)};
    for (int k : flipped)
        out.flipped_obs[k] ^= 1;
    return out;
}

} // namespace decoder
