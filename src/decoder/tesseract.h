/*
 *  author: OpenAI GPT-6-Luna
 *  date:   23 Sep 2026
 *  purpose: Isolate the Tesseract decoder wrapper from the surface-code decoders.
 * */

#ifndef DECODER_TESSERACT_h
#define DECODER_TESSERACT_h

#include "decoder/common.h"

#include <stim.h>
#include <tesseract.h>

namespace decoder
{

/*
 * `Tesseract` wraps Google's Tesseract decoder (vendored under `deps/tesseract`,
 * see https://github.com/quantumlib/tesseract-decoder). Tesseract is an A*-style
 * search decoder: it explores subsets of error mechanisms consistent with the
 * observed syndrome and returns the lowest-cost (most likely) explanation.
 *
 * This wrapper only exposes the plain single-shot decode path
 * (`tesseract_decoder::TesseractDecoder::decode`); Tesseract's complementary-gap /
 * low-confidence signals, multi-pass decoding, and visualization are out of scope
 * for now, so `result_type::matching_data` is left at its default value.
 * */
class Tesseract
{
public:
    const size_t num_detectors;
    const size_t num_observables;
public:
    Tesseract(const stim::DetectorErrorModel&);

    result_type decode(SyndromeRef, ObsRef);

    void print_stats(std::ostream&) const {}
    void mpi_accumulate() {}
private:
    tesseract_decoder::TesseractDecoder decoder_;
};

} // namespace decoder

#endif // DECODER_TESSERACT_h
