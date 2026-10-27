/*
 *  contributor: OpenAI Codex (decoder types migration).
 *  author: OpenAI GPT-6-Luna
 *  date:   23 Sep 2026
 *  purpose: Isolate the Tesseract decoder wrapper from the surface-code decoders.
 * */

#ifndef DECODER_TESSERACT_h
#define DECODER_TESSERACT_h

#include "decoder/types.h"

#include <stim.h>
#include <tesseract.h>

namespace dec
{

/*
 * `Tesseract` wraps Google's Tesseract decoder (vendored under `deps/tesseract`,
 * see https://github.com/quantumlib/tesseract-decoder). Tesseract is an A*-style
 * search decoder: it explores subsets of error mechanisms consistent with the
 * observed syndrome and returns the lowest-cost (most likely) explanation.
 *
 * OpenAI Codex: Optional diagnostics expose the final selected hyperedges in
 * Result::debug_info. Search verbosity also prints these diagnostics. Matching
 * data is left at its default value because this is a hypergraph search.
 * */
class Tesseract
{
public:
    const size_t num_detectors;
    const size_t num_observables;
public:
    Tesseract(const stim::DetectorErrorModel&);
    // OpenAI GPT-6-Sol: Let controlled comparisons supply paper benchmark settings.
    Tesseract(tesseract_decoder::TesseractConfig);

    Result decode(SyndromeRef, ObsRef);

    // OpenAI Codex: Enable final-solution diagnostics without verbose search output.
    void set_debug(bool enabled=true) { debug_=enabled; }

    // OpenAI GPT-6-Sol: Expose Tesseract's per-decode confidence signal for comparisons.
    bool last_low_confidence() const { return decoder_.low_confidence_flag; }

    void print_stats(std::ostream&) const {}
    void mpi_accumulate() {}
private:
    tesseract_decoder::TesseractDecoder decoder_;
    bool debug_{false};
};

} // namespace dec

#endif // DECODER_TESSERACT_h
