// Author: OpenAI Codex. Shared Tesseract options for qudec and qudec_chains.
#ifndef QUDEC_TESSERACT_OPTIONS_H
#define QUDEC_TESSERACT_OPTIONS_H

#include <tesseract.h>
#include <cstdint>
#include <iosfwd>
#include <string_view>

class ARGPARSE;

struct TesseractOptions
{
    int64_t beam{tesseract_decoder::DEFAULT_DET_BEAM};
    int64_t pqlimit{tesseract_decoder::DEFAULT_PQLIMIT};
    int64_t num_orders{1},order_seed{0};
    bool beam_climbing{false},allow_revisit{false},no_merge{false},verbose{false};
    bool order_bfs{false},order_index{false},order_coordinate{false};
    double det_penalty{0};
    bool sparsify{false};
    int64_t base_degree{-1},max_degree{-1},reactivate_limit{-1};

    void add_arguments(ARGPARSE&);
    void validate(std::string_view decoder,int argc,char* const argv[]) const;
    tesseract_decoder::TesseractConfig make_config(const stim::DetectorErrorModel&) const;
    void print(std::ostream&) const;
};

#endif
