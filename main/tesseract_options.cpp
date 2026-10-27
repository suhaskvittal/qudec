// Author: OpenAI Codex. CLI-to-decoder configuration without changing defaults.
#include "tesseract_options.h"
#include <argparse/argparse.h>
#include <cmath>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <vector>

void TesseractOptions::add_arguments(ARGPARSE& parser)
{
    parser.optional("", "--tesseract-beam", "Tesseract: detector beam (-1 = unrestricted sentinel)", beam, tesseract_decoder::DEFAULT_DET_BEAM)
        .optional("", "--tesseract-beam-climbing", "Tesseract: cycle beam values from zero to the configured beam", beam_climbing, false)
        .optional("", "--tesseract-pqlimit", "Tesseract: priority-queue push limit (-1 = unlimited)", pqlimit, tesseract_decoder::DEFAULT_PQLIMIT)
        .optional("", "--tesseract-num-det-orders", "Tesseract: total detector orders, shared across selected methods", num_orders, 1)
        .optional("", "--tesseract-det-order-seed", "Tesseract: fixed detector-order seed, shared across MPI ranks", order_seed, 0)
        .optional("", "--tesseract-det-order-bfs", "Tesseract: include BFS detector orders", order_bfs, false)
        .optional("", "--tesseract-det-order-index", "Tesseract: include index orders (default when none selected)", order_index, false)
        .optional("", "--tesseract-det-order-coordinate", "Tesseract: include coordinate detector orders", order_coordinate, false)
        .optional("", "--tesseract-allow-revisit-dets", "Tesseract: allow revisiting previously processed detectors", allow_revisit, false)
        .optional("", "--tesseract-no-merge-errors", "Tesseract: disable merging of identical error mechanisms", no_merge, false)
        .optional("", "--tesseract-det-penalty", "Tesseract: extra nonnegative cost per residual detector", det_penalty, 0.0)
        .optional("", "--tesseract-verbose", "Tesseract: print per-search diagnostics", verbose, false)
        .optional("", "--tesseract-sparsify-errors", "Tesseract: enable error sparsification (requires base degree)", sparsify, false)
        .optional("", "--tesseract-sparsify-base-degree", "Tesseract: always retain errors of degree at most this value", base_degree, -1)
        .optional("", "--tesseract-sparsify-max-degree", "Tesseract: maximum optional error degree (-1 = unlimited)", max_degree, -1)
        .optional("", "--tesseract-sparsify-reactivate-limit", "Tesseract: optional errors to reactivate (-1 = automatic)", reactivate_limit, -1);
}

void TesseractOptions::validate(std::string_view decoder,int argc,char* const argv[]) const
{
    for (int j=1;j<argc;j++)
        if (std::string_view(argv[j]).starts_with("--tesseract-") && decoder!="tesseract")
            throw std::invalid_argument("--tesseract-* options require the tesseract decoder");
    if (beam < -1 || beam>int64_t(tesseract_decoder::INF_DET_BEAM))
        throw std::invalid_argument("--tesseract-beam must be -1 or between 0 and 65535");
    if (pqlimit < -1 || pqlimit==0)
        throw std::invalid_argument("--tesseract-pqlimit must be positive or -1 for unlimited");
    if (num_orders<1 || num_orders>std::numeric_limits<int>::max() || order_seed<0)
        throw std::invalid_argument("Tesseract detector-order count must be positive and fit int; seed must be nonnegative");
    const int methods=int(order_bfs)+int(order_index)+int(order_coordinate);
    if (num_orders<methods)
        throw std::invalid_argument("--tesseract-num-det-orders must cover every selected ordering method");
    if (!std::isfinite(det_penalty) || det_penalty<0)
        throw std::invalid_argument("--tesseract-det-penalty must be finite and nonnegative");
    for (auto value:{base_degree,max_degree,reactivate_limit})
        if (value < -1 || value>std::numeric_limits<int>::max())
            throw std::invalid_argument("Tesseract sparsification values must be -1 or nonnegative and fit int");
    if (!sparsify && (base_degree!=-1 || max_degree!=-1 || reactivate_limit!=-1))
        throw std::invalid_argument("Tesseract sparsification parameters require --tesseract-sparsify-errors");
    if (sparsify && (base_degree<=0 || (max_degree!=-1 && max_degree<base_degree)))
        throw std::invalid_argument("Tesseract sparsification requires base degree > 0 and max degree >= base (or -1)");
}

tesseract_decoder::TesseractConfig
TesseractOptions::make_config(const stim::DetectorErrorModel& dem) const
{
    validate("tesseract",0,nullptr);
    using Order=tesseract_decoder::DetectorOrder;
    std::vector<Order::Method> methods;
    if (order_bfs) methods.push_back(Order::Method::BFS);
    if (order_index) methods.push_back(Order::Method::Index);
    if (order_coordinate) methods.push_back(Order::Method::Coordinate);
    if (methods.empty()) methods.push_back(Order::Method::Index);
    tesseract_decoder::TesseractConfig config{.dem=dem};
    config.det_beam=beam==-1 ? tesseract_decoder::INF_DET_BEAM:static_cast<int>(beam);
    config.beam_climbing=beam_climbing;
    config.pqlimit=pqlimit==-1 ? std::numeric_limits<size_t>::max():static_cast<size_t>(pqlimit);
    config.no_revisit_dets=!allow_revisit;
    config.merge_errors=!no_merge;
    config.verbose=verbose;
    config.det_penalty=det_penalty;
    config.detector_orders.clear();
    config.detector_orders.reserve(num_orders);
    for (int64_t j=0;j<num_orders;j++)
        config.detector_orders.emplace_back(methods[j%methods.size()],static_cast<uint64_t>(order_seed)+j);
    config.sparsify_errors=sparsify;
    config.sparsify_base_degree=static_cast<int>(base_degree);
    config.sparsify_max_degree=static_cast<int>(max_degree);
    config.sparsify_reactivate_limit=static_cast<int>(reactivate_limit);
    return config;
}

void TesseractOptions::print(std::ostream& out) const
{
    out<<" beam="<<beam<<" beam_climbing="<<beam_climbing<<" pqlimit="<<pqlimit
       <<" detector_orders="<<num_orders<<" order_methods=";
    bool separator=false;
    auto method=[&](const char* name) {if (separator) out<<',';out<<name;separator=true;};
    if (order_bfs) method("bfs");
    if (order_index || (!order_bfs && !order_coordinate)) method("index");
    if (order_coordinate) method("coordinate");
    out<<" order_seed="<<order_seed<<" no_revisit_dets="<<!allow_revisit<<" merge_errors="<<!no_merge
       <<" det_penalty="<<det_penalty<<" verbose="<<verbose<<" sparsify_errors="<<sparsify
       <<" sparsify_base_degree="<<base_degree<<" sparsify_max_degree="<<max_degree
       <<" sparsify_reactivate_limit="<<reactivate_limit;
}
