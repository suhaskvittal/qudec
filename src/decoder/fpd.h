/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 * */

#ifndef DECODER_FPD_h
#define DECODER_FPD_h

#include "decoding_graph.h"
#include "decoder/common.h"

#include <limits>
#include <unordered_map>
#include <vector>

// Global debug configuration variable
extern bool GL_DEBUG_DECODER;

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

struct FPD_CONFIG
{
    size_t cache_chain_limit{3};
    bool   do_not_predecode_if_any_without_pref{true};
};

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

template <class LL_DECODER>
class FPD
{
public:
    constexpr static GRAPH_COMPONENT_ID NO_PREF{-1};

    struct error_info
    {
        int8_t             length;
        syndrome_type      flipped_obs;
    };

    struct preference_entry_type
    {
        GRAPH_COMPONENT_ID pref{NO_PREF};
        uint8_t            count{0};
        int8_t             length{std::numeric_limits<int8_t>::max()};
    };

    using ec_cache_array = std::unordered_map<GRAPH_COMPONENT_ID, error_info>;

    const FPD_CONFIG conf;
private:
    // cache error chains up-to `cache_chain_limit`
    std::vector<ec_cache_array> ec_cache;

    LL_DECODER ll_decoder;

    GRAPH_COMPONENT_ID boundary_index;
public:
    FPD(const stim::Circuit&, LL_DECODER&&, FPD_CONFIG={});

    DECODER_RESULT decode(std::vector<GRAPH_COMPONENT_ID>, std::ostream& debug_strm) const;
private:
    void init_ec_cache_array(SC_DECODING_GRAPH*, GRAPH_COMPONENT_ID);

    std::vector<preference_entry_type> compute_prefs(const std::vector<GRAPH_COMPONENT_ID>&) const;
};

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

#include "fpd.tpp"

#endif  // DECODER_FPD_h
