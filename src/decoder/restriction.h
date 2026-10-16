/*
 *  author: Suhas Vittal
 *  date:   16 October 2026
 * */

#ifndef DEC_RESTRICTION_h
#define DEC_RESTRICTION_h

#include "hypergraph.h"
#include "decoder/types.h"

#include <stim.h>

#include <memory>
#include <vector>

namespace dec
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace restriction
{

constexpr Color C_RED{0}, C_GREEN{1}, G_BLUE{2};

using det_id_type = hg::id_type;

struct v_data_type
{
    Color color;
    det_id_type canonical_detector;
    bool is_boundary;
    bool is_flag;

    /*
     * Corresponding detector IDs in the restricted lattices
     * */
    det_id_type rg_id,
                rb_id;
};

struct e_data_type
{
    double pr;
};

struct Match
{
    det_id_type x,
                y;
    Color c;
};

struct ConnectedComponent
{
    std::vector<det_id_type> edges;
    Color color;
};

/*
 * Used in connected components algorithm:
 * */
namespace dfs
{
struct v_data_type { det_id_type d; };
struct e_data_type { };
using graph_type = Hypergraph<v_data_type, e_data_type, 2>;
}


constexpr det_id_type 
boundary_id(Color c) 
{ 
    constexpr det_id_type BASE{0xdeadbeef};
    return BASE ^ (det_id_type{c} << 24);
}

constexpr bool
is_boundary_id(det_id_type d)
{
    return d == boundary_id(C_RED) || d == boundary_id(C_BLUE) || d == boundary_id(C_GREEN);
}

} // namespace restriction

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*
 * Color code decoder:
 * */

template <class Matcher>
class RestrictionDecoder
{
public:
    using det_id_type = restriction::det_id_type;
    using hg_type = HyperGraph<restriction::v_data_type, restriction::e_data_type, 3>;
    using hg_ptr = std::unique_ptr<hg_type>;

    using ma_ptr = std::unique_ptr<Matcher>;

    const size_t num_detectors,
                 num_observables;
    const stim::DetectorErrorModel& dem;

private:
    /*
     * Represent the overall DEM as a hypergraph, which helps in lifting and such
     * */
    hg_ptr hg_; 
    
    /*
     * And we need an instance of PyMatching for the RG and RB restricted lattices
     * */
    ma_ptr rg_matcher_,
           rb_matcher_;

    /*
     * Data structures used in each call:
     * */
    mutable std::vector<det_id_type> detectors_;
    mutable std::unordered_map<size_t, lift_orbit> orbits_;
public:
    RestrictionDecoder(const stim::DetectorErrorModel&);

    Result decode(SyndromeRef, ObsRef) const;
private:
    void matching(ma_ptr&, Color ignored_color) const;
    void build_connected_components(int8_t lattice_id, std::vector<restriction::Match>);
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace dec

#endif // DEC_RESTRICTION_h
