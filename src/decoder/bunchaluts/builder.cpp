/*
 *  author: Suhas Vittal
 *  date:   21 September 2026
 * */

#include "decoder/bunchaluts/builder.h"

namespace decoder
{
namespace bal
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace
{

template <class GraphType> 
struct LUTBuildMeta
{
    GraphType  gr;
    CoordTable coord_map;
    std::vector<LUT> lut_levels_; // 0-indexed, so level 1 = idx 0
};

/*
 * Maps geometric coordinates to detectors. Note that the hypergraph
 * object should have the reverse mapping for each detector.
 * */
using CoordTable = std::unordered_map<hg::coord_type, hg::id_type>;

/*
 * Builds a `CoordTable` object given the DEM by simply reading
 * the arguments of each DEM instruction.
 * */
CoordTable _build_coord_table(const stim::DetectorErrorModel&);

/*
 * Computes the global frame flip information for `LUTEntry` using
 * `LUTKey`, coordinate information in `CoordTable`, and other information
 * in `G` which is a hypergraph.
 * */
template <class G>
void _elucidate_lut_entry_frame_flips(const LUTBuildMeta<G>&, const LUTKey&, LUTEntry&);

/*
 * Builds an LUT for the given error chain size (`k`). The LUT is built by
 * using the previous level's LUT as a jumping-off point. For each entry in
 * the previous level's LUT, we try and create new syndromes and their
 * information to build a new LUT.
 * */
template <class G>
LUT _build_lut_level(LUTBuildMeta<G>&, size_t k);

} // anon

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace
{

} // anon

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace bal
} // namespace decoder
