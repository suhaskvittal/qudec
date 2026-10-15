/*
 *  author: Suhas Vittal
 *  date    15 October 2026
 * */

#ifndef DECODER_PHANTOM_h
#define DECODER_PHANTOM_h

#include "decoder/types.h"
#include "decoder/surface_code.h"
#include "preprocessor/affinity.h"
#include "stats.h"

#include <stim.h>

#include <iosfwd>
#include <limits>
#include <memory>
#include <vector>

namespace dec
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

class Phantom
{
public:
    using det_id_type = size_t;
    using aff_type = pp::AffinityTable::value_type;

    /* `code_distance` and `d` are the same thing. */
    const size_t code_distance;
    const size_t d; 

    const stim::DetectorErrorModel& dem;
    const size_t num_detectors,
                 num_observables;

    /* Stats: */
private:
    std::unique_ptr<PyMatching> internal_matcher_; 
    pp::AffinityTable aff_;

    mutable std::vector<det_id_type> detectors_;
    mutable std::vector<bool> is_matched_;
    mutable std::vector<aff_type> pairwise_affinity_;
public:
    Phantom(const stim::DetectorErrorModel&, size_t code_distance);

    Result decode(SyndromeRef, ObsRef);

    void print_stats(std::ostream&) const;
    void mpi_accumulate();

    size_t remaining() const { return detectors_.size() - std::count(is_matched_.begin(), is_matched_.end(), true); }
private:
    Result filter(size_t w);
    std::vector<std::vector<det_id_type>> build_clusters(size_t max_aff);
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace dec

#endif
