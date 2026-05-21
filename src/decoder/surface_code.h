/*
 *  author: Suhas Vittal
 *  date:   16 May 2026
 * */

#ifndef DECODER_SURFACE_CODE_h
#define DECODER_SURFACE_CODE_h

#include "decoder/common.h"
#include "decoder/logger.h"

#include <stim.h>
#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"

namespace decoder
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

class PYMATCHING
{
public:
    const size_t num_detectors;
    const size_t num_observables;
private:
    pm::Mwpm    mwpm_;
public:
    PYMATCHING(const stim::DetectorErrorModel&);

    result_type decode(syndrome_ref, LOGGER&);
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace decoder

#endif // DECODER_SURFACE_CODE_h
