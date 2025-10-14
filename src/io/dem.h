/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 * */

#ifndef IO_DEM_h
#define IO_DEM_h

#include "decoding_graph.h"

#include <stim/dem/detector_error_model.h>

namespace io
{

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

constexpr size_t DEM_COLOR_COORD_IDX{0};
constexpr size_t DEM_FLAG_COORD_IDX{1}; 

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

struct DEM_BLOCK_INFO
{
    constexpr static size_t MAX_COORD{6};

    // this structure fits in 8 + 4*6 = 32 bytes as of now (half a cacheline)
    int64_t                      id_shift{0};
    std::array<float, MAX_COORD> coord_shift{};
};

struct DEM_READ_RESULT
{
    // first part is detectors, second part is the actual data
    using detector_type = std::pair<int64_t, DETECTOR_DATA>;
    using error_type = std::pair<std::vector<int64_t>, DECODER_ERROR_DATA>;

    std::vector<detector_type> detectors;
    std::vector<error_type>    errors;
};

DEM_READ_RESULT read_dem_block(const stim::DetectorErrorModel& dem);

void read_dem_block_helper(const stim::DetectorErrorModel& dem, DEM_READ_RESULT&, DEM_BLOCK_INFO&);

std::vector<DEM_READ_RESULT::detector_type> read_detector_decl(const stim::DemInstruction&, DEM_BLOCK_INFO& info);
std::vector<DEM_READ_RESULT::error_type>    read_dem_error(const stim::DemInstruction&, DEM_BLOCK_INFO& info);

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

}   // namespace io

#endif
