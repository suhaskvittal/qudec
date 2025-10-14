/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 * */

#include "io/dem.h"

namespace io
{

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

DEM_READ_RESULT
read_dem_block(const stim::DetectorErrorModel& dem)
{
    DEM_BLOCK_INFO info;
    DEM_READ_RESULT result;
    read_dem_block_helper(dem, result, info);
    return result;
}

void
read_dem_block_helper(const stim::DetectorErrorModel& dem, DEM_READ_RESULT& result, DEM_BLOCK_INFO& info)
{
    for (const auto& inst : dem.instructions)
    {
        if (inst.type == stim::DemInstructionType::DEM_ERROR)
        {
            auto errors = read_dem_error(inst, info);
            result.errors.insert(result.errors.end(), errors.begin(), errors.end());
        }
        else if (inst.type == stim::DemInstructionType::DEM_DETECTOR)
        {
            auto decl = read_detector_decl(inst, info);
            result.detectors.insert(result.detectors.end(), decl.begin(), decl.end());
        }
        else if (inst.type == stim::DemInstructionType::DEM_REPEAT_BLOCK)
        {
            const auto& b = inst.repeat_block_body(dem);
            size_t num_reps = inst.repeat_block_rep_count();
            for (size_t i = 0; i < num_reps; ++i)
                read_dem_block_helper(b, result, info);
        }
        else if (inst.type == stim::DemInstructionType::DEM_SHIFT_DETECTORS)
        {
            for (size_t i = 0; i < inst.arg_data.size(); ++i)
                info.coord_shift[i] += inst.arg_data[i];
            info.id_shift += inst.target_data[0].data;
        }
    }
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

std::vector<DEM_READ_RESULT::detector_type>
read_detector_decl(const stim::DemInstruction& inst, DEM_BLOCK_INFO& info)
{
    // copy coordinates
    std::array<float, DEM_BLOCK_INFO::MAX_COORD> coords(info.coord_shift);
    for (size_t i = 0; i < inst.arg_data.size(); ++i)
        coords[i] += inst.arg_data[i];

    // get information from coordinates
    int color_id = std::round(coords[DEM_COLOR_COORD_IDX]);
    bool is_flag = std::round(coords[DEM_FLAG_COORD_IDX]) > 0.0f;

    // compute id:
    DETECTOR_DATA data
    {
        static_cast<DETECTOR_DATA::COLOR>(color_id),
        is_flag
    };

    std::vector<DEM_READ_RESULT::detector_type> decl(inst.target_data.size());
    std::transform(inst.target_data.begin(), inst.target_data.end(), decl.begin(),
                [&data, s=info.id_shift] (const stim::DemTarget& t) 
                { 
                    return std::make_pair(t.data + s, data);
                });

    return decl;
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

std::vector<DEM_READ_RESULT::error_type>
read_dem_error(const stim::DemInstruction& inst, DEM_BLOCK_INFO& info)
{
    std::vector<DEM_READ_RESULT::error_type> errors;

    double error_prob = inst.arg_data[0];
    inst.for_separated_targets(
        [&errors, &info, error_prob] (std::span<const stim::DemTarget> tg)
        {
            std::vector<int64_t> detectors;
            DECODER_ERROR_DATA ed{error_prob};

            // identify detectors and observables:
            for (const auto& t : tg)
            {
                if (t.is_observable_id())
                    ed.flipped_observables.insert(t.data);
                else
                    detectors.push_back(t.data + info.id_shift);
            }

            errors.emplace_back(detectors, ed);
        });

    return errors;
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

}   // namespace io

