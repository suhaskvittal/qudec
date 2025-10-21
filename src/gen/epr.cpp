/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 * */

#include "gen/epr.h"
#include "gen/scheduling.h"

namespace gen
{

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

stim::Circuit
epr_generation(const EPR_GEN_CONFIG& config, size_t rounds, size_t distance)
{
    // merged surface code has asymmetric distance for ZZ measurement
    SC_SCHEDULE_INFO sc1(distance, distance, 0);
    SC_SCHEDULE_INFO sc2(distance, distance, 2*distance*distance-1);

    // we need to pin extra ancilla qubits for the ZZ measurement
    // these pins belong to 
    std::vector<PINNED_QUBIT> pins;
    pins.reserve(distance);
    stim_qubit_type pinned_qubit{4*distance*distance-1};  // physical qubit id

    // pins are on last row of first code, and first row of second code
    for (size_t c = 0; c < distance; c++)
    {
        pins.push_back({pinned_qubit, 1, sc2.check_matrix.at(0).at(c)});
        pinned_qubit++;
    }

    stim::Circuit circuit;


}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

}   // namespace gen   
