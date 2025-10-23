/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 * */

#ifndef GEN_SCHEDULING_h
#define GEN_SCHEDULING_h

#include "gen.h"

#include <utility>

namespace gen
{

////////////////////////////////////////////////////
////////////////////////////////////////////////////

struct SC_SCHEDULE_INFO
{
    using qubit_coord_type = std::pair<double,double>;

    const size_t dx;
    const size_t dz;
    const bool   is_dual;

    std::vector<stim_qubit_type> data_qubits;
    std::vector<stim_qubit_type> x_check_qubits;
    std::vector<stim_qubit_type> z_check_qubits;
    std::unordered_map<stim_qubit_type, std::vector<stim_qubit_type>> check_cx_order;

    std::unordered_set<stim_qubit_type> x_check_set;

    std::vector<stim_qubit_type> x_obs;
    std::vector<stim_qubit_type> z_obs;

    // a useful structure for looking up check qubits by row + column. If there is no
    // entry at a given location, then the entry is `NO_QUBIT` (will occur on boundary).
    std::vector<std::vector<stim_qubit_type>> check_matrix;

    // utility for crumble setup:
    std::map<stim_qubit_type, qubit_coord_type> qubit_coords;

    // dual simply swaps the x and z check locations -- more useful for even distance codes
    SC_SCHEDULE_INFO(size_t dx, size_t dz, bool dual=false);
private:
    stim_qubit_type init_bulk_checks(stim_qubit_type start_qubit);
    stim_qubit_type init_left_right_boundary_checks(stim_qubit_type start_qubit);
    stim_qubit_type init_top_bottom_boundary_checks(stim_qubit_type start_qubit);

    void add_check_decl(stim_qubit_type, size_t r, size_t c, std::vector<stim_qubit_type>&& cx_order, bool is_x_check);
};

////////////////////////////////////////////////////
////////////////////////////////////////////////////

std::vector<stim_qubit_type> surface_code_cx_order(stim_qubit_type nw, 
                                                    stim_qubit_type ne, 
                                                    stim_qubit_type sw,
                                                    stim_qubit_type se,
                                                    bool is_x_check);

////////////////////////////////////////////////////
////////////////////////////////////////////////////

}   // namespace gen

#endif  // GEN_SCHEDULING_h
