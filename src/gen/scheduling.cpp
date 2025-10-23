/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 * */

#include "gen/scheduling.h"
#include "gen/utils.h"

#include <iomanip>
#include <iostream>
#include <numeric>

namespace gen
{

////////////////////////////////////////////////////
////////////////////////////////////////////////////

SC_SCHEDULE_INFO::SC_SCHEDULE_INFO(size_t _dx, size_t _dz, bool _is_dual)
    :dx(_dx),
    dz(_dz),
    is_dual(_is_dual),
    data_qubits(dx*dz),
    x_obs(dz),
    z_obs(dx),
    check_matrix(dz+1, std::vector<stim_qubit_type>(dx+1, NO_QUBIT))  // +1 on both dimensions due to boundary
{
    // initialize data structures describing syndrome extraction circuit:
    std::iota(data_qubits.begin(), data_qubits.end(), 0);

    const size_t num_checks = dx*dz-1;
    x_check_qubits.reserve(num_checks/2);
    z_check_qubits.reserve(num_checks/2);

    // start with bulk checks first:
    stim_qubit_type check_qubit{dx*dz};
    check_qubit = init_bulk_checks(check_qubit);
    check_qubit = init_left_right_boundary_checks(check_qubit);
    check_qubit = init_top_bottom_boundary_checks(check_qubit);

    x_check_set = std::unordered_set<stim_qubit_type>(x_check_qubits.begin(), x_check_qubits.end());

    // x observable is a column of the code (choose leftmost)
    for (size_t i = 0; i < dz; i++)
        x_obs[i] = i*dx;
    // z observable is a row of the code (choose top)
    for (size_t i = 0; i < dx; i++)
        z_obs[i] = i;

    util::validate_check_cx_order(check_cx_order);

    // create coordinates:
    for (size_t r = 0; r < dz; r++)
    {
        for (size_t c = 0; c < dx; c++)
            qubit_coords.insert({r*dx + c, std::make_pair(c+0.5, r+0.5)});
    }

    for (size_t r = 0; r <= dz; r++)
    {
        for (size_t c = 0; c <= dx; c++)
        {
            if (check_matrix[r][c] == NO_QUBIT)
                continue;
            qubit_coords.insert({check_matrix[r][c], std::make_pair(c, r)});
        }
    }
}

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////

stim_qubit_type
SC_SCHEDULE_INFO::init_bulk_checks(stim_qubit_type check_qubit)
{
    for (size_t r = 0; r < dz-1; r++)
    {
        for (size_t c = 0; c < dx-1; c++)
        {
            bool is_x_check = (((r+c) & 1) == 0) ^ is_dual;

            stim_qubit_type nw = r*dx + c;
            stim_qubit_type ne = r*dx + c+1;
            stim_qubit_type sw = (r+1)*dx + c;
            stim_qubit_type se = (r+1)*dx + c+1;

            std::vector<stim_qubit_type> cx_order = surface_code_cx_order(nw, ne, sw, se, is_x_check);
            add_check_decl(check_qubit, r+1, c+1, std::move(cx_order), is_x_check);
            check_qubit++;
        }
    }
    return check_qubit;
}

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////

stim_qubit_type
SC_SCHEDULE_INFO::init_left_right_boundary_checks(stim_qubit_type check_qubit)
{
    const bool checks_are_x_type = is_dual;

    if (dz & 1)
    {
        for (size_t r = 0; r < dz-1; r++)
        {
            std::vector<stim_qubit_type> check_order;
            if (r & 1)  // right boundary
            {
                stim_qubit_type q1 = (r+1)*dx - 1,
                                q2 = (r+2)*dx - 1;
                check_order = surface_code_cx_order(q1, NO_QUBIT, q2, NO_QUBIT, checks_are_x_type);
            }
            else  // left boundary
            {
                stim_qubit_type q1 = r*dx,
                                q2 = (r+1)*dx;
                check_order = surface_code_cx_order(NO_QUBIT, q1, NO_QUBIT, q2, checks_are_x_type);
            }
            add_check_decl(check_qubit, r+1, (r & 1) ? dx : 0, std::move(check_order), checks_are_x_type);
            check_qubit++;
        }
    }
    else
    {
        for (size_t r = 0; r < dz-1; r += 2)
        {
            // two checks in a row:
            stim_qubit_type q1 = r*dx,
                            q2 = (r+1)*dx,
                            q3 = (r+1)*dx - 1,
                            q4 = (r+2)*dx - 1;

            auto check_order_1 = surface_code_cx_order(NO_QUBIT, q1, NO_QUBIT, q2, checks_are_x_type);
            auto check_order_2 = surface_code_cx_order(q3, NO_QUBIT, q4, NO_QUBIT, checks_are_x_type);
            add_check_decl(check_qubit, r+1, 0, std::move(check_order_1), checks_are_x_type);
            add_check_decl(check_qubit+1, r+1, dx, std::move(check_order_2), checks_are_x_type);
            check_qubit += 2;
        }
    }

    return check_qubit;
}

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////

stim_qubit_type
SC_SCHEDULE_INFO::init_top_bottom_boundary_checks(stim_qubit_type check_qubit)
{
    const bool dx_is_odd = (dx & 1);
    const bool checks_are_x_type = (dx_is_odd && !is_dual) || (!dx_is_odd && is_dual);

    if (dx_is_odd)
    {
        for (size_t c = 0; c < dx-1; c++)
        {
            std::vector<stim_qubit_type> check_order;
            if (c & 1)  // upper boundary
            {
                stim_qubit_type q1 = c,
                                q2 = c+1;
                check_order = surface_code_cx_order(NO_QUBIT, NO_QUBIT, q1, q2, checks_are_x_type);
            }
            else  // bottom boundary
            {
                stim_qubit_type q1 = (dz-1)*dx + c,
                                q2 = (dz-1)*dx + c+1;
                check_order = surface_code_cx_order(q1, q2, NO_QUBIT, NO_QUBIT, checks_are_x_type);
            }

            add_check_decl(check_qubit, (c & 1) ? 0 : dz, c+1, std::move(check_order), checks_are_x_type);
            check_qubit++;
        }
    }
    else
    {
        for (size_t c = 0; c < dx-1; c += 2)
        {
            stim_qubit_type q1 = c,
                            q2 = c+1,
                            q3 = (dz-1)*dx + c,
                            q4 = (dz-1)*dx + c+1;

            auto check_order_1 = surface_code_cx_order(NO_QUBIT, NO_QUBIT, q1, q2, checks_are_x_type);
            auto check_order_2 = surface_code_cx_order(q3, q4, NO_QUBIT, NO_QUBIT, checks_are_x_type);

            add_check_decl(check_qubit, 0, c+1, std::move(check_order_1), checks_are_x_type);
            add_check_decl(check_qubit+1, dz, c+1, std::move(check_order_2), checks_are_x_type);
            check_qubit += 2;
        }
    }
    return check_qubit;
}

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////

void
SC_SCHEDULE_INFO::add_check_decl(
        stim_qubit_type check_qubit, 
        size_t r, 
        size_t c, 
        std::vector<stim_qubit_type>&& cx_order, 
        bool is_x_check)
{
    check_cx_order[check_qubit] = std::move(cx_order);
    if (is_x_check)
        x_check_qubits.push_back(check_qubit);
    else
        z_check_qubits.push_back(check_qubit);
    check_matrix[r][c] = check_qubit;
}

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////

std::vector<stim_qubit_type>
surface_code_cx_order(stim_qubit_type nw, stim_qubit_type ne, stim_qubit_type sw, stim_qubit_type se, bool is_x_check)
{
    if (is_x_check) 
        return {nw, ne, sw, se};
    else
        return {nw, sw, ne, se};
}

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////

}   // namespace gen
