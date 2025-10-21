
#include "gen/scheduling.h"

namespace gen
{

////////////////////////////////////////////////////
////////////////////////////////////////////////////

std::vector<stim_qubit_type>
surface_code_cx_order(stim_qubit_type nw, stim_qubit_type ne, stim_qubit_type sw, stim_qubit_type se, bool is_x_check)
{
    if (is_x_check) 
        return {nw, ne, sw, se};
    else
        return {nw, sw, ne, se};
}

SC_SCHEDULE_INFO::SC_SCHEDULE_INFO(size_t dx, size_t dy, stim_qubit_type start_qubit)
    :data_qubits(dx*dy),
    x_obs(dz),
    z_obs(dx),
    check_matrix(dx+1, std::vector<stim_qubit_type>(dy+1, NO_QUBIT))  // +1 on both dimensions due to boundary
{
    // initialize data structures describing syndrome extraction circuit:
    std::iota(data_qubits.begin(), data_qubits.end(), start_qubit);

    x_check_qubits.reserve((dx*dy-1)/2);
    z_check_qubits.reserve((dx-dy-1)/2);

    // start with bulk checks first:
    stim_qubit_type check_qubit{start_qubit + dx*dy};
    for (size_t r = 0; r < dx-1; r++)
    {
        for (size_t c = 0; c < dy-1; c++)
        {
            bool is_x_check = ((r+c) & 1) == 0;

            stim_qubit_type nw = start_qubit + r*distance + c;
            stim_qubit_type ne = start_qubit + r*distance + c+1;
            stim_qubit_type sw = start_qubit + (r+1)*distance + c;
            stim_qubit_type se = start_qubit + (r+1)*distance + c+1;

            std::vector<stim_qubit_type> cx_order = surface_code_cx_order(nw, ne, sw, se, is_x_check);
            if (is_x_check)
                x_check_qubits.push_back(check_qubit);
            else
                z_check_qubits.push_back(check_qubit);

            check_cx_order[check_qubit] = cx_order;
            check_matrix[r+1][c+1] = check_qubit;
            check_qubit++;
        }
    }

    // now handle boundary checks:
    // Z boundary checks:
    for (size_t r = 0; r < dx-1; r++)
    {
        std::vector<stim_qubit_type> check_order;
        if (r & 1)  // right boundary
        {
            stim_qubit_type q1 = start_qubit + (r+1)*dy - 1,
                            q2 = start_qubit + (r+2)*dy - 1;
            check_order = surface_code_cx_order(q1, NO_QUBIT, q2, NO_QUBIT, false);
        }
        else  // left boundary
        {
            stim_qubit_type q1 = start_qubit + r*dy,
                            q2 = start_qubit + (r+1)*dy;
            check_order = surface_code_cx_order(NO_QUBIT, q1, NO_QUBIT, q2, false);
        }

        z_check_qubits.push_back(check_qubit);
        check_cx_order[check_qubit] = check_order;
        
        size_t c = (r & 1) ? dy : 0;
        check_matrix[r+1][c] = check_qubit;

        check_qubit++;
    }

    // X boundary checks:
    for (size_t c = 0; c < dy-1; c++)
    {
        std::vector<stim_qubit_type> check_order;
        if (c & 1)  // upper boundary
        {
            stim_qubit_type q1 = start_qubit + c,
                            q2 = start_qubit + c+1;
            check_order = surface_code_check_order(NO_QUBIT, NO_QUBIT, q1, q2, true);
        }
        else  // bottom boundary
        {
            stim_qubit_type q1 = start_qubit + (dx-1)*dy + c,
                            q2 = start_qubit + (dx-1)*dy + c+1;
            check_order = surface_code_check_order(q1, q2, NO_QUBIT, NO_QUBIT, true);
        }

        x_check_qubits.push_back(check_qubit);
        check_cx_order[check_qubit] = check_order;

        size_t r = (c & 1) ? 0 : dx;
        check_matrix[r][c+1] = check_qubit;

        check_qubit++;
    }

    x_check_set = std::unordered_set<stim_qubit_type>(x_check_qubits.begin(), x_check_qubits.end());

    // x observable is a column of the code (choose leftmost)
    for (size_t i = 0; i < dx; i++)
        x_obs[i] = i*dy;
    // z observable is a row of the code (choose top)
    for (size_t i = 0; i < dy; i++)
        z_obs[i] = i;

#if defined(GEN_VERIFY_CX_SCHEDULE)
    // make sure no qubit is scheduled at the same timestep:
    for (size_t t = 0; t < 4; t++)
    {
        std::unordered_set<stim_qubit_type> qubits_in_this_timestep;
        for (const auto& [__unused_check, cx_order] : check_cx_order)
        {
            stim_qubit_type q = cx_order.at(t);
            if (q == NO_QUBIT)
                continue;
            if (qubits_in_this_timestep.count(q))
                throw std::runtime_error("qubit " + std::to_string(q) + " is scheduled at the same timestep twice");
            qubits_in_this_timestep.insert(q);
        }
    }
#endif
}

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////

}   // namespace gen
