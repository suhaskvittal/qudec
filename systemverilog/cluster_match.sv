// author:  Suhas Vittal
// date:    28 May 2026
//
// Top level module implementation of Cluster-Match decoder

`include "cluster_match/astrea.sv"
`include "cluster_match/common.svh"

////////////////////////////////////////////
////////////////////////////////////////////

module CLUSTER_MATCH
#(
    parameter int MAX_CLUSTERS = 10,
    parameter int D_ROUND_BUFFER_SIZE = 16
(
    // `frame_flips` is the output data. `obs_assert` indicates `frame_flips`
    // is being driven by the decoder.
    output obs_type frame_flips,
    output logic obs_assert,
    // `syndrome` correpsonds to the data send in after measurements
    // in a syndrome extraction round. This is generally less than
    // `D_COUNT`. `syndrome_valid` points to one bit past
    // the last syndrome bit that is valid.
    input logic[`D_COUNT_PER_ROUND-1:0] syndrome,
    input logic[`D_COUNT_PER_ROUND_BIT_WIDTH-1:0] syndrome_valid,
    // `syndrome_assert` is active-high and indicates that `syndrome` and
    // `syndrome_valid` can be read.
    input logic syndrome_assert,
    // `clk` and `reset_n` are self-explanatory. `reset_n` is active-low and
    // must be asserted on the positive edge of `clk`
    input logic clk,
    input logic reset_n
);
    // adjacency matrix memory:
    adj_list_type adj_mem[`D_COUNT-1:0];

    // We will process syndromes as they come in per round. We need
    // to buffer the identified detectors before giving the `FILTER`
    // submodule.

    typedef logic[`D_COUNT_PER_ROUND_BIT_WIDTH-1:0] d_buf_entry;
    d_buf_type bf_prev_round[D_ROUND_BUFFER_SIZE-1:0];
    d_buf_type bf_curr_round[D_ROUND_BUFFER_SIZE-1:0];
    d_buf_type bf_next_round[D_ROUND_BUFFER_SIZE-1:0];


    // `cl_frame_flips` are the frame flips specific to each
    // cluster (matching subproblem). `cl_edges` are the
    // MWPM edges computed using a distance finding algorithm,
    // and `cl_severity` indicates which version of Astrea to use.
    obs_type    cl_frame_flips[MAX_CLUSTERS-1:0];
    edge_type   cl_edges[MAX_CLUSTERS-1:0][14:0];
    logic[1:0]  cl_severity[MAX_CLUSTERS-1:0];
    logic       cl_valid[MAX_CLUSTERS-1:0];

    // initial block for instantiating structures + begin validation
    initial
    begin
        $readmemh("adjacency_matrix.mem", adj_mem);
    end

    // `frame_flips` will be driven by XOR of all `cl_frame_flips`:
    always_comb
    begin
        frame_flips = '0;
        for (int i = 0; i < MAX_CLUSTERS; i++)
        begin
            if (cl_valid[i])
                frame_clips ^= cl_frame_flips[i];
        end
    end
    
    // `cl_frame_flips` will be driven by Astrea (one instance per
    // cluster)
    genvar i;
    generate
        for (i = 0; i < MAX_CLUSTERS; i++)
        begin : gen_astrea
            ASTREA_MUX astrea( .frame_flips(cl_frame_flips[i]),
                                .edges(cl_edges[i]),
                                .severity(cl_sev[i]) );
        end
    endgenerate

endmodule

////////////////////////////////////////////
////////////////////////////////////////////

module ASTREA_MUX
(
    output obs_type frame_flips,
    input edge_type edges[14:0],
    input logic[1:0] severity
);  
    obs_type hw2_f;
    obs_type hw4_f;
    obs_type hw6_f;

    always_comb
    begin
        if (severity == 0)  // HW = 2
            frame_flips = hw2_f;
        else if (severity == 1) // HW = 4
            frame_flips = hw4_f;
        else
            frame_flips = hw6_f;
    end

    // `hw2_f` is simple -- only one edge is valid:
    // `hw4_f` and `hw6_f` require running Astrea sub-modules:
    assign hw2_f = edges[0].f;
    ASTREA_HW4 a4( .frame_flips(hw4_f), .edges(edges[5:0]) );
    ASTREA_HW6 a6( .frame_flips(hw6_f), .edges(edges[14:0]) );
endmodule

////////////////////////////////////////////
////////////////////////////////////////////
