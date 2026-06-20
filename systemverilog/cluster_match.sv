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
    input detector_type  d[`D_MAX-1:0],
    input logic[`D_MAX-1:0] d_valid,

    // `clk` and `reset_n` are self-explanatory. `reset_n` is active-low and
    // must be asserted on the positive edge of `clk`
    input logic clk,
    input logic reset_n
);
    // `mem` is the data used for decoding
    logic[] mem[`D_COUNT-1:0]; 

endmodule

////////////////////////////////////////////
////////////////////////////////////////////
