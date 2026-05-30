// author:  Suhas Vittal
// date:    28 May 2026
//
// Frontend of Cluster-Match decoder

`include "cluster_match/extract_and_append.sv"

////////////////////////////////////////////
////////////////////////////////////////////

module FRONTEND
(
    input logic clk,
    input logic reset_n
);
    // The frontend will maintain three round-level buffers.
    //  buf0 = round R-2
    //  buf1 = round R-1
    //  buf2 = round R
    // 
    // We can determine the number of active neighbors for any
    // detector in round R-1. We will speculatively start filtering
    // out neighbors. If a new (nonzero) syndrome is asserted, then
    // the filtered syndrome is discarded.

    d_rbuf_entry rbf0;  // non-speculative always
    d_rbuf_entry rbf1_ns;  // non-speculative `rbf1`
    d_rbuf_entry rbf2_ns;  // non-speculative `rbf2`
    d_rbuf_entry rbf1_s;  // speculative `rbf1`
    d_rbuf_entry rbf2_s;  // speculative `rbf2`
    
    EXTRACT_AND_APPEND exta( .rbfd(rbf2_s.d),
                                .d(),
                                .d_valid(),
                                .ptr(rbf2_s.ptr),
                                .syndrome(),
                                .syndrome_valid_last(),
                                .clk(clk),
                                .reset_n(reset_n) );
endmodule

////////////////////////////////////////////
////////////////////////////////////////////
