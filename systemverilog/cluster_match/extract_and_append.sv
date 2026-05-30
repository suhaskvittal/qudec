// author:  Suhas Vittal
// date:    28 May 2026

`include "cluster_match/extract_active_bits.sv"

////////////////////////////////////////////
////////////////////////////////////////////

module EXTRACT_AND_APPEND
#(
    localparam int OUT = `EXTRACT_ACTIVE_OUT_COUNT,
    localparam int LG_OUT = $clog2(OUT)
(
    // `rbfd` is the data buffer of the round buffer.
    output d_rbf_entry rbfd[OUT-1:0],
    // we also want to pipe through `d` and `d_valid` which are driven by
    // `EXTRACT_ACTIVE_BITS` so we can access BRAM in the frontend while
    output d_rbf_entry    d[OUT-1:0],
    output logic[OUT-1:0] d_valid,
    // `ptr` points to the first unoccupied entry in `rbfd`. This module
    // assumes `rbfd` is contiguous in terms of occupancy.
    inout logic[`LG_RBF_SIZE-1:0] ptr,
    // syndrome data to drive `EXTRACT_ACTIVE_BITS`
    input logic[`D_COUNT_PER_ROUND-1:0] syndrome,
    input logic[`D_COUNT_PER_ROUND_BIT_WIDTH-1:0] syndrome_valid_last,
    input logic clk,
    input logic reset_n
);
    // assert that `ptr` is not larger than `RBF_SIZE`
    always_comb
    begin
        assert(ptr < `RBF_SIZE)
            else $error("OOB access to RBF: ptr = %0d", ptr);
    end

    logic[LG_OUT-1:0] d_offset[OUT-1:0];

    // Move `d` into `rbfd`. First determine the offsets for each
    // insertion:
    always_comb
    begin
        d_offset[0] = '0;
        for (int i = 1; i < OUT; i++)
            d_offset[i] = d_offset[0] + d_valid[i-1];
    end

    // Now drive `rbfd`
    always_ff @(posedge clk)
    begin
        for (int i = 0; i < OUT; i++)
        begin
            if (d_valid[i])
                rbfd[ptr+d_offset[i]] <= d[i];
        end

        ptr <= d_offset[OUT-1];
    end

    // drive `d` and `d_valid`
    EXTRACT_ACTIVE_BITS exta( .d(d),
                                .d_valid(d_valid),
                                .syndrome(syndrome),
                                .syndrome_valid_last(syndrome_valid_last),
                                .clk(clk),
                                .reset_n(reset_n) );
endmodule

////////////////////////////////////////////
////////////////////////////////////////////
