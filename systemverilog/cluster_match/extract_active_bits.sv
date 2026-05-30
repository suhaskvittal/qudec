// author:  Suhas Vittal
// date:    28 May 2026

`include "cluster_match/common.svh"

////////////////////////////////////////////
////////////////////////////////////////////

module EXTRACT_ACTIVE_BITS
#(
    localparam int W         = `EXTRACT_ACTIVE_CHUNK_WIDTH,
    localparam int OUT_COUNT = `EXTRACT_ACTIVE_OUT_COUNT
)
(
    // Every cycle, this module will output at most `OUT_COUNT` detectors.
    // These are not guaranteed to be the first set bit, second set, etc.
    // as these detectors are extracted from different chunks. But eventually,
    // this module will complete extracting all bits.
    //
    // If a detector is invalid, the `d_valid` is set.
    output logic[`D_COUNT_PER_ROUND_BIT_WIDTH-1:0] d[OUT_COUNT-1:0],
    output logic[OUT_COUNT-1:0] d_valid,
    // `syndrome` and `syndrome_valid` indicate whether a syndrome should
    // be read, and what is the last valid bit in the syndrome.
    input logic[`D_COUNT_PER_ROUND-1:0] syndrome,
    input logic[`D_COUNT_PER_ROUND_BIT_WIDTH-1:0] syndrome_valid_last,
    // other signals
    input logic clk,
    input logic reset_n
);
    typedef logic[W-1:0] mask_type;

    mask_type mask[OUT_COUNT-1:0];
    logic[W-1:0] chunks[OUT_COUNT-1:0];
    logic[OUT_COUNT-1:0] d_penc_valid;

    // `d` is found using the priority encoder. We pass in `chunks[i]`
    // filtered out using `mask[i]`
    genvar ii;
    generate
        for (ii = 0; ii < OUT_COUNT; ii++)
        begin : gen_priority_enc
            PRIORITY_ENCODER #(.W(W)) enc( .idx(d[ii]),
                                            .valid(d_penc_valid[ii]),
                                            .chunk(chunks[ii] & ~mask[ii]) );
        end
    endgenerate

    // `d_valid` and `chunks` signals driven here:
    always_comb
    begin
        for (int i = 0; i < OUT_COUNT; i++)
        begin
            // `d_valid` is just `d_penc_valid` but with the check that the
            // detector id is beneath `syndrome_valid_last`
            d_valid[i] = d_penc_valid[i] & (d[i] < syndrome_valid_last);
            // `chunks` is just the syndrome divided accordingly:
            chunks[i] = syndrome[i*W +: W];
        end
    end

    always_ff @(posedge clk)
    begin
        // only reset `mask` on `reset_n`
        if (!reset_n)
        begin
            for (int i = 0; i < OUT_COUNT; i++)
                mask[i] <= '0;
        end

        // update `mask[i]` by setting bit `d[i]`
        for (int i = 0; i < OUT_COUNT; i++)
            mask[i][d[i]] <= 1'b1;
    end
endmodule

////////////////////////////////////////////
////////////////////////////////////////////

// LSB priority encoder
module PRIORITY_ENCODER
#(
    parameter int W = 32,
    localparam int LG_W = $clog2(W)
)
(
    output logic[LG_W-1:0] idx,
    output logic           valid,
    input  logic[W-1:0]    chunk
);
    always_comb
    begin
        idx = '0;
        valid = 1'b0;
        for (int i = 0; i < W; i++)
        begin
            if (chunk[i])
            begin
                idx = i;
                valid = 1'b1;
                break;
            end
        end
    end
endmodule

////////////////////////////////////////////
////////////////////////////////////////////
