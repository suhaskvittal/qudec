// author: Suhas Vittal
// date: 19 June 2026
//
// Implementation of logic to filter out isolated edges from
// syndrome

`include "cluster_match/common.svh"

////////////////////////////////////////////
////////////////////////////////////////////

module FILTER
(
    output logic[`D_MAX-1:0]    s_mask,
    output logic[`LG_D_MAX-1:0] d_processed_cnt,

    input logic         d_v,
    input detector_type d,
    input logic         is_degree_one,
    input detector_type a,

    input clk,
    input reset_n
);
    // Whenever a detector (`d`) and its associated data comes in,
    // if its degree is one, we should buffer `a`. That way, if `a` comes in
    // the future, then we can update `s_mask`
    //
    // We can simplify this design by simply assuming data comes in-order,
    // so we can just increment a counter whenever we get new data.

    // In the worst case, every detection event has degree one. Since we will
    // filter on a match, we will have at most `D_MAX/2 installs.
    localparam int BUFFER_SIZE = `D_MAX / 2;
    localparam int LG_BUFFER_SIZE = $clog2(BUFFER_SIZE);

    typedef struct packed
    {
        logic                v;
        logic[`LG_D_MAX-1:0] src_idx;
        detector_type        dst_match;
    } lookup_buf_entry;

    // This is a mirror of the `s_mask` output signal:
    logic[`D_MAX-1:0] _s_mask;

    // counter + lookup buffer
    logic[`LG_D_MAX-1:0] ctr;
    lookup_buf_entry[BUFFER_SIZE-1:0] lu_buffer;
    logic[LG_BUFFER_SIZE-1:0] buf_ptr;

    logic                     buf_match;
    logic[LG_BUFFER_SIZE-1:0] buf_match_idx;

    // Drive output signals:
    assign s_mask = _s_mask;
    assign d_processed_cnt = ctr;

    always_ff @(posedge clk)
    begin
        if (!reset_n)
        begin
            ctr <= '0;
            buf_ptr <= '0;
            for (int i = 0; i < BUFFER_SIZE; i++)
                lu_buffer[i].v <= 1'b0;
            _s_mask <= '1;
        end
        else
        begin
            // increment counter
            if (d_v)
            begin
                ctr <= ctr+1;

                // only do anything if `is_degree_one` is active high (otherwise
                // there is nothing to filter out).
                if (is_degree_one)
                begin
                    if (buf_match)
                    begin
                        // update mask and invalidate buffer entry
                        _s_mask[lu_buffer[buf_match_idx].src_idx] <= 1'b0;
                        _s_mask[ctr] <= 1'b0;
                        lu_buffer[buf_match_idx].v <= 1'b0;
                    end
                    else
                    begin
                        // install into next buffer entry
                        lu_buffer[buf_ptr].v <= 1'b1;
                        lu_buffer[buf_ptr].src_idx <= ctr;
                        lu_buffer[buf_ptr].dst_match <= a;
                        buf_ptr <= buf_ptr+1;
                    end
                end
            end
        end
    end

    always_comb
    begin
        // drive `buf_match` and `buf_match_idx`
        buf_match = 1'b0;
        buf_match_idx = '0;
        for (int i = 0; i < BUFFER_SIZE; i++)
        begin
            if (lu_buffer[i].v && lu_buffer[i].dst_match == d)
            begin
                buf_match = 1'b1;
                buf_match_idx = i;
                break;
            end
        end
    end
endmodule;

////////////////////////////////////////////
////////////////////////////////////////////
