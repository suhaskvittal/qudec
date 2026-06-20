// author: Suhas Vittal
// date: 19 June 2026
//
// Implementation of logic to load neighbor data

`include "cluster_match/common.svh"

////////////////////////////////////////////
////////////////////////////////////////////

module INITIALIZE_NEIGHBORS
(
    output logic         requestor_out_v,
    output detector_type requestor_out,
    output adj_type      adj_list,
    output logic         is_degree_one,
    output detector_type a,

    input logic         requestor_in_v,
    input detector_type requestor_in,

    input syndrome_type syn,

    input adj_type mem[`D_COUNT-1:0],

    input clk,
    input reset_n
);
    localparam int STAGES = 3;

    typedef struct packed
    {
        logic                      v;
        detector_type              requestor;
        adj_type                   adj_list;
        logic[`ADJ_MAX_DEGREE-1:0] present_in_syndrome;
        logic[`LG_D_MAX:0]         degree;
        detector_type              a; 
    } stage_type;

    stage_type st[STAGES-1:0];

    // drive output signals
    assign requestor_out_v = st[STAGES-1].v;
    assign requestor_out = st[STAGES-1].requestor;
    assign adj_list = st[STAGES-1].adj_list;
    assign is_degree_one = (st[STAGES-1].degree == 1);
    assign a = st[STAGES-1].a;

    always_ff @(posedge clk)
    begin
        if (!reset_n)
        begin
            for (int i = 0; i < STAGES; i++)
                st[i].v <= 1'b0;
        end
        else
        begin
            // passthrough logic:
            st[0].v <= requestor_in_v;
            for (int i = 1; i < STAGES; i++)
            begin
                st[i].v <= st[i-1].v;
                st[i].requestor <= st[i-1].requestor;
                st[i].adj_list <= st[i-1].adj_list;
            end

            // Pipeline stage 0: read adjacency list from `mem`
            st[0].requestor <= requestor_in;
            st[0].adj_list <= mem[requestor_in];

            // Pipeline stage 1: fully associative lookup with `syn` to
            // check which neighbors are present in the syndrome
            for (int i = 0; i < `ADJ_MAX_DEGREE; i++)
            begin
                for (int j = 0; j < `D_MAX; j++)
                begin
                    st[1].present_in_syndrome[i] <= st[1].present_in_syndrome
                                                    | (st[0].adj_list[i].d == syn[j]);
                end
            end

            // Pipeline stage 2: get degree of `present_in_syndrome` and use
            // priority encoder to set `a`
            st[2].degree <= $countones(st[1].present_in_syndrome);
            st[2].a <= '0;
            for (int i = 0; i < `ADJ_MAX_DEGREE; i++)
            begin
                if (st[1].present_in_syndrome[i])
                begin
                    st[2].a <= st[1].adj_list[i];
                    break;
                end
            end
        end
    end
endmodule;

////////////////////////////////////////////
////////////////////////////////////////////
