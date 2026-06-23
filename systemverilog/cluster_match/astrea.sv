//
//  author: Suhas Vittal
//  date:   28 May 2026
//
// Implementations of Astrea decoder.
//  Expected ordering of input edges:
//      01, 02, 03, ..., 0n, 12, 13, 14, ..., 1n
//  Where "xy" is an edge between detector x and detector y.
//  You will notice that `edge_type` does not contain
//  any data about the detectors connected by the edge.
//  This is unimportant for calculating the frame changes
//  via matching. We only need the weight and frame change.
// 
// All implementations of Astrea will largely follow
// the same input/output semantics.

`include "cluster_match/common.svh"

/* verilator lint_off MULTITOP */

////////////////////////////////////////////
////////////////////////////////////////////

// Overall Astrea Decoder:
module ASTREA
(
    output weight_type    m_weight,
    output obs_type       frame_flips,
    input edge_type[14:0] edges,
    input logic[14:0]     edges_v
);
    // `_m_weight` and `_frame_flips` mirror the outputs 
    // and contains the results of running each variant of Astrea.
    weight_type[2:0] _m_weight;
    obs_type[2:0] _frame_flips;

    // depending on the number of edges, we will use the outcome
    // of a different version of Astrea
    logic[3:0] edge_count;

    // HW2:
    assign _m_weight[0] = edges[0].w;
    assign _frame_flips[0] = edges[0].f;
    // HW4:
    ASTREA_HW4 a4(.m_weight(_m_weight[1]), .frame_flips(_frame_flips[1]), .edges(edges[5:0]));
    // HW6:
    ASTREA_HW6 a6(.m_weight(_m_weight[2]), .frame_flips(_frame_flips[2]), .edges(edges[14:0]));

    assign edge_count = $countones(edges_v);

    always_comb
    begin
        m_weight = '0;
        frame_flips = '0;

        if (edge_count == 1) // hw 2
        begin
            m_weight = _m_weight[0];
            frame_flips = _frame_flips[0];
        end
        else if (edge_count == 6) // hw 4
        begin
            m_weight = _m_weight[1];
            frame_flips = _frame_flips[1];
        end
        else if (edge_count == 15) // hw 6
        begin
            m_weight = _m_weight[2];
            frame_flips = _frame_flips[2];
        end
        else
        begin
            $fatal(1, "ASTREA: unexpected edge count: %0d\n", edge_count);
        end
    end
endmodule

////////////////////////////////////////////
////////////////////////////////////////////

// Hamming weight 4 decoder
module ASTREA_HW4
(
    output weight_type   m_weight,
    output obs_type      frame_flips,
    input  edge_type[5:0] edges
);
    edge_type[2:0] m;
    edge_type cmp_win_round_one;
    edge_type winner;

    assign m_weight = winner.w;
    assign frame_flips = winner.f;

    // comparison logic for `m` into `cmp_win_round_one`
    assign cmp_win_round_one = (m[0].w < m[1].w) ? m[0] : m[1];

    // then `winner` is the better one of `cmp_win_round_one`
    // and `m[2]` 
    assign winner = (cmp_win_round_one.w < m[2].w) ? cmp_win_round_one : m[2];

    // Combinatorial logic for matching combinations:
    // m[0] = 01 + 23 = idx 0 + idx 5
    // m[1] = 02 + 13 = idx 1 + idx 4
    // m[2] = 03 + 12 = idx 2 + idx 3
    _edge_acc2 a1(m[0], edges[0], edges[5]);
    _edge_acc2 a2(m[1], edges[1], edges[4]);
    _edge_acc2 a3(m[2], edges[2], edges[3]);
endmodule

////////////////////////////////////////////
////////////////////////////////////////////

// Hamming weight 6 decoder
module ASTREA_HW6
(
    output weight_type    m_weight,
    output obs_type       frame_flips,
    input  edge_type[14:0] edges
);
    edge_type[14:0] m;
    edge_type[7:0]  cmp_win_l1;  // result of pairwise comparison -- 7 winners (+1 which is m[14])
    edge_type[3:0]  cmp_win_l2;  // four winners in L2
    edge_type[1:0]  cmp_win_l3;  // two winners in L3
    edge_type winner;

    // final winner's frame flips are output:
    assign m_weight = winner.w;
    assign frame_flips = winner.f;

    assign winner = (cmp_win_l3[0].w < cmp_win_l3[1].w) ? cmp_win_l3[0] : cmp_win_l3[1];

    // assign last winner of L1 comparisons to be `m[14]` (odd one out)
    assign cmp_win_l1[7] = m[14];

    // L1, L2, and L3 winners:
    MASS_EDGE_COMPARE #( .N(7) ) cmp_l1(cmp_win_l1[6:0], m[13:0]);
    MASS_EDGE_COMPARE #( .N(4) ) cmp_l2(cmp_win_l2, cmp_win_l1);
    MASS_EDGE_COMPARE #( .N(2) ) cmp_l3(cmp_win_l3, cmp_win_l2);

    // Combinatorial logic for matching combinations:
    // idx: 01=0  02=1  03=2  04=3  05=4
    //      12=5  13=6  14=7  15=8
    //      23=9  24=10 25=11
    //      34=12 35=13
    //      45=14
    _edge_acc3 a1(m[0],  edges[0],  edges[9],  edges[14]); // 01+23+45
    _edge_acc3 a2(m[1],  edges[0],  edges[10], edges[13]); // 01+24+35
    _edge_acc3 a3(m[2],  edges[0],  edges[11], edges[12]); // 01+25+34
    _edge_acc3 a4(m[3],  edges[1],  edges[6],  edges[14]); // 02+13+45
    _edge_acc3 a5(m[4],  edges[1],  edges[7],  edges[13]); // 02+14+35
    _edge_acc3 a6(m[5],  edges[1],  edges[8],  edges[12]); // 02+15+34
    _edge_acc3 a7(m[6],  edges[2],  edges[5],  edges[14]); // 03+12+45
    _edge_acc3 a8(m[7],  edges[2],  edges[7],  edges[11]); // 03+14+25
    _edge_acc3 a9(m[8],  edges[2],  edges[8],  edges[10]); // 03+15+24
    _edge_acc3 a10(m[9],  edges[3],  edges[5],  edges[13]); // 04+12+35
    _edge_acc3 a11(m[10], edges[3],  edges[6],  edges[11]); // 04+13+25
    _edge_acc3 a12(m[11], edges[3],  edges[8],  edges[9]);  // 04+15+23
    _edge_acc3 a13(m[12], edges[4],  edges[5],  edges[12]); // 05+12+34
    _edge_acc3 a14(m[13], edges[4],  edges[6],  edges[10]); // 05+13+24
    _edge_acc3 a15(m[14], edges[4],  edges[7],  edges[9]);  // 05+14+23
endmodule

////////////////////////////////////////////
////////////////////////////////////////////

// Module that sums N edges together.
// The weight (`w`) is added as an integer. The
// frame (`f`) is XOR'd.
module EDGE_ACC #(parameter int N = 2)
(
    output edge_type          sum,
    input  edge_type[N-1:0]   edges
);
    always_comb begin
        sum.w = '0;
        sum.f = '0;
        for (int i = 0; i < N; i++)
        begin
            sum.w = sum.w + edges[i].w;
            sum.f = sum.f ^ edges[i].f;
        end
    end
endmodule

module _edge_acc2(output edge_type s, input edge_type a, input edge_type b);
    EDGE_ACC #(.N(2)) adder( .sum(s), .edges('{a,b}) );
endmodule

module _edge_acc3(output edge_type s, input edge_type a, input edge_type b, input edge_type c);
    EDGE_ACC #(.N(3)) adder( .sum(s), .edges('{a,b,c}) );
endmodule

////////////////////////////////////////////
////////////////////////////////////////////

// Module used by Astrea that compares pairs of
// edges from a list of `2*N` edges and puts the
// smaller weighted edge to `out`
module MASS_EDGE_COMPARE #(parameter int N)
(
    output edge_type[N-1:0]   out,
    input  edge_type[2*N-1:0] cmp
);
    always_comb
        for (int i = 0; i < N; i++)
            out[i] = (cmp[2*i].w < cmp[2*i+1].w) ? cmp[2*i] : cmp[2*i+1];
endmodule

////////////////////////////////////////////
////////////////////////////////////////////
