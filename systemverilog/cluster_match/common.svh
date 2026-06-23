//  author: Suhas Vittal
//  date:   28 May 2026

`include "cluster_match/_generated.svh"

////////////////////////////////////////////
////////////////////////////////////////////

typedef logic[`D_BIT_WIDTH-1:0]    detector_type;
typedef logic[`W_BIT_WIDTH-1:0]    weight_type;
typedef logic[`OBS_BIT_WIDTH-1:0]  obs_type;

typedef detector_type[`D_MAX-1:0] syndrome_type;

////////////////////////////////////////////
////////////////////////////////////////////

typedef struct packed
{
    detector_type d;
    weight_type   w;
    obs_type      f;
} neighbor_data_type;

typedef neighbor_data_type[`ADJ_MAX_DEGREE-1:0]  adj_type;

////////////////////////////////////////////
////////////////////////////////////////////

// Data structure that implements the data for an edge
// passed into `ASTREA`
typedef struct packed
{
    obs_type f;
    weight_type w;
} edge_type;

////////////////////////////////////////////
////////////////////////////////////////////
