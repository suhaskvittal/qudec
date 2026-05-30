//  author: Suhas Vittal
//  date:   28 May 2026

`include "cluster_match/_generated.svh"

////////////////////////////////////////////
////////////////////////////////////////////

typedef logic[`D_BIT_WIDTH-1:0]  detector_id_type;
typedef logic[`W_BIT_WIDTH:0]    weight_type;
typedef logic[`OBS_BIT_WIDTH:0]  obs_type;

////////////////////////////////////////////
////////////////////////////////////////////

// `d_rbf_entry` and `d_fbf_entry` correspond to
// entries in frontend structures that buffer active
// detectors for a single round of syndromes and 
// all round of syndromes, respectively.
//
// R = round, F = full
//
// `d_rbf_type` and `d_fbf_type` are the buffers
// contains these entries.
typedef logic[`D_COUNT_PER_ROUND_BIT_WIDTH-1:0] d_rbf_entry;
typedef logic[`D_COUNT_BIT_WIDTH-1:0]           d_fbf_entry;

typedef struct packed
{
    d_rbf_entry             d[`RBF_SIZE-1:0];
    logic[`LG_RBF_SIZE-1:0] ptr;
} d_rbf_type;

typedef struct packed
{
    d_fbf_entry             d[`FBF_SIZE-1:0];
    logic[`LG_FBF_SIZE-1:0] ptr;
} d_fbf_type;

////////////////////////////////////////////
////////////////////////////////////////////

typedef struct packed
{
    detector_id_type d;
    weight_type      w;
    obs_type         f;
} neighbor_data_type;

typedef struct packed
{
    logic[`ADJ_MAX_DEGRE-1:0]  v;  // valid bits
    neighbor_data_type         n[`ADJ_MAX_DEGREE-1:0];
} adj_list_type;

////////////////////////////////////////////
////////////////////////////////////////////

// Data structure that implements the data for an edge
// passed into `ASTREA`
typedef struct packed
{
    weight_type w;
    obs_type f;
} edge_type;

////////////////////////////////////////////
////////////////////////////////////////////
