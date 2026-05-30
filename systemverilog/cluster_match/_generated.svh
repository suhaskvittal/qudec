`define D_COUNT                      7800
`define ADJ_MAX_DEGREE               12

`define D_COUNT_PER_ROUND            312
`define D_COUNT_PER_ROUND_BIT_WIDTH  9

`define D_BIT_WIDTH                  13
`define W_BIT_WIDTH                  8
`define OBS_BIT_WIDTH                1

`define EXTRACT_ACTIVE_CHUNK_WIDTH   32
`define EXTRACT_ACTIVE_OUT_COUNT     10   // ceiling of `D_COUNT_PER_ROUND/EXTRACT_ACTIVE_CHUNK_WIDTH`

`define RBF_SIZE    16
`define LG_RBF_SIZE 4

`define FBF_SIZE    256
`define LG_FBF_SIZE 8
