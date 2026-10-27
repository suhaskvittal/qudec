// Author: OpenAI Codex. Two-phase matching and grouped-sampler regressions.
#include "circuit_generator.h"
#include "decoder/concat_mwpm.h"
#include "decoder/concat_mwpm/circuit.h"
#include "sampler.h"
#include <stim/util_top/circuit_to_dem.h>
#include <stim/simulators/dem_sampler.h>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <tuple>
namespace {
void require(bool ok,const char* message) {if (!ok) throw std::runtime_error(message);}
void partitions() {
    dec::concat::Model model{{0,1,2,0,1,2},
        {{.01,{0,1,2},0},{.02,{3,4,5},1},{.01,{0,3},0},{.01,{1,4},0},{.01,{2,5},0},
         {.01,{0,1},1},{.01,{3,4},0},{.01,{0},1},{.01,{1},0},{.01,{2},0},
         {.01,{3},0},{.01,{4},0},{.01,{5},1}}};
    for (bool boundary:{false,true}) for (size_t ensemble:{1u,5u}) {
        dec::concat::Decoder serial(model,false,boundary,{ensemble,134});
        for (size_t workers:{1u,2u,3u,8u}) {
            std::vector<std::unique_ptr<dec::concat::Decoder>> parallel;
            for (size_t r=0;r<workers;r++) parallel.push_back(std::make_unique<dec::concat::Decoder>(model,false,boundary,dec::concat::EnsembleConfig{ensemble,134}));
            for (size_t mask=0;mask<64;mask++) {
                std::vector<uint8_t> syndrome(6);
                for (size_t v=0;v<6;v++) syndrome[v]=(mask>>v)&1;
                const auto expected=serial.decode_candidates(syndrome,true);
                std::vector<dec::concat::FirstStageResult> first;
                for (size_t r=0;r<workers;r++) {
                    auto local=parallel[r]->decode_first_stage(syndrome,r,workers,true);
                    first.insert(first.end(),local.begin(),local.end());
                }
                auto key=[](const auto& v){return std::tuple{v.member,v.color,v.boundary_mask};};
                std::sort(first.begin(),first.end(),[&](const auto& a,const auto& b){return key(a)<key(b);});
                std::vector<dec::concat::Result> candidates;
                size_t second_calls=0;
                for (size_t r=0;r<workers;r++) {
                    auto local=parallel[r]->decode_second_stage(syndrome,first,r,workers,true);
                    second_calls+=parallel[r]->statistics.second_calls;
                    candidates.insert(candidates.end(),local.begin(),local.end());
                }
                std::sort(candidates.begin(),candidates.end(),[&](const auto& a,const auto& b){return key(a)<key(b);});
                require(candidates.size()==expected.size(),"Distributed coalescing changed candidate count");
                require(second_calls==serial.statistics.second_calls,"Distributed decoding repeated second-stage solves");
                for (size_t j=0;j<candidates.size();j++) {
                    const auto &a=candidates[j],&b=expected[j];
                    require(key(a)==key(b) && a.observable==b.observable && a.weight==b.weight &&
                        a.events==b.events && a.restricted_edges==b.restricted_edges && a.second_stage_paths==b.second_stage_paths,
                        "Parallel partition differs from serial candidate");
                }
            }
        }
    }
}
#if defined(ENABLE_MPI)
void grouped(int rank,int size) {
    if (size!=6) return;
    const auto dem=stim::circuit_to_dem(cc666_code_capacity(3,.03,false),{});
    dec::concat::CircuitDecoder diagnostic(dem,false);
    for (size_t cores:{2u,3u}) for (bool boundary:{false,true}) {
        dec::ConcatMWPM parallel(dem,{.boundary_gap=boundary,.ensemble_size=5,.debug=true,.decoding_cores=cores});
        dec::ConcatMWPM serial(dem,{.boundary_gap=boundary,.ensemble_size=5});
        const bool root=rank%cores==0;
        for (size_t shot=0;shot<40;shot++) {
            const auto bits=diagnostic.fault_syndrome({(shot+rank/cores)%diagnostic.faults.size()});
            dec::SyndromeType syndrome(bits.size());dec::ObsType truth(1);
            if (root) for (size_t v=0;v<bits.size();v++) syndrome[v]=bits[v];
            auto actual=parallel.decode(syndrome,truth);
            if (root) {
                require(actual.flipped_obs==serial.decode(syndrome,truth).flipped_obs,"MPI frame differs from serial");
                require(actual.debug_info.find("Stage 1 matching")!=std::string::npos &&
                    actual.debug_info.find("Stage 2 matching")!=std::string::npos &&
                    actual.debug_info.find("Winner:")!=std::string::npos,"MPI debug paths missing");
            } else require(actual.debug_info.empty() && parallel.statistics().shots==0,"Helper returned/counts a sampled result");
        }
    }
    // Real circuit faults exercise flag-selected enabled edges, probabilities,
    // virtual IDs and logical frames across repeated persistent graph updates.
    const auto circuit_dem=stim::circuit_to_dem(cc666_si1000(5,5,.001,false),{});
    for (size_t cores:{2u,3u}) for (bool boundary:{false,true}) for (bool renormalize:{false,true}) {
        dec::ConcatMWPM parallel(circuit_dem,{.boundary_gap=boundary,.ensemble_size=5,
            .flag_renormalization=renormalize,.physical_error_rate=.001,.decoding_cores=cores});
        dec::ConcatMWPM serial(circuit_dem,{.boundary_gap=boundary,.ensemble_size=5,
            .flag_renormalization=renormalize,.physical_error_rate=.001});
        SyndromeTable shots(0,0),truths(0,0);
        const bool root=rank%cores==0;
        if (root) {
            stim::DemSampler<64> sampler(circuit_dem,std::mt19937_64(481+rank/cores),64);
            sampler.resample(false);shots=sampler.det_buffer.transposed();truths=sampler.obs_buffer.transposed();
        }
        dec::SyndromeType helper(circuit_dem.count_detectors());dec::ObsType helper_truth(1);
        for (size_t shot=0;shot<64;shot++) {
            auto syndrome=root ? shots[shot]:dec::SyndromeRef(helper);
            auto truth=root ? truths[shot]:dec::ObsRef(helper_truth);
            const auto result=parallel.decode(syndrome,truth);
            if (root) require(result.flipped_obs==serial.decode(syndrome,truth).flipped_obs,
                "Flag-conditioned MPI prediction differs from serial");
        }
    }
    for (size_t cores:{2u,3u}) for (int64_t shots:{1,5,2051}) {
        dec::ConcatMWPM decoder(dem,{.ensemble_size=3,.decoding_cores=cores});
        ExperimentConfig conf;conf.monte_carlo.max_samples=shots;conf.monte_carlo.stop_at_error_count=10000;conf.verbosity=1;
        estimate_logical_error_rate(dem,decoder,conf,[&](auto,auto,auto){require(rank%cores==0,"Helper invoked MC callback");});
        uint64_t count=decoder.statistics().shots;
        MPI_Allreduce(MPI_IN_PLACE,&count,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
        require(count==uint64_t(shots),"Grouped MC global count incorrect");
    }
    for (size_t cores:{2u,3u}) {
        dec::ConcatMWPM decoder(dem,{.decoding_cores=cores});
        ExperimentConfig conf;conf.method="rare_event";conf.rare_event.max_level=dem.count_errors()+1;
        conf.rare_event.samples_per_level=5;conf.rare_event.max_errors_per_level=10000;
        conf.rare_event.relative_tail_tolerance=0;conf.verbosity=1;
        estimate_logical_error_rate(dem,decoder,conf,[&](auto,auto,auto){require(rank%cores==0,"Helper invoked rare callback");});
        uint64_t count=decoder.statistics().shots;
        MPI_Allreduce(MPI_IN_PLACE,&count,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
        require(count==5*(dem.count_errors()+1),"Grouped rare-event global count incorrect");
    }
}
#endif
}
int main(int argc,char** argv) {
    int rank=0,size=1;
#if defined(ENABLE_MPI)
    MPI_Init(&argc,&argv);MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
#endif
    try {
        if (size==1) partitions();
#if defined(ENABLE_MPI)
        grouped(rank,size);
#endif
        if (!rank) std::cout<<"Two-phase matching, pre-second-stage coalescing and grouped samplers passed\n";
    } catch (const std::exception& e) {
        std::cerr<<"rank "<<rank<<": "<<e.what()<<'\n';
#if defined(ENABLE_MPI)
        MPI_Abort(MPI_COMM_WORLD,1);
#endif
        return 1;
    }
#if defined(ENABLE_MPI)
    MPI_Finalize();
#endif
}
