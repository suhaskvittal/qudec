/*
 *  contributor: OpenAI Codex
 * */

#include "decoder/libra.h"
#include <pymatching/sparse_blossom/driver/user_graph.h>
#include <pymatching/sparse_blossom/driver/mwpm_decoding.h>
#include <stim.h>
#include <iostream>
#include <cmath>
#include <array>
#include <limits>
#include <stdexcept>
#include <sstream>
int main(){
 // Two CSS graphs with correlated mechanisms and both logical boundaries.
 stim::DetectorErrorModel dem(R"(
 error(0.025) D0 L0
 error(0.04) D1
 error(0.015) D0 D1
 error(0.035) D2
 error(0.022) D3
 error(0.013) D2 D3
 error(0.009) D0 L0 ^ D2
 error(0.011) D1 ^ D3
 error(0.004) D0 D1 ^ D2 D3
 )");
 struct E{unsigned dets;bool obs;double w;};std::vector<E> es;
 dem.iter_flatten_error_instructions([&](const auto& x){E e{};e.w=std::log((1-x.arg_data[0])/x.arg_data[0]);for(auto t:x.target_data){if(t.is_relative_detector_id())e.dets^=1U<<t.val();else if(t.is_observable_id())e.obs=!e.obs;}es.push_back(e);});
 std::array<std::array<double,2>,16> best,mass{};for(auto& a:best)a.fill(std::numeric_limits<double>::infinity());
 for(unsigned m=0;m<(1U<<es.size());++m){unsigned s=0;bool o=false;double w=0;for(size_t i=0;i<es.size();++i)if(m>>i&1){s^=es[i].dets;o^=es[i].obs;w+=es[i].w;}best[s][o]=std::min(best[s][o],w);mass[s][o]+=std::exp(-w);}
 dec::Libra::Config cfg;cfg.validate=true;cfg.ensemble_size=60;cfg.gap_threshold_db=1000;cfg.degeneracy=false;cfg.latency_histogram=true;dec::Libra libra(dem,cfg);
 for(unsigned s=1;s<16;++s){dec::SyndromeType syn(4);for(unsigned d=0;d<4;++d)syn[d]=(s>>d)&1;auto r=libra.decode_detailed(syn);for(int o=0;o<2;++o)if(std::abs(r.weights[o]-best[s][o])>1e-6)throw std::runtime_error("small hypergraph minimum-weight oracle failed");}
 std::cout<<"exhaustive 512 configurations, 15 syndromes, both classes passed\n";
 cfg.degeneracy=true;dec::Libra deg(dem,cfg);
 for(unsigned s=1;s<16;++s){dec::SyndromeType syn(4);for(unsigned d=0;d<4;++d)syn[d]=(s>>d)&1;auto r=deg.decode_detailed(syn);for(int o=0;o<2;++o){if(r.log_masses[o]< -r.weights[o]-1e-8||r.log_masses[o]>std::log(mass[s][o])+1e-8)throw std::runtime_error("degeneracy probability oracle failed");}}
 std::cout<<"degeneracy masses are valid subsets of exact class masses\n";
 // The actual observable must never be used to decide a correction.
 dec::SyndromeType syn(4),o0(1),o1(1);syn[0]=1;o1[0]=1;auto a=libra.decode(syn,o0),b=libra.decode(syn,o1);if(a.flipped_obs[0]!=b.flipped_obs[0])throw std::runtime_error("observable leakage");
 std::ostringstream report;libra.print_stats(report);
 if(report.str().find("LIBRA_DECODE_LATENCY_US =================================")==std::string::npos ||
    report.str().find("1000 logarithmic bins")==std::string::npos ||
    report.str().find("count = 2")==std::string::npos)
  throw std::runtime_error("Libra latency histogram is missing timed shots");
 cfg.latency_histogram=false;dec::Libra unprofiled(dem,cfg);
 unprofiled.decode(syn,o0);
 std::ostringstream unprofiled_report;unprofiled.print_stats(unprofiled_report);
 if(unprofiled_report.str().find("LIBRA_DECODE_LATENCY_US")!=std::string::npos)
  throw std::runtime_error("Disabled Libra latency histogram was printed");
 std::cout<<"truth-independence passed\n";
 // A zero gap threshold always takes the paper's unconstrained correlated
 // fallback. Its prediction must equal that matcher's boundary logical parity.
 cfg.gap_threshold_db=0;dec::Libra fallback(dem,cfg);
 auto graph=pm::detector_error_model_to_mwpm(dem,pm::NUM_DISTINCT_WEIGHTS,true,true);
 for(unsigned s=1;s<16;++s){
  dec::SyndromeType ss(4);std::vector<uint64_t> events;
  for(unsigned d=0;d<4;++d)if(s>>d&1){ss[d]=1;events.push_back(d);}
  std::vector<int64_t> edges;pm::decode_detection_events_to_edges_with_edge_correlations(graph,events,edges);
  bool prediction=false;for(size_t k=0;k<edges.size();k+=2)prediction^=(edges[k]==0&&edges[k+1]==-1)||(edges[k]==-1&&edges[k+1]==0);
  auto r=fallback.decode_detailed(ss);if(!r.fallback_used||r.ensemble_used||r.prediction!=prediction)throw std::runtime_error("unconstrained correlated fallback oracle failed");
 }
 std::cout<<"paper's large-gap correlated fallback passed\n";
}
