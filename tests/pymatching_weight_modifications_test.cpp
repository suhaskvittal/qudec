// Author: OpenAI Codex. Generic edge-weight API and paired correlated matching regression.
#include "circuit_generator.h"
#include <pymatching/sparse_blossom/driver/mwpm_decoding.h>
#include <pymatching/sparse_blossom/driver/user_graph.h>
#include <stim/util_top/circuit_to_dem.h>
#include <stim/simulators/dem_sampler.h>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
using namespace pm;
void require(bool ok,const char* message) {if(!ok)throw std::runtime_error(message);}
std::vector<weight_int> weights(const Mwpm& decoder) {
 std::vector<weight_int> out;
 for(const auto& n:decoder.flooder.graph.nodes)out.insert(out.end(),n.neighbor_weights.begin(),n.neighbor_weights.end());
 for(const auto& n:decoder.search_flooder.graph.nodes)out.insert(out.end(),n.neighbor_weights.begin(),n.neighbor_weights.end());
 return out;
}
// Literal two-pass flow from the frozen original mwpm_decoding.cc. Only second
// uncorrelated extraction calls the shared primitive; reweight/undo stay legacy.
MatchingResult legacy(Mwpm& decoder,const std::vector<uint64_t>& syndrome,
 std::vector<weight_int>* changed=nullptr) {
 std::vector<int64_t> edges;decode_detection_events_to_edges(decoder,syndrome,edges);
 decoder.flooder.graph.reweight_for_edges(edges);decoder.search_flooder.graph.reweight_for_edges(edges);
 if(changed)*changed=weights(decoder);
 auto result=decode_detection_events_for_up_to_64_observables(decoder,syndrome,false);
 decoder.flooder.graph.undo_reweights();decoder.search_flooder.graph.undo_reweights();return result;
}
std::vector<int64_t> legacy_edges(Mwpm& decoder,const std::vector<uint64_t>& syndrome) {
 std::vector<int64_t> edges;decode_detection_events_to_edges(decoder,syndrome,edges);
 decoder.flooder.graph.reweight_for_edges(edges);decoder.search_flooder.graph.reweight_for_edges(edges);
 edges.clear();decode_detection_events_to_edges(decoder,syndrome,edges);
 decoder.flooder.graph.undo_reweights();decoder.search_flooder.graph.undo_reweights();return edges;
}
std::vector<std::pair<int64_t,int64_t>> canonical(const std::vector<int64_t>& edges) {
 std::vector<std::pair<int64_t,int64_t>> out;
 for(size_t i=0;i<edges.size();i+=2)out.push_back(std::minmax(edges[i],edges[i+1]));
 std::sort(out.begin(),out.end());return out;
}
void api() {
 const stim::DetectorErrorModel dem("error(0.1) D0 D1\nerror(0.2) D0 L0\nerror(0.3) D1\n");
 auto decoder=detector_error_model_to_mwpm(dem,NUM_DISTINCT_WEIGHTS,true,true);
 const auto base=weights(decoder);EdgeWeightModifications outer(decoder);
 const auto id=outer.edge_id(0,1);require(id==outer.edge_id(1,0),"Directed edge identity");
 const auto original=outer.edge_weight(id);outer.modify_edge_weight(id,original+100);
 require(outer.modified_edge_count()==1,"Unexpected touched count");
 {
  EdgeWeightModifications inner(decoder);inner.modify_edge_weight(1,0,0);
  require(outer.edge_weight(id)==0,"Nested override missing");
  inner.clear_edge_weight_modifications();require(outer.edge_weight(id)==original+100,"Nested clear lost enclosing weights");
 }
 outer.modify_edge_weight(id,original+200);require(outer.modified_edge_count()==1,"Repeated override duplicated undo");
 outer.modify_edge_weight(0,-1,0);outer.clear_edge_weight_modifications();require(weights(decoder)==base,"Clear did not restore every adjacency copy");
 for(auto bad:{-2ll,1ll,2147483648ll}) {
  bool rejected=false;try{outer.modify_edge_weight(id,bad);}catch(const std::invalid_argument&){rejected=true;}
  require(rejected,"Invalid integer override accepted");
 }
 // Selectively change paths and verify both prediction extraction and shortest
 // path reconstruction use the same overridden undirected weights.
 outer.modify_edge_weight(0,1,0);
 std::vector<int64_t> chosen;decode_detection_events_to_edges(decoder,{0,1},chosen);
 require(canonical(chosen)==std::vector<std::pair<int64_t,int64_t>>{{0,1}},"Path reconstruction ignored override");
 outer.clear_edge_weight_modifications();require(weights(decoder)==base,"Clear after solve failed");
 // Callback failure must roll back the nested correlated layer and preserve outer.
 outer.modify_edge_weight(0,1,original+100);const auto before=weights(decoder);
 uint8_t result=0;total_weight_int w=0;bool threw=false;
 try{decode_detection_events(decoder,{0,1},&result,w,true,[](const Mwpm&){throw std::runtime_error("callback");});}
 catch(const std::runtime_error&){threw=true;}
 require(threw && weights(decoder)==before,"Failed correlated decode leaked modifications");
 outer.clear_edge_weight_modifications();require(weights(decoder)==base,"Outer restoration failed");
 // The public API does not change signed-edge preprocessing.
 const stim::DetectorErrorModel signed_dem("error(0.7) D0 L0\nerror(0.2) D1\nerror(0.1) D0 D1\n");
 auto signed_decoder=detector_error_model_to_mwpm(signed_dem,NUM_DISTINCT_WEIGHTS,true,false);
 EdgeWeightModifications signed_updates(signed_decoder);bool rejected=false;
 try{signed_updates.modify_edge_weight(0,1,0);}catch(const std::invalid_argument&){rejected=true;}
 require(rejected,"Negative graph accepted without signed bookkeeping support");
 auto signed_old=detector_error_model_to_mwpm(signed_dem,NUM_DISTINCT_WEIGHTS,true,false);
 for(const auto& syndrome:std::vector<std::vector<uint64_t>>{{},{0},{1},{0,1}})
  require(decode_detection_events_for_up_to_64_observables(signed_old,syndrome,false)==decode_detection_events_for_up_to_64_observables(signed_decoder,syndrome,false),"Signed ordinary compatibility changed");
 // Exercise path-based observable extraction when more than64 frames exist.
 const stim::DetectorErrorModel wide_dem("error(0.01) D0 L64 ^ D1 L0\nerror(0.1) D0\nerror(0.2) D1\nerror(0.03) D0 D1\n");
 auto wide_old=detector_error_model_to_mwpm(wide_dem,NUM_DISTINCT_WEIGHTS,true,true);
 auto wide_new=detector_error_model_to_mwpm(wide_dem,NUM_DISTINCT_WEIGHTS,true,true);
 for(size_t repeat=0;repeat<4;repeat++)for(const auto& syndrome:std::vector<std::vector<uint64_t>>{{},{0},{1},{0,1}}) {
  std::vector<int64_t> edges;decode_detection_events_to_edges(wide_old,syndrome,edges);
  wide_old.flooder.graph.reweight_for_edges(edges);wide_old.search_flooder.graph.reweight_for_edges(edges);
  std::vector<uint8_t> a(65),b(65);total_weight_int wa=0,wb=0;
  decode_detection_events(wide_old,syndrome,a.data(),wa,false);
  wide_old.flooder.graph.undo_reweights();wide_old.search_flooder.graph.undo_reweights();
  decode_detection_events(wide_new,syndrome,b.data(),wb,true);
  require(a==b && wa==wb && weights(wide_old)==weights(wide_new),"Wide-observable correlated path differs");
 }
 std::cout<<"API tests passed\n";
}
int main(int argc,char**argv) try {
 api();const size_t shots=argc>1?std::stoul(argv[1]):256;
 std::cout<<"d,p,shots,legacy_failures,new_failures,prediction_disagreements,weight_disagreements,correction_checks\n";
 for(size_t d:{3,5,7})for(double p:{.001,.005}) {
  const auto circuit=sc_si1000(d,d,p,false,true);
  const auto dem=stim::circuit_to_dem(circuit,{.decompose_errors=true});
  auto old=detector_error_model_to_mwpm(dem,NUM_DISTINCT_WEIGHTS,true,true);
  auto actual=detector_error_model_to_mwpm(dem,NUM_DISTINCT_WEIGHTS,true,true);
  const auto base=weights(actual);size_t f0=0,f1=0,dis=0,wd=0,checked=0;
  stim::DemSampler<64> sampler(dem,std::mt19937_64(20261027+d),shots);sampler.resample(false);
  const auto det=sampler.det_buffer.transposed(),obs=sampler.obs_buffer.transposed();
  for(size_t shot=0;shot<shots;shot++) {
   std::vector<uint64_t> syndrome;for(size_t q=0;q<dem.count_detectors();q++)if(det[shot][q])syndrome.push_back(q);
   std::vector<weight_int> expected_weights;
   const auto a=legacy(old,syndrome,shot%100==0?&expected_weights:nullptr);
   MatchingResult b;
   if(shot%100==0) {
    std::vector<uint8_t> prediction(dem.count_observables());total_weight_int cost=0;
    decode_detection_events(actual,syndrome,prediction.data(),cost,true,[&](const Mwpm& g){require(weights(g)==expected_weights,"Second-pass graph weights differ");});
    b={bit_vector_to_obs_mask(prediction),cost};
    std::vector<int64_t> new_edges;decode_detection_events_to_edges_with_edge_correlations(actual,syndrome,new_edges);
    require(canonical(legacy_edges(old,syndrome))==canonical(new_edges),"Corrections differ");checked++;
    require(weights(actual)==base && weights(old)==base,"Weights not restored after decoding");
   } else b=decode_detection_events_for_up_to_64_observables(actual,syndrome,true);
   const bool truth=obs[shot][0];f0+=(bool(a.obs_mask&1)!=truth);f1+=(bool(b.obs_mask&1)!=truth);
   dis+=a.obs_mask!=b.obs_mask;wd+=a.weight!=b.weight;
  }
  require(!dis && !wd,"Paired logical predictions or weights disagree");
  std::cout<<d<<','<<p<<','<<shots<<','<<f0<<','<<f1<<','<<dis<<','<<wd<<','<<checked<<std::endl;
 }
} catch(const std::exception& e){std::cerr<<"FAILED: "<<e.what()<<'\n';return 1;}
