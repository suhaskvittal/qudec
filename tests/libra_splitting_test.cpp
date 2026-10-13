/*
 *  contributor: OpenAI Codex
 * */

#include "libra_splitting.h"
#include <iostream>
#include <array>
int main(){stim::DetectorErrorModel dem("error(0.02) D0 L0\nerror(0.11) D1\nerror(0.23) D2\n");LibraSplittingModel model(dem);
 auto event=[](dec::SyndromeRef syn,dec::ObsRef){return LibraSplittingModel::Labels{bool(syn[0]),bool(syn[0])};};
 std::mt19937_64 rng(441);auto s=model.make_state({0},event);double theta=0.5,to=0.;double ratio_sum=0;std::array<double,3> incidence{};
 for(int i=0;i<300000;++i){model.step(s,theta,rng,event);if(i<10000)continue;ratio_sum+=std::exp(model.log_ratio(theta,to,s.active.size()));for(size_t j=0;j<3;++j)incidence[j]+=s.position[j]>=0;}
 for(size_t j=0;j<3;++j){double p=j==0?1.:std::exp(theta)*model.odds[j]/(1+std::exp(theta)*model.odds[j]);if(std::abs(incidence[j]/290000-p)>0.005)throw std::runtime_error("conditional MCMC distribution check failed");}
 double ptheta=std::exp(theta)*model.odds[0]/(1+std::exp(theta)*model.odds[0]);double exact_ratio=model.probabilities[0]/ptheta;
 if(std::abs(ratio_sum/290000-exact_ratio)>0.005)throw std::runtime_error("splitting ratio oracle failed");std::cout<<"conditional distribution and likelihood-ratio oracle passed\n";}
