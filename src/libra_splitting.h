#ifndef QUDEC_LIBRA_SPLITTING_H
#define QUDEC_LIBRA_SPLITTING_H

#include "decoder/types.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>

// Exact Metropolis kernel for the independent DEM Bernoulli distribution,
// restricted to the event that at least one of the two decoders fails.
// Intermediate distributions multiply every original odds ratio by exp(theta).
// Their likelihood ratios depend only on the number of active instructions.
class LibraSplittingModel {
public:
    struct Labels {bool pymatching=false,libra=false;bool any()const{return pymatching||libra;}};
    struct Pattern {std::vector<uint32_t> detectors;bool observable=false;};
    struct State {
        std::vector<size_t> active;
        std::vector<int32_t> position;
        dec::SyndromeType syndrome;
        dec::ObsType truth;
        Labels labels;
        double active_odds=0;
        uint64_t attempts=0,physical_accepts=0,event_accepts=0;
        State(size_t errors,size_t detectors):position(errors,-1),syndrome(detectors),truth(1){}
    };
    const size_t num_detectors;
    std::vector<double> probabilities,odds,cumulative;
    std::vector<Pattern> patterns;
    double total_odds=0;

    explicit LibraSplittingModel(const stim::DetectorErrorModel& dem):num_detectors(dem.count_detectors()) {
        if(dem.count_observables()!=1)throw std::invalid_argument("Splitting requires one observable");
        dem.iter_flatten_error_instructions([&](const auto& inst){
            double p=inst.arg_data[0];if(!(p>=0&&p<1))throw std::invalid_argument("Invalid Bernoulli probability");
            probabilities.push_back(p);odds.push_back(p/(1-p));total_odds+=odds.back();cumulative.push_back(total_odds);
            Pattern pat;for(auto t:inst.target_data){if(t.is_relative_detector_id())pat.detectors.push_back(t.val());else if(t.is_observable_id())pat.observable=!pat.observable;}
            std::sort(pat.detectors.begin(),pat.detectors.end());std::vector<uint32_t> ds;
            for(size_t i=0;i<pat.detectors.size();){size_t j=i+1;while(j<pat.detectors.size()&&pat.detectors[i]==pat.detectors[j])++j;if((j-i)&1)ds.push_back(pat.detectors[i]);i=j;}
            pat.detectors=std::move(ds);patterns.push_back(std::move(pat));
        });
        if(total_odds<=0)throw std::invalid_argument("Splitting needs nonzero error channels");
    }
    double log_normaliser(double theta)const {
        double out=0,t=std::expm1(theta);for(double p:probabilities)out+=std::log1p(p*t);return out;
    }
    double log_ratio(double from,double to,size_t count)const {
        return (to-from)*count+log_normaliser(from)-log_normaliser(to);
    }
    template<class Decide> State make_state(const std::vector<size_t>& ids,const Decide& decide)const {
        State s(odds.size(),num_detectors);
        for(size_t i:ids){if(i>=odds.size()||s.position[i]>=0)throw std::invalid_argument("Invalid initial error set");s.position[i]=s.active.size();s.active.push_back(i);s.active_odds+=odds[i];flip(s,i);}
        s.labels=decide(s.syndrome,s.truth);if(!s.labels.any())throw std::invalid_argument("Initial splitting state is not a decoder failure");return s;
    }
    template<class Decide> bool step(State& s,double theta,std::mt19937_64& rng,const Decide& decide)const {
        ++s.attempts;std::uniform_real_distribution<double> unit(0,1);
        bool add=(rng()&1);size_t id;double acceptance;
        if(add){
            if(s.active.size()==odds.size())return false;
            double inactive=total_odds-s.active_odds;if(inactive<=0)return false;
            // Rejection from the static odds distribution gives the exact
            // distribution odds_i / sum_inactive_odds over inactive channels.
            do {double x=unit(rng)*total_odds;id=std::upper_bound(cumulative.begin(),cumulative.end(),x)-cumulative.begin();if(id==odds.size())id=odds.size()-1;}while(s.position[id]>=0||odds[id]==0);
            acceptance=std::exp(theta)*inactive/(s.active.size()+1);
        }else{
            if(s.active.empty())return false;
            std::uniform_int_distribution<size_t> which(0,s.active.size()-1);id=s.active[which(rng)];
            acceptance=s.active.size()/(std::exp(theta)*(total_odds-s.active_odds+odds[id]));
        }
        if(unit(rng)>=std::min(1.,acceptance))return false;
        ++s.physical_accepts;flip(s,id);auto labels=decide(s.syndrome,s.truth);
        if(!labels.any()){flip(s,id);return false;}
        ++s.event_accepts;s.labels=labels;
        if(add){s.position[id]=s.active.size();s.active.push_back(id);s.active_odds+=odds[id];}
        else {size_t j=s.position[id],last=s.active.back();s.active[j]=last;s.position[last]=j;s.active.pop_back();s.position[id]=-1;s.active_odds-=odds[id];}
        return true;
    }
private:
    void flip(State& s,size_t id)const {for(uint32_t d:patterns[id].detectors)s.syndrome[d]^=1;s.truth[0]^=patterns[id].observable;}
};
#endif
