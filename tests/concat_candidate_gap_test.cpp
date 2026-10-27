// Author: OpenAI Codex. Candidate-frame gap correctness and truth-independent selection.
#include "circuit_generator.h"
#include "decoder/concat_mwpm.h"
#include "decoder/concat_mwpm/circuit.h"
#include <stim/util_top/circuit_to_dem.h>
#include <stim/simulators/dem_sampler.h>
#include <cmath>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>

void require(bool ok,const char* message)
{ if (!ok) throw std::runtime_error(message); }
int main() try
{
    const auto dem=stim::circuit_to_dem(cc666_code_capacity(5,.08,false),{});
    dec::concat::CircuitDecoder reference(dem,false,true,16,false,{10,20261024},{},false,20,true);
    dec::ConcatMWPM decoder(dem,{.ensemble_size=10,.belief_propagation=true});
    stim::DemSampler<64> sampler(dem,std::mt19937_64(78263),256);
    sampler.resample(false);
    const auto shots=sampler.det_buffer.transposed();
    size_t available=0,missing=0;
    for (size_t shot=0;shot<256;shot++)
    {
        std::vector<uint8_t> bits(dem.count_detectors());
        for (size_t d=0;d<bits.size();d++) bits[d]=shots[shot][d];
        const auto inspection=reference.inspect(bits);
        const auto& best=inspection.correction;
        double expected=INFINITY;
        for (const auto& c:inspection.candidates) if (c.observable!=best.observable)
            expected=std::min(expected,std::max(0.0,c.weight-best.weight));
        dec::ObsType correct(1),wrong(1);
        correct[0]=best.observable;wrong[0]=!best.observable;
        const auto a=decoder.decode(shots[shot],correct),b=decoder.decode(shots[shot],wrong);
        require(a.flipped_obs[0]==best.observable && b.flipped_obs[0]==best.observable,"Truth changed prediction");
        require(a.has_candidate_gap==std::isfinite(expected),"Missing frame treated as finite gap");
        require(a.has_candidate_gap==b.has_candidate_gap,"Truth changed gap availability");
        if (a.has_candidate_gap)
        {
            available++;
            const double db=expected*10/std::log(10.0);
            require(std::abs(a.candidate_gap_db-db)<1e-9,"Incorrect candidate gap or dB conversion");
            require(std::abs(a.matching_data.gap-db)<1e-9 && std::abs(b.matching_data.gap+db)<1e-9,"Incorrect signed gap");
        }
        else missing++;
    }
    require(available && missing,"Test must exercise both gap availability states");
    require(decoder.statistics().candidate_gap_available==2*available &&
            decoder.statistics().candidate_gap_unavailable==2*missing,"Incorrect gap counters");
    require(std::abs(decoder.s_signed_gap.mean())<1e-10,"Signed gap histogram did not balance correct and wrong shots");
    require(std::abs(decoder.s_unsigned_gap.mean()-decoder.s_unsigned_gap_errors.mean())<1e-10,"Error histogram differs from paired gap magnitudes");
    // Single-frame vacuum and history: a missing alternative must clear the previous gap.
    dec::SyndromeType zero(dem.count_detectors());dec::ObsType truth(1);
    require(!decoder.decode(zero,truth).has_candidate_gap,"Vacuum should not imply infinite confidence");
    std::ostringstream report;decoder.print_stats(report);
    require(report.str().find("CONCAT_BP_SIGNED_GAP_DB")!=std::string::npos,"Gap histogram missing from output");
    // The gated path must use only candidates that actually ran.
    dec::ConcatMWPM gated(dem,{.ensemble_size=10,.gap_gating=true,.gap_gating_threshold_db=0,.belief_propagation=true});
    dec::ConcatMWPM single(dem,{.belief_propagation=true});
    for (size_t i=0;i<32;i++)
    {
        const auto a=gated.decode(shots[i],truth),b=single.decode(shots[i],truth);
        require(a.flipped_obs[0]==b.flipped_obs[0] && a.has_candidate_gap==b.has_candidate_gap,"Gating used unseen candidates");
        if (a.has_candidate_gap) require(std::abs(a.candidate_gap_db-b.candidate_gap_db)<1e-9,"Gated candidate gap changed");
    }
    // Real flagged circuit: posterior weights, debug path, and both matching stages.
    const auto circuit_dem=stim::circuit_to_dem(cc666_si1000(5,5,.001,false),{});
    dec::concat::CircuitDecoder bp_reference(circuit_dem,false,true,16,false,{3,20261024},{},false,20,true);
    dec::ConcatMWPM bp_decoder(circuit_dem,{.ensemble_size=3,.debug=true,.belief_propagation=true});
    stim::DemSampler<64> circuit_sampler(circuit_dem,std::mt19937_64(92815),32);
    circuit_sampler.resample(false);
    const auto circuit_shots=circuit_sampler.det_buffer.transposed();
    for (size_t i=0;i<32;i++)
    {
        std::vector<uint8_t> bits(circuit_dem.count_detectors());
        for (size_t d=0;d<bits.size();d++) bits[d]=circuit_shots[i][d];
        const auto inspected=bp_reference.inspect(bits);
        double expected=INFINITY;
        for (const auto& c:inspected.candidates) if (c.observable!=inspected.correction.observable)
            expected=std::min(expected,std::max(0.0,c.weight-inspected.correction.weight));
        const auto result=bp_decoder.decode(circuit_shots[i],truth);
        require(result.flipped_obs[0]==inspected.correction.observable,"BP debug altered prediction");
        require(result.has_candidate_gap==std::isfinite(expected),"Flag-conditioned alternative missing");
        if (result.has_candidate_gap)
            require(std::abs(result.candidate_gap_db-expected*10/std::log(10.0))<1e-9,"BP posterior gap incorrect");
        require(result.debug_info.find("Candidate frame gap")!=std::string::npos,"Debug confidence missing");
    }
    std::cout<<"Candidate gap tests passed: available="<<available<<" missing="<<missing<<'\n';
}
catch (const std::exception& e) { std::cerr<<e.what()<<'\n';return 1; }
