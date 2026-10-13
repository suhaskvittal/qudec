#include "circuit_generator.h"
#include "decoder/libra_variant.h"
#include "decoder/surface_code.h"
#include <stim.h>
#include <stim/util_top/circuit_to_dem.h>
#include <algorithm>
#include <atomic>
#include <bit>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <random>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <vector>

namespace {
using Clock=std::chrono::steady_clock;
struct Counts {
    uint64_t shots=0, original=0, alternative=0, pm=0;
    uint64_t both=0, original_only=0, alternative_only=0;
    uint64_t original_ensemble=0, alternative_ensemble=0;
    uint64_t original_ns=0, alternative_ns=0, pm_ns=0;
    Counts& operator+=(const Counts& b) {
        shots+=b.shots; original+=b.original; alternative+=b.alternative; pm+=b.pm;
        both+=b.both; original_only+=b.original_only; alternative_only+=b.alternative_only;
        original_ensemble+=b.original_ensemble; alternative_ensemble+=b.alternative_ensemble;
        original_ns+=b.original_ns; alternative_ns+=b.alternative_ns; pm_ns+=b.pm_ns;
        return *this;
    }
};
struct Worker {
    Counts total;
    std::map<size_t,Counts> syndrome, faults;
    std::map<std::pair<size_t,size_t>,Counts> joint;
    std::vector<std::string> failures;
};
uint64_t elapsed(Clock::time_point t) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now()-t).count();
}
size_t hw(stim::simd_bits_range_ref<64> bits,size_t n) {
    size_t count=0;
    for(size_t k=0;k<n/64;++k)count+=std::popcount(bits.u64[k]);
    if(n%64)count+=std::popcount(bits.u64[n/64]&((uint64_t(1)<<(n%64))-1));
    return count;
}
std::string indices(stim::simd_bits_range_ref<64> bits,size_t n) {
    std::ostringstream s;bool first=true;
    for(size_t k=0;k<(n+63)/64;++k){uint64_t v=bits.u64[k];
        while(v){size_t i=64*k+std::countr_zero(v);v&=v-1;
            if(i<n){if(!first)s<<';';s<<i;first=false;}
        }
    }
    return s.str();
}
void csv_header(std::ostream& out,const std::string& key) {
    out<<key<<",shots,original_errors,alternative_errors,pymatching_errors,both_libra_errors,original_only_errors,alternative_only_errors,original_ensemble,alternative_ensemble,original_decode_ns,alternative_decode_ns,pymatching_decode_ns\n";
}
void csv_counts(std::ostream& out,const Counts& c) {
    out<<','<<c.shots<<','<<c.original<<','<<c.alternative<<','<<c.pm<<','<<c.both
       <<','<<c.original_only<<','<<c.alternative_only<<','<<c.original_ensemble
       <<','<<c.alternative_ensemble<<','<<c.original_ns<<','<<c.alternative_ns<<','<<c.pm_ns<<'\n';
}
std::ofstream create(const std::string& path) {
    if(std::filesystem::exists(path))throw std::runtime_error("Refusing to overwrite "+path);
    std::ofstream out(path);if(!out)throw std::runtime_error("Cannot create "+path);return out;
}
}

int main(int argc,char** argv) {
    try {
        if(argc!=9&&argc!=10)throw std::invalid_argument(
            "usage: libra_variant_evaluate fallback|correlation DIST WORKERS MAX_SHOTS SAMPLE_SEED PREFIX ERROR_TARGET KEEP_FALLBACK(0|1) [ENSEMBLE_SEED]");
        const std::string mode=argv[1],prefix=argv[6];
        const int distance=std::stoi(argv[2]),workers=std::stoi(argv[3]);
        const uint64_t max_shots=std::stoull(argv[4]),seed=std::stoull(argv[5]),target=std::stoull(argv[7]);
        const uint64_t ensemble_seed=argc==10?std::stoull(argv[9]):seed;
        const bool keep_fallback=std::stoi(argv[8])!=0;
        if((mode!="fallback"&&mode!="correlation")||distance<3||workers<1||workers>64||!max_shots||!target)
            throw std::invalid_argument("Invalid mode, distance, workers, shots, or target");
        for(const auto& suffix:{".txt","_hw.csv","_fault_hw.csv","_joint_hw.csv","_failures.csv","_thresholds.csv"})
            if(std::filesystem::exists(prefix+suffix))throw std::runtime_error("Output already exists: "+prefix+suffix);
        const auto circuit=sc_si1000(distance,distance,1e-3,false,true);
        const auto dem=stim::circuit_to_dem(circuit,{.decompose_errors=true});
        const size_t nd=dem.count_detectors(),ne=dem.count_errors();
        const size_t batch=1024;
        std::atomic<uint64_t> ticket{0},processed{0},original_errors{0},alternative_errors{0};
        std::atomic<int> live{workers};std::atomic<bool> failed{false};
        std::mutex mutex;std::string error;
        std::vector<Worker> results(workers);std::vector<std::thread> threads;
        const auto started=Clock::now();
        for(int w=0;w<workers;++w)threads.emplace_back([&,w]{
            try {
                dec::LibraVariant::Config cfg;cfg.seed=ensemble_seed;
                cfg.large_gap_fallback=mode=="fallback"?true:keep_fallback;
                dec::LibraVariant original(dem,cfg);
                cfg.large_gap_fallback=mode=="fallback"?false:keep_fallback;
                cfg.correlated_matching=mode!="correlation";
                dec::LibraVariant alternative(dem,cfg);
                dec::PyMatching pm(dem,false);
                original.prepare_ensemble();alternative.prepare_ensemble();
                std::mt19937_64 rng(seed+0x9e3779b97f4a7c15ULL*uint64_t(w+1));
                stim::DemSampler<64> sampler(dem,std::move(rng),batch);
                auto& r=results[w];uint64_t local_shot=0;
                while(!failed.load()) {
                    if(original_errors.load()>=target&&alternative_errors.load()>=target)break;
                    const auto start=ticket.fetch_add(batch);if(start>=max_shots)break;
                    const size_t take=std::min<uint64_t>(batch,max_shots-start);
                    sampler.resample(false);
                    const auto dets=sampler.det_buffer.transposed(),obs=sampler.obs_buffer.transposed();
                    const auto faults=sampler.err_buffer.transposed();
                    for(size_t i=0;i<take&&!failed.load();++i,++local_shot){
                        const size_t dhw=hw(dets[i],nd),fhw=hw(faults[i],ne);
                        Counts c;c.shots=1;
                        dec::LibraVariant::DetailedResult a,b;
                        // Alternate call order to reduce systematic timing bias.
                        if(local_shot&1){auto t=Clock::now();b=alternative.decode_detailed(dets[i]);c.alternative_ns=elapsed(t);
                            t=Clock::now();a=original.decode_detailed(dets[i]);c.original_ns=elapsed(t);}
                        else {auto t=Clock::now();a=original.decode_detailed(dets[i]);c.original_ns=elapsed(t);
                            t=Clock::now();b=alternative.decode_detailed(dets[i]);c.alternative_ns=elapsed(t);}
                        auto t=Clock::now();const auto p=pm.decode(dets[i],obs[i]);c.pm_ns=elapsed(t);
                        c.original=a.prediction!=obs[i][0];c.alternative=b.prediction!=obs[i][0];c.pm=p.flipped_obs[0]!=obs[i][0];
                        c.both=c.original&&c.alternative;c.original_only=c.original&&!c.alternative;
                        c.alternative_only=c.alternative&&!c.original;
                        c.original_ensemble=a.ensemble_used;c.alternative_ensemble=b.ensemble_used;
                        if(mode=="fallback") {
                            if(a.gap_db!=b.gap_db)throw std::runtime_error("Fallback changed initial gap");
                            if(a.gap_db<20&&(a.prediction!=b.prediction||a.errors!=b.errors||a.log_masses[0]!=b.log_masses[0]||a.log_masses[1]!=b.log_masses[1]))
                                throw std::runtime_error("Fallback ablation changed small-gap result");
                            if(a.gap_db>=20&&b.prediction!=b.initial_prediction)
                                throw std::runtime_error("No-fallback did not choose cheaper constrained result");
                        }
                        r.total+=c;r.syndrome[dhw]+=c;r.faults[fhw]+=c;r.joint[{dhw,fhw}]+=c;
                        if(c.original||c.alternative||c.pm){
                            std::ostringstream f;f<<std::setprecision(17)<<w<<','<<local_shot<<','<<dhw<<','<<fhw<<','<<bool(obs[i][0])
                                <<','<<a.prediction<<','<<b.prediction<<','<<bool(p.flipped_obs[0])<<','<<a.gap_db<<','<<b.gap_db
                                <<','<<a.ensemble_used<<','<<b.ensemble_used<<','<<a.initial_prediction<<','<<b.initial_prediction
                                <<','<<indices(dets[i],nd)<<','<<indices(faults[i],ne);
                            r.failures.push_back(f.str());
                        }
                        processed.fetch_add(1);original_errors.fetch_add(c.original);alternative_errors.fetch_add(c.alternative);
                        if(original_errors.load()>=target&&alternative_errors.load()>=target)break;
                    }
                }
            }catch(const std::exception& e){std::lock_guard<std::mutex> lock(mutex);if(error.empty())error=e.what();failed.store(true);}
            live.fetch_sub(1);
        });
        auto last=started;
        while(live.load()){
            std::this_thread::sleep_for(std::chrono::seconds(1));const auto now=Clock::now();
            if(now-last>=std::chrono::seconds(10)){
                std::cerr<<"mode="<<mode<<" d="<<distance<<" shots="<<processed.load()<<" original_errors="<<original_errors.load()
                    <<" alternative_errors="<<alternative_errors.load()<<" live="<<live.load()<<" elapsed_s="<<std::chrono::duration<double>(now-started).count()<<'\n';last=now;
            }
        }
        for(auto& t:threads)t.join();if(failed.load())throw std::runtime_error(error);
        Worker all;
        for(auto& r:results){all.total+=r.total;for(auto& [k,c]:r.syndrome)all.syndrome[k]+=c;
            for(auto& [k,c]:r.faults)all.faults[k]+=c;for(auto& [k,c]:r.joint)all.joint[k]+=c;}
        auto out=create(prefix+".txt");const auto& c=all.total;
        out<<std::setprecision(17)<<"mode="<<mode<<"\ndistance="<<distance<<"\nrounds="<<distance
            <<"\nphysical_error_rate=0.001\ninclude_opposite_basis_detectors=true\nmemory_basis=Z\nworkers="<<workers
            <<"\nseed="<<seed<<"\nensemble_seed="<<ensemble_seed<<"\nensemble_size=60\ngap_threshold_db=20\nlarge_gap_fallback="<<(mode=="fallback"?1:keep_fallback)
            <<"\nerror_target="<<target<<"\nmaximum_shots="<<max_shots<<"\nshots="<<c.shots
            <<"\noriginal_errors="<<c.original<<"\nalternative_errors="<<c.alternative<<"\npymatching_errors="<<c.pm
            <<"\nboth_libra_errors="<<c.both<<"\noriginal_only_errors="<<c.original_only<<"\nalternative_only_errors="<<c.alternative_only
            <<"\noriginal_ler="<<double(c.original)/c.shots<<"\nalternative_ler="<<double(c.alternative)/c.shots
            <<"\npymatching_ler="<<double(c.pm)/c.shots<<"\noriginal_ensemble_shots="<<c.original_ensemble
            <<"\nalternative_ensemble_shots="<<c.alternative_ensemble<<"\noriginal_mean_decode_ns="<<double(c.original_ns)/c.shots
            <<"\nalternative_mean_decode_ns="<<double(c.alternative_ns)/c.shots<<"\npymatching_mean_decode_ns="<<double(c.pm_ns)/c.shots
            <<"\nwall_time_seconds="<<std::chrono::duration<double>(Clock::now()-started).count()<<'\n';
        auto h=create(prefix+"_hw.csv");csv_header(h,"detector_hamming_weight");
        for(auto& [k,c]:all.syndrome){h<<k;csv_counts(h,c);}
        auto f=create(prefix+"_fault_hw.csv");csv_header(f,"dem_fault_hamming_weight");
        for(auto& [k,c]:all.faults){f<<k;csv_counts(f,c);}
        auto j=create(prefix+"_joint_hw.csv");csv_header(j,"detector_hamming_weight,dem_fault_hamming_weight");
        for(auto& [k,c]:all.joint){j<<k.first<<','<<k.second;csv_counts(j,c);}
        auto failures=create(prefix+"_failures.csv");
        failures<<"worker,local_shot,detector_hw,dem_fault_hw,truth,original_prediction,alternative_prediction,pymatching_prediction,original_gap_db,alternative_gap_db,original_ensemble,alternative_ensemble,original_initial_prediction,alternative_initial_prediction,detector_ids,dem_error_ids\n";
        for(auto& r:results)for(auto& s:r.failures)failures<<s<<'\n';
        auto thresholds=create(prefix+"_thresholds.csv");
        thresholds<<"threshold,alternative_shots_below,alternative_fraction,original_errors_below,alternative_errors_below,pymatching_errors_below,original_only_errors_below,alternative_only_errors_below,original_ler_below,alternative_ler_below,pymatching_ler_below,adaptive_total_errors,adaptive_ler,mean_adaptive_decode_ns\n";
        thresholds<<std::setprecision(17);Counts below;
        for(auto& [k,c]:all.syndrome){below+=c;
            const auto errors=below.alternative+all.total.original-below.original;
            const auto time=below.alternative_ns+all.total.original_ns-below.original_ns;
            thresholds<<k<<','<<below.shots<<','<<double(below.shots)/all.total.shots<<','<<below.original<<','<<below.alternative
                <<','<<below.pm<<','<<below.original_only<<','<<below.alternative_only
                <<','<<double(below.original)/below.shots<<','<<double(below.alternative)/below.shots<<','<<double(below.pm)/below.shots
                <<','<<errors<<','<<double(errors)/all.total.shots<<','<<double(time)/all.total.shots<<'\n';
        }
        std::cout<<"completed "<<prefix<<" shots="<<c.shots<<" original="<<c.original<<" alternative="<<c.alternative<<" pm="<<c.pm<<'\n';
        return 0;
    }catch(const std::exception& e){std::cerr<<"error: "<<e.what()<<'\n';return 1;}
}
