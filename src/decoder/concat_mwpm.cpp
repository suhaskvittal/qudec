/*
 *  author: OpenAI Codex
 *  date:   22 October 2026
 * */
#include "decoder/concat_mwpm.h"
#include "decoder/concat_mwpm/circuit.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <bit>
#include <climits>
#include <sstream>
#include <tuple>
#include <ostream>
#include <stdexcept>

#if defined(ENABLE_MPI)
#include <mpi.h>
#endif

namespace dec
{

ConcatMWPM::ConcatMWPM(const stim::DetectorErrorModel& dem)
    :ConcatMWPM(dem,config_type{})
{
}

ConcatMWPM::ConcatMWPM(const stim::DetectorErrorModel& dem, config_type conf)
    :num_detectors(dem.count_detectors()), num_observables(dem.count_observables()),
     decoder_(std::make_unique<concat::CircuitDecoder>(dem,conf.memory_x,
                                                       conf.use_flags,conf.flag_cache_capacity,
                                                       conf.boundary_gap,concat::EnsembleConfig{conf.ensemble_size,conf.ensemble_seed},
                                                       concat::FlagRenormalization{conf.flag_renormalization,conf.physical_error_rate},
                                                       conf.gap_gating,conf.gap_gating_threshold_db,conf.belief_propagation)),
     debug_(conf.debug),belief_propagation_(conf.belief_propagation),decoding_cores_(conf.decoding_cores),gap_gating_(conf.gap_gating),
     ensemble_size_(conf.ensemble_size),gap_gating_threshold_db_(conf.gap_gating_threshold_db)
{
    if (!decoding_cores_ || decoding_cores_>INT_MAX || num_detectors>INT_MAX)
        throw std::invalid_argument("Invalid concatenated MWPM decoding core count or detector count");
    if (!std::isfinite(conf.gap_gating_threshold_db) || conf.gap_gating_threshold_db<0)
        throw std::invalid_argument("Concat gap-gating threshold must be finite and nonnegative");
    if (gap_gating_ && conf.boundary_gap)
        throw std::invalid_argument("--concat-gap-gating cannot be combined with --concat-boundary-gap");
    if (gap_gating_ && decoding_cores_>1)
        throw std::invalid_argument("--concat-gap-gating does not support grouped MPI decoding");
    if (decoding_cores_>1)
    {
#if defined(ENABLE_MPI)
        int initialized=0,rank=0,size=1;
        MPI_Initialized(&initialized);
        if (!initialized) throw std::invalid_argument("Parallel ConcatMWPM requires initialized MPI");
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
        if (size%decoding_cores_) throw std::invalid_argument("MPI world size must be divisible by decoding_core_count()");
        MPI_Comm_split(MPI_COMM_WORLD,rank/decoding_cores_,rank,&decoding_comm_);
        MPI_Comm_rank(decoding_comm_,&decoding_rank_);
#else
        throw std::invalid_argument("Parallel ConcatMWPM requires an ENABLE_MPI build");
#endif
    }
}

ConcatMWPM::~ConcatMWPM()
{
#if defined(ENABLE_MPI)
    int finalized=0;MPI_Finalized(&finalized);
    if (!finalized && decoding_comm_!=MPI_COMM_NULL) MPI_Comm_free(&decoding_comm_);
#endif
}

#if defined(ENABLE_MPI)
namespace
{
// Make every group member leave a failing phase together, instead of leaving
// helpers blocked in the following collective.
void collective_error(const std::string& local,MPI_Comm comm)
{
    int rank=0;MPI_Comm_rank(comm,&rank);
    int failing=local.empty() ? INT_MAX:rank;
    MPI_Allreduce(MPI_IN_PLACE,&failing,1,MPI_INT,MPI_MIN,comm);
    if (failing==INT_MAX) return;
    int length=rank==failing ? std::min<size_t>(local.size(),INT_MAX):0;
    MPI_Bcast(&length,1,MPI_INT,failing,comm);
    std::string message=rank==failing ? local.substr(0,length):std::string(length,' ');
    MPI_Bcast(message.data(),length,MPI_CHAR,failing,comm);
    throw std::runtime_error("Parallel ConcatMWPM worker "+std::to_string(failing)+": "+message);
}
std::vector<uint64_t> gather_words(const std::vector<uint64_t>& local,MPI_Comm comm)
{
    collective_error(local.size()>INT_MAX ? "MPI payload exceeds INT_MAX":"",comm);
    int size=1;MPI_Comm_size(comm,&size);
    std::vector<int> sizes(size),offsets(size);
    int count=local.size();MPI_Allgather(&count,1,MPI_INT,sizes.data(),1,MPI_INT,comm);
    int64_t total=0;
    for (int r=0;r<size;r++)
    {
        if (total>INT_MAX) throw std::overflow_error("Combined MPI payload exceeds INT_MAX");
        offsets[r]=total;total+=sizes[r];
    }
    if (total>INT_MAX) throw std::overflow_error("Combined MPI payload exceeds INT_MAX");
    std::vector<uint64_t> all(total);
    MPI_Allgatherv(local.data(),count,MPI_UINT64_T,all.data(),sizes.data(),offsets.data(),MPI_UINT64_T,comm);
    return all;
}
template<class T> void put_vector(std::vector<uint64_t>& out,const std::vector<T>& v)
{
    out.push_back(v.size());for (auto x:v) out.push_back(static_cast<uint64_t>(x));
}
template<class T> void get_vector(const std::vector<uint64_t>& in,size_t& at,std::vector<T>& v)
{
    const auto count=in.at(at++);v.reserve(count);
    for (size_t j=0;j<count;j++) v.push_back(static_cast<T>(in.at(at++)));
}
template<class T> void put_paths(std::vector<uint64_t>& out,const std::vector<std::vector<T>>& paths)
{
    out.push_back(paths.size());for (const auto& p:paths) put_vector(out,p);
}
template<class T> void get_paths(const std::vector<uint64_t>& in,size_t& at,std::vector<std::vector<T>>& paths)
{
    const auto count=in.at(at++);paths.resize(count);
    for (auto& p:paths) get_vector(in,at,p);
}
}

uint64_t ConcatMWPM::decode_parallel(std::vector<uint8_t>& bits,std::string* debug)
{
    MPI_Bcast(bits.data(),bits.size(),MPI_UINT8_T,0,decoding_comm_);
    std::vector<concat::FirstStageResult> local_first;
    std::string error;
    try {local_first=decoder_->decode_first_stage(bits,decoding_rank_,decoding_cores_,debug!=nullptr);}
    catch (const std::exception& e) {error=e.what();}
    collective_error(error,decoding_comm_);
    std::vector<uint64_t> words;
    for (const auto& first:local_first)
    {
        words.insert(words.end(),{uint64_t(first.color),uint64_t(first.boundary_mask),first.member});
        put_vector(words,first.selected);put_paths(words,first.paths);put_paths(words,first.nodes);
    }
    const auto first_words=gather_words(words,decoding_comm_);
    std::vector<concat::FirstStageResult> first;
    for (size_t at=0;at<first_words.size();)
    {
        concat::FirstStageResult v{int(first_words.at(at)),int(first_words.at(at+1)),first_words.at(at+2),{},{},{}};at+=3;
        get_vector(first_words,at,v.selected);get_paths(first_words,at,v.paths);get_paths(first_words,at,v.nodes);
        first.push_back(std::move(v));
    }
    auto order=[](const auto& v){return std::tuple{v.member,v.color,v.boundary_mask};};
    std::sort(first.begin(),first.end(),[&](const auto& a,const auto& b){return order(a)<order(b);});
    std::vector<concat::Result> candidates;
    std::string trace;
    error.clear();
    try {candidates=decoder_->decode_second_stage(bits,first,decoding_rank_,decoding_cores_,debug ? &trace:nullptr);}
    catch (const std::exception& e) {error=e.what();}
    collective_error(error,decoding_comm_);
    words.clear();
    for (const auto& v:candidates)
    {
        words.insert(words.end(),{v.observable,std::bit_cast<uint64_t>(v.weight),v.member,uint64_t(v.color),uint64_t(v.boundary_mask)});
        put_vector(words,v.events);
    }
    const auto final_words=gather_words(words,decoding_comm_);
    candidates.clear();
    for (size_t at=0;at<final_words.size();)
    {
        concat::Result v{final_words.at(at),std::bit_cast<double>(final_words.at(at+1)),int(final_words.at(at+3)),{},{},{}};
        v.member=final_words.at(at+2);v.boundary_mask=final_words.at(at+4);at+=5;
        get_vector(final_words,at,v.events);candidates.push_back(std::move(v));
    }
    if (candidates.empty()) throw std::runtime_error("No feasible parallel concatenated matching candidate");
    std::sort(candidates.begin(),candidates.end(),[&](const auto& a,const auto& b){return order(a)<order(b);});
    auto best=candidates.begin();
    for (auto it=candidates.begin();it!=candidates.end();++it) if (it->weight<best->weight-1e-9) best=it;
    parallel_candidate_gap_.reset();
    for (const auto& candidate:candidates) if (candidate.observable!=best->observable)
    {
        const double gap=std::max(0.0,candidate.weight-best->weight);
        if (!parallel_candidate_gap_ || gap<*parallel_candidate_gap_) parallel_candidate_gap_=gap;
    }
    if (debug)
    {
        words.assign(trace.begin(),trace.end());
        const auto traces=gather_words(words,decoding_comm_);
        if (!decoding_rank_)
        {
            debug->assign(traces.begin(),traces.end());
            std::ostringstream winner;
            winner<<"Winner: !"<<"RGB"[best->color]<<" member="<<best->member<<" boundary_mask="<<best->boundary_mask
                  <<" weight="<<best->weight<<" observable="<<best->observable<<" correction:";
            for (auto event:best->events) winner<<" E"<<event;
            *debug+=winner.str()+"\n";
        }
    }
    return decoding_rank_ ? 0:best->observable;
}
#endif

Result
ConcatMWPM::decode(SyndromeRef syndrome, ObsRef truth)
{
    if (decoding_cores_==1 && syndrome.num_bits_padded() < num_detectors)
        throw std::invalid_argument("Concatenated MWPM syndrome is too short");
    const auto start = std::chrono::steady_clock::now();
    std::vector<uint8_t> bits(num_detectors);
    if (!decoding_rank_ && syndrome.num_bits_padded()>=num_detectors)
        for (size_t v = 0; v < num_detectors; v++) bits[v] = syndrome[v];
    Result result;
    uint64_t frame=0;
    bool evaluated_gap=false,expanded=false;
    double gap_db=0;
#if defined(ENABLE_MPI)
    if (decoding_cores_>1)
    {
        collective_error(!decoding_rank_ && syndrome.num_bits_padded()<num_detectors ? "Syndrome is too short":"",decoding_comm_);
        frame=decode_parallel(bits,debug_ ? &result.debug_info:nullptr);
    }
    else
#endif
    {
        if (gap_gating_ && ensemble_size_>1)
        {
            const auto gated=decoder_->decode_gated(bits,gap_gating_threshold_db_,debug_ ? &result.debug_info:nullptr);
            frame=gated.result.observable;gap_db=gated.gap_db;expanded=gated.used_ensemble;evaluated_gap=true;
        }
        else frame=decoder_->decode(bits,debug_ ? &result.debug_info:nullptr);
    }
    const auto candidate_gap=decoding_cores_>1 ? parallel_candidate_gap_:decoder_->candidate_gap();
    if (belief_propagation_ && candidate_gap)
    {
        result.has_candidate_gap=true;
        result.candidate_gap_db=*candidate_gap*(10/std::log(10.0));
    }
    result.flipped_obs = ObsType(num_observables);
    result.cluster_gap_db=gap_db;
    result.has_cluster_gap=evaluated_gap;
    result.flipped_obs[0] = frame & 1;
    result.matching_data.frame_flips = result.flipped_obs;
    result.decode_latency_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now()-start).count();
    if (!decoding_rank_)
    {
        if (belief_propagation_)
        {
            if (result.has_candidate_gap)
            {
                stats_.candidate_gap_available++;
                const double gap=result.candidate_gap_db;
                s_unsigned_gap.add(gap);
                const bool error=truth.num_bits_padded() && truth[0]!=result.flipped_obs[0];
                if (error) s_unsigned_gap_errors.add(gap);
                if (truth.num_bits_padded()) s_signed_gap.add(error ? -gap:gap);
                result.matching_data.gap=error ? -gap:gap;
            }
            else stats_.candidate_gap_unavailable++;
            if (debug_) result.debug_info+=result.has_candidate_gap ?
                "Candidate frame gap (dB): "+std::to_string(result.candidate_gap_db)+"\n":
                "Candidate frame gap unavailable: no opposing-frame candidate found.\n";
        }
        stats_.shots++;
        stats_.decode_total_ns += result.decode_latency_ns;
        stats_.decode_min_ns = std::min(stats_.decode_min_ns,result.decode_latency_ns);
        stats_.decode_max_ns = std::max(stats_.decode_max_ns,result.decode_latency_ns);
        if (evaluated_gap) {stats_.gap_evaluations++;stats_.gap_expansions+=expanded;}
    }
    return result;
}

ConcatMWPM::Statistics
ConcatMWPM::statistics() const
{
    auto out = stats_;
    if (!out.shots) out.decode_min_ns = 0;
    return out;
}

void
ConcatMWPM::print_stats(std::ostream& out) const
{
    const auto stats = statistics();
    out << "ConcatMWPM shots=" << stats.shots
        << " average_decode_us=" << (stats.shots ? double(stats.decode_total_ns)/stats.shots/1000 : 0)
        << " min_decode_us=" << double(stats.decode_min_ns)/1000
        << " max_decode_us=" << double(stats.decode_max_ns)/1000
        << " gap_evaluations=" << stats.gap_evaluations
        << " gap_expansions=" << stats.gap_expansions << '\n';
    /*
    if (belief_propagation_)
        out << "BP candidate_frame_gap available=" << stats.candidate_gap_available
            << " unavailable=" << stats.candidate_gap_unavailable << '\n'
            << s_signed_gap.to_string_full() << '\n'
            << s_unsigned_gap.to_string_full() << '\n'
            << s_unsigned_gap_errors.to_string_full() << '\n';
    */
}

void
ConcatMWPM::mpi_accumulate()
{
#if defined(ENABLE_MPI)
    uint64_t sums[]{stats_.shots,stats_.decode_total_ns};
    MPI_Allreduce(MPI_IN_PLACE,sums,2,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
    stats_.shots = sums[0];
    stats_.decode_total_ns = sums[1];
    MPI_Allreduce(MPI_IN_PLACE,&stats_.decode_min_ns,1,MPI_UINT64_T,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&stats_.decode_max_ns,1,MPI_UINT64_T,MPI_MAX,MPI_COMM_WORLD);
    uint64_t gates[]{stats_.gap_evaluations,stats_.gap_expansions};
    MPI_Allreduce(MPI_IN_PLACE,gates,2,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
    stats_.gap_evaluations=gates[0];stats_.gap_expansions=gates[1];
    uint64_t confidence[]{stats_.candidate_gap_available,stats_.candidate_gap_unavailable};
    MPI_Allreduce(MPI_IN_PLACE,confidence,2,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
    stats_.candidate_gap_available=confidence[0];stats_.candidate_gap_unavailable=confidence[1];
    s_signed_gap.mpi_accumulate();s_unsigned_gap.mpi_accumulate();s_unsigned_gap_errors.mpi_accumulate();
#endif
}
} // namespace dec
