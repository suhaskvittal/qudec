// Libra matching/fallback variants. Default options reproduce the preserved
// PyMatching auxiliary-assignment implementation exactly. Separate class and
// files allow paired experiments without changing existing project files.
#include "decoder/libra_variant.h"
#include <pymatching/sparse_blossom/driver/mwpm_decoding.h>
#include <pymatching/sparse_blossom/driver/user_graph.h>
#include "decoder/libra_assignment_pymatching.h"
#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <unordered_map>

namespace dec {
namespace {
using Ids = std::vector<uint32_t>;
using Clock = std::chrono::steady_clock;
uint64_t ns(Clock::time_point t) {return std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now()-t).count();}
struct Dsu {
    std::vector<uint32_t> p;
    explicit Dsu(size_t n):p(n){std::iota(p.begin(),p.end(),0);}
    uint32_t root(uint32_t x){return p[x]==x?x:p[x]=root(p[x]);}
    void join(uint32_t a,uint32_t b){a=root(a);b=root(b);if(a!=b)p[b]=a;}
};
Ids sx(const Ids& a,const Ids& b){Ids c;c.reserve(a.size()+b.size());std::set_symmetric_difference(a.begin(),a.end(),b.begin(),b.end(),std::back_inserter(c));return c;}
Ids parity(Ids a){std::sort(a.begin(),a.end());Ids b;for(size_t i=0;i<a.size();){size_t j=i+1;while(j<a.size()&&a[j]==a[i])++j;if((j-i)&1)b.push_back(a[i]);i=j;}return b;}
uint64_t key(uint32_t a,uint32_t b){if(a>b)std::swap(a,b);return (uint64_t(a)<<32)|b;}
double logadd(double a,double b){if(a<b)std::swap(a,b);return a+std::log1p(std::exp(b-a));}
}

struct LibraVariant::Impl {
    struct Hyperedge {Ids dets; bool obs; double p,w; Ids components;};
    struct GraphEdge {uint32_t a,b; bool obs; uint32_t single=UINT32_MAX; Ids pairs;};
    struct Candidate {Ids errors; bool obs=false; double weight=0;};
    struct Cycle {Ids errors; double delta; bool obs=false;};
    struct Member {std::unique_ptr<pm::Mwpm> graph;std::vector<double> weights;};
    Config config;
    size_t nd;
    std::vector<Hyperedge> hyper;
    std::vector<GraphEdge> gedges;
    std::unordered_map<uint64_t,uint32_t> edge_ids;
    stim::DetectorErrorModel canonical;
    std::unique_ptr<pm::Mwpm> complement;
    std::unique_ptr<pm::Mwpm> original;
    std::vector<Member> members;
    Statistics stats;
    size_t rules=0;
    std::vector<int32_t> edge_position,det_owner;
    std::vector<uint32_t> touched_dets;

    Impl(const stim::DetectorErrorModel& dem,Config conf):config(conf),nd(dem.count_detectors()) {
        if(dem.count_observables()!=1)throw std::invalid_argument("LibraVariant requires exactly one logical observable");
        if(nd>=UINT32_MAX-1)throw std::invalid_argument("LibraVariant detector index overflow");
        if(!std::isfinite(conf.gap_threshold_db)||conf.gap_threshold_db<0||conf.positive_cycle_capacity>30)
            throw std::invalid_argument("Invalid LibraVariant gap/cycle capacity");
        Dsu basis(nd);uint32_t logical_basis=UINT32_MAX;
        dem.iter_flatten_error_instructions([&](const auto& x){x.for_separated_targets([&](const auto& g){Ids ds;for(auto t:g)if(t.is_relative_detector_id())ds.push_back(t.val());if(ds.size()>2)throw std::invalid_argument("LibraVariant requires a decomposed DEM");if(ds.size()==2)basis.join(ds[0],ds[1]);});});
        dem.iter_flatten_error_instructions([&](const auto& x){x.for_separated_targets([&](const auto& g){Ids ds;bool o=false;for(auto t:g){if(t.is_relative_detector_id())ds.push_back(t.val());else if(t.is_observable_id())o=!o;}if(o){if(ds.size()!=1)throw std::invalid_argument("LibraVariant complementary matching requires observables on boundary edges");uint32_t b=basis.root(ds[0]);if(logical_basis!=UINT32_MAX&&logical_basis!=b)throw std::invalid_argument("Logical observable crosses multiple CSS components");logical_basis=b;}});});
        if(logical_basis==UINT32_MAX)throw std::invalid_argument("Missing observable-crossing errors");
        // Stim separator choices can include redundant pieces. Recover each
        // channel's total parity, then split it by the two CSS graph components.
        std::map<std::pair<Ids,bool>,double> merged;
        dem.iter_flatten_error_instructions([&](const auto& x){Ids ds;bool o=false;for(auto t:x.target_data){if(t.is_relative_detector_id())ds.push_back(t.val());else if(t.is_observable_id())o=!o;}ds=parity(std::move(ds));double p=x.arg_data[0];if(p==0)return;if(!(p>0&&p<0.5))throw std::invalid_argument("LibraVariant requires 0 < error probabilities < 1/2");if(ds.empty()){if(o)throw std::invalid_argument("Undetectable logical channel is unsupported");return;}auto& v=merged[{ds,o}];v=v+p-2*v*p;});
        for(auto& [pat,p]:merged){
            Hyperedge h{pat.first,pat.second,p,std::log((1-p)/p),{}};
            std::map<uint32_t,Ids> parts;for(uint32_t d:h.dets)parts[basis.root(d)].push_back(d);
            if(parts.size()>2)throw std::invalid_argument("Hyperedge spans more than two CSS components");
            std::vector<stim::DemTarget> ts;
            for(auto& [b,ds]:parts){
                if(ds.size()>2)throw std::invalid_argument("Error has more than two detectors in one CSS basis");
                bool obs=h.obs&&b==logical_basis;
                if(obs&&ds.size()!=1)throw std::invalid_argument("Canonical observable component is not a boundary edge");
                uint32_t a=ds[0],z=ds.size()==2?ds[1]:UINT32_MAX;
                auto [it,inserted]=edge_ids.emplace(key(a,z),gedges.size());
                if(inserted)gedges.push_back({a,z,obs,UINT32_MAX,{}});
                else if(gedges[it->second].obs!=obs)throw std::invalid_argument("Parallel graph edges have different logical actions");
                h.components.push_back(it->second);
                if(!ts.empty())ts.push_back(stim::DemTarget::separator());
                for(uint32_t d:ds)ts.push_back(stim::DemTarget::relative_detector_id(d));
                if(obs)ts.push_back(stim::DemTarget::observable_id(0));
            }
            if(h.obs&&!parts.count(logical_basis))throw std::invalid_argument("Logical error has no logical-basis component");
            canonical.append_error_instruction(p,ts,"");
            uint32_t id=hyper.size();hyper.push_back(std::move(h));
            if(hyper.back().components.size()==1){auto& e=gedges[hyper.back().components[0]];if(e.single==UINT32_MAX||hyper[id].w<hyper[e.single].w)e.single=id;}
            else for(uint32_t e:hyper.back().components)gedges[e].pairs.push_back(id);
        }
        for(const auto& e:gedges)if(e.single==UINT32_MAX)throw std::invalid_argument("A CSS component has no independent single-edge channel");
        // Turn the logical boundary into an explicit detector to constrain parity.
        stim::DetectorErrorModel folded;
        canonical.iter_flatten_error_instructions([&](const auto& x){std::vector<stim::DemTarget> ts;for(auto t:x.target_data)ts.push_back(t.is_observable_id()?stim::DemTarget::relative_detector_id(nd):t);folded.append_error_instruction(x.arg_data[0],ts,"");});
        complement=std::make_unique<pm::Mwpm>(pm::detector_error_model_to_mwpm(folded,pm::NUM_DISTINCT_WEIGHTS,true,config.correlated_matching));
        if(config.large_gap_fallback)original=std::make_unique<pm::Mwpm>(pm::detector_error_model_to_mwpm(canonical,pm::NUM_DISTINCT_WEIGHTS,true,config.correlated_matching));
        for(const auto& node:complement->flooder.graph.nodes)
            for(const auto& implied:node.neighbor_implied_weights)rules+=implied.size();
        members.resize(config.ensemble_size);edge_position.resize(gedges.size(),-1);det_owner.resize(nd,-1);
    }

    void prepare(size_t i) {
        auto& m=members[i];if(m.graph)return;
        std::mt19937_64 rng(config.seed+0x9e3779b97f4a7c15ULL*(i+1));
        std::normal_distribution<double> normal(0,(i&1)?std::log(4.):std::log(2.));
        stim::DetectorErrorModel perturbed;size_t h=0;
        canonical.iter_flatten_error_instructions([&](const auto& x){double p=std::min(0.499999999,x.arg_data[0]*std::exp(normal(rng)));p=std::max(p,std::numeric_limits<double>::min());perturbed.append_error_instruction(p,x.target_data,"");m.weights.push_back(std::log((1-p)/p));++h;});
        m.graph=std::make_unique<pm::Mwpm>(pm::detector_error_model_to_mwpm(perturbed,pm::NUM_DISTINCT_WEIGHTS,true,config.correlated_matching));
    }

    Candidate candidate(pm::Mwpm& graph,const std::vector<uint64_t>& events,bool folded,const std::vector<double>* ws=nullptr) {
        std::vector<int64_t> es;
        if(config.correlated_matching)
            pm::decode_detection_events_to_edges_with_edge_correlations(graph,events,es);
        else
            pm::decode_detection_events_to_edges(graph,events,es);
        Ids chosen;chosen.reserve(es.size()/2);
        for(size_t i=0;i<es.size();i+=2){uint32_t a=es[i],b=es[i+1]==-1?UINT32_MAX:uint32_t(es[i+1]);if(folded){if(a==nd)a=UINT32_MAX;if(b==nd)b=UINT32_MAX;}auto it=edge_ids.find(key(a,b));if(it==edge_ids.end())throw std::runtime_error("Matched edge missing from original hypergraph");chosen.push_back(it->second);}
        chosen=parity(std::move(chosen));
        return assign(chosen,ws);
    }

    Candidate assign(const Ids& chosen,const std::vector<double>* ws) {
        auto w=[&](uint32_t h){return ws?(*ws)[h]:hyper[h].w;};
        for(size_t i=0;i<chosen.size();++i)edge_position[chosen[i]]=i;
        struct Pair {uint32_t a,b,h;};std::vector<Pair> pairs;Dsu components(chosen.size());
        Ids degree(chosen.size());
        for(size_t i=0;i<chosen.size();++i)for(uint32_t h:gedges[chosen[i]].pairs){const auto& cs=hyper[h].components;uint32_t other=cs[0]^cs[1]^chosen[i];int32_t j=edge_position[other];
            // Replacing two independent singles by a more expensive pair cannot
            // be part of an optimal assignment. Keep only profitable pairings.
            if(j>int32_t(i)&&w(h)<w(gedges[chosen[i]].single)+w(gedges[other].single)){
                pairs.push_back({uint32_t(i),uint32_t(j),h});components.join(i,j);++degree[i];++degree[j];
            }
        }
        Candidate out;std::map<uint32_t,Ids> blocks;
        for(size_t i=0;i<chosen.size();++i){if(degree[i])blocks[components.root(i)].push_back(i);else out.errors.push_back(gedges[chosen[i]].single);}
        std::vector<int32_t> loc(chosen.size(),-1);
        for(auto& [root,vs]:blocks){
            size_t n=vs.size();
            if(n==1){out.errors.push_back(gedges[chosen[vs[0]]].single);continue;}
            if(n==2){uint32_t sa=gedges[chosen[vs[0]]].single,sb=gedges[chosen[vs[1]]].single;double cost=w(sa)+w(sb);uint32_t best=UINT32_MAX;for(auto p:pairs)if(components.root(p.a)==root&&w(p.h)<cost){cost=w(p.h);best=p.h;}if(best==UINT32_MAX){out.errors.push_back(sa);out.errors.push_back(sb);}else out.errors.push_back(best);continue;}
            for(size_t j=0;j<n;++j)loc[vs[j]]=j;
            std::vector<double> singles;singles.reserve(n);
            for(size_t j=0;j<n;++j)singles.push_back(w(gedges[chosen[vs[j]]].single));
            std::vector<libra_assignment::Pair> assignment_pairs;
            std::vector<uint32_t> channels;
            for(auto p:pairs)if(components.root(p.a)==root){
                assignment_pairs.push_back({uint32_t(loc[p.a]),uint32_t(loc[p.b]),w(p.h)});
                channels.push_back(p.h);
            }
            const auto assignment=libra_assignment::solve(singles,assignment_pairs);
            for(size_t j=0;j<n;++j)if(assignment.mate[j]<0)
                out.errors.push_back(gedges[chosen[vs[j]]].single);
            for(size_t i:assignment.pair_indices)out.errors.push_back(channels[i]);
        }
        for(uint32_t e:chosen)edge_position[e]=-1;
        out.errors=parity(std::move(out.errors));for(uint32_t h:out.errors){out.weight+=hyper[h].w;out.obs^=hyper[h].obs;}
        return out;
    }

    double relative(const Ids& cy,const Candidate& base) const {double delta=0;for(uint32_t h:cy)delta+=(std::binary_search(base.errors.begin(),base.errors.end(),h)?-1:1)*hyper[h].w;return delta;}
    std::vector<Cycle> difference(const Candidate& base,const Candidate& other) {
        Ids diff=sx(base.errors,other.errors);Dsu dsu(diff.size());touched_dets.clear();
        for(size_t i=0;i<diff.size();++i)for(uint32_t d:hyper[diff[i]].dets){if(det_owner[d]<0){det_owner[d]=i;touched_dets.push_back(d);}else dsu.join(i,det_owner[d]);}
        for(uint32_t d:touched_dets)det_owner[d]=-1;
        std::map<uint32_t,Cycle> cs;
        for(size_t i=0;i<diff.size();++i){auto& c=cs[dsu.root(i)];uint32_t h=diff[i];c.errors.push_back(h);c.obs^=hyper[h].obs;}
        std::vector<Cycle> out;for(auto& [r,c]:cs){c.delta=relative(c.errors,base);out.push_back(std::move(c));}return out;
    }
    void save_cycle(std::vector<Cycle>& heap,Cycle c) {
        if(!config.degeneracy||!config.positive_cycle_capacity||c.errors.empty()||c.delta< -1e-9)return;
        c.delta=std::max(0.,c.delta);for(const auto& old:heap)if(old.errors==c.errors)return;
        auto cmp=[](const Cycle& a,const Cycle& b){return a.delta<b.delta;};
        if(heap.size()==config.positive_cycle_capacity){if(c.delta>=heap.front().delta)return;std::pop_heap(heap.begin(),heap.end(),cmp);heap.pop_back();}
        heap.push_back(std::move(c));std::push_heap(heap.begin(),heap.end(),cmp);
    }
    void synthesize(Candidate& base,const Candidate& other,std::vector<Cycle>& positive) {
        auto cs=difference(base,other);Ids flip;std::vector<size_t> logical,neg;
        for(size_t i=0;i<cs.size();++i){if(cs[i].obs){logical.push_back(i);if(cs[i].delta< -1e-9)neg.push_back(i);}else if(cs[i].delta< -1e-9){flip=sx(flip,cs[i].errors);++stats.improving_cycles;}}
        // Two disjoint logical operators form an observable-neutral cycle.
        // Find the least-cost even subset of the logical components.
        if(neg.size()&1){size_t remove=neg[0];for(size_t i:neg)if(cs[i].delta>cs[remove].delta)remove=i;size_t add=SIZE_MAX;for(size_t i:logical)if(cs[i].delta>= -1e-9&&(add==SIZE_MAX||cs[i].delta<cs[add].delta))add=i;if(add!=SIZE_MAX&&cs[add].delta< -cs[remove].delta)neg.push_back(add);else neg.erase(std::find(neg.begin(),neg.end(),remove));}
        double ld=0;for(size_t i:neg)ld+=cs[i].delta;if(ld< -1e-9){for(size_t i:neg)flip=sx(flip,cs[i].errors);stats.improving_cycles+=neg.size()/2;}
        if(!flip.empty()) {base.errors=sx(base.errors,flip);base.weight=0;for(uint32_t h:base.errors)base.weight+=hyper[h].w;positive.clear();}
        // Recompute retained cycles relative to the new representative.
        for(auto c:cs)if(!c.obs){c.delta=relative(c.errors,base);save_cycle(positive,std::move(c));}
        for(size_t a=0;a<logical.size();++a)for(size_t b=a+1;b<logical.size();++b){Cycle c{sx(cs[logical[a]].errors,cs[logical[b]].errors),0,false};c.delta=relative(c.errors,base);save_cycle(positive,std::move(c));}
    }
    double log_mass(const Candidate& base,const std::vector<Cycle>& cycles) const {
        Dsu ds(cycles.size());std::unordered_map<uint32_t,uint32_t> owner;
        for(size_t i=0;i<cycles.size();++i)for(uint32_t h:cycles[i].errors){auto[it,ok]=owner.emplace(h,i);if(!ok)ds.join(i,it->second);}
        std::map<uint32_t,Ids> groups;for(size_t i=0;i<cycles.size();++i)groups[ds.root(i)].push_back(i);
        double mass=-base.weight;
        for(auto& [r,g]:groups){std::set<Ids> span{{}};
            if(g.size()<10){for(uint32_t i:g){std::vector<Ids> next;next.reserve(span.size());for(const auto& v:span)next.push_back(sx(v,cycles[i].errors));for(auto& v:next)span.insert(std::move(v));}}
            else {for(uint32_t i:g)span.insert(cycles[i].errors);for(size_t i=0;i<g.size();++i)for(size_t j=i+1;j<g.size();++j)span.insert(sx(cycles[g[i]].errors,cycles[g[j]].errors));}
            double logsum=-std::numeric_limits<double>::infinity();for(const auto& v:span)logsum=logadd(logsum,-relative(v,base));mass+=logsum;
        }return mass;
    }
    std::pair<SyndromeType,bool> syndrome(const Ids& errors) const {SyndromeType syn(nd);bool obs=false;for(uint32_t h:errors){if(h>=hyper.size())throw std::out_of_range("LibraVariant hyperedge id");for(uint32_t d:hyper[h].dets)syn[d]^=1;obs^=hyper[h].obs;}return {std::move(syn),obs};}
    void validate(const Candidate& c,SyndromeRef given) const {auto [syn,o]=syndrome(c.errors);if(o!=c.obs)throw std::runtime_error("LibraVariant observable invariant failed");for(size_t d=0;d<nd;++d)if(syn[d]!=given[d])throw std::runtime_error("LibraVariant syndrome invariant failed");}
    DetailedResult decode(SyndromeRef syn) {
        ++stats.shots;std::vector<uint64_t> events;
        for(size_t k=0;k<(nd+63)/64;++k){uint64_t bits=syn.u64[k];while(bits){size_t d=k*64+std::countr_zero(bits);if(d<nd)events.push_back(d);bits&=bits-1;}}
        DetailedResult out;if(events.empty())return out;
        auto time=Clock::now();Candidate reps[2];reps[0]=candidate(*complement,events,true);events.push_back(nd);reps[1]=candidate(*complement,events,true);events.pop_back();
        if(reps[0].obs||!reps[1].obs)throw std::runtime_error("LibraVariant complementary class invariant failed");
        out.gap_db=std::abs(reps[0].weight-reps[1].weight)*10/std::log(10.);
        out.initial_prediction=reps[1].weight<reps[0].weight;
        if(config.profile)stats.complementary_ns+=ns(time);
        // Section V: when the gap is large, return the original unconstrained
        // correlated matcher. Constrained correlated matching is approximate,
        // so its cheaper representative need not predict the same class.
        if(config.large_gap_fallback&&out.gap_db>=config.gap_threshold_db){
            time=Clock::now();auto fallback=candidate(*original,events,false);
            out.fallback_used=true;out.prediction=fallback.obs;++stats.fallback_shots;
            reps[fallback.obs]=std::move(fallback);
            if(config.profile)stats.fallback_ns+=ns(time);
        }
        std::vector<Cycle> positive[2];
        if(out.gap_db<config.gap_threshold_db&&config.ensemble_size){
            out.ensemble_used=true;++stats.ensemble_shots;time=Clock::now();std::vector<Candidate> ensemble;ensemble.reserve(members.size()+2);ensemble.push_back(reps[0]);ensemble.push_back(reps[1]);
            for(size_t i=0;i<members.size();++i){prepare(i);auto c=candidate(*members[i].graph,events,false,&members[i].weights);
                // Consecutive identical assignments are redundant. Preserve
                // repetitions after intervening candidates: the base can change.
                if(ensemble.back().errors!=c.errors)ensemble.push_back(std::move(c));
            }
            if(config.profile)stats.ensemble_ns+=ns(time);
            time=Clock::now();
            for(size_t pass=0;pass<(config.degeneracy?2:1);++pass)for(const auto& c:ensemble)for(size_t cls=0;cls<2;++cls)synthesize(reps[cls],c,positive[cls]);
            if(config.profile)stats.synthesis_ns+=ns(time);
        }
        time=Clock::now();for(size_t i=0;i<2;++i){out.weights[i]=reps[i].weight;out.log_masses[i]=config.degeneracy?log_mass(reps[i],positive[i]):-reps[i].weight;if(config.validate)validate(reps[i],syn);}
        if(config.profile)stats.degeneracy_ns+=ns(time);
        if(!out.fallback_used)out.prediction=out.log_masses[1]>out.log_masses[0];
        out.errors=std::move(reps[out.prediction].errors);return out;
    }
};

LibraVariant::LibraVariant(const stim::DetectorErrorModel& dem):LibraVariant(dem,Config{}){}
LibraVariant::LibraVariant(const stim::DetectorErrorModel& dem,Config config):num_detectors(dem.count_detectors()),num_observables(dem.count_observables()),impl_(std::make_unique<Impl>(dem,config)){}
LibraVariant::~LibraVariant()=default;
LibraVariant::LibraVariant(LibraVariant&&) noexcept=default;
LibraVariant::DetailedResult LibraVariant::decode_detailed(SyndromeRef syn){return impl_->decode(syn);}
result_type LibraVariant::decode(SyndromeRef syn,ObsRef){auto r=decode_detailed(syn);result_type out;out.flipped_obs=ObsType(1);out.flipped_obs[0]=r.prediction;out.matching_data.gap=r.gap_db;out.matching_data.probability=std::exp(r.log_masses[r.prediction]);return out;}
void LibraVariant::prepare_ensemble(){for(size_t i=0;i<impl_->members.size();++i)impl_->prepare(i);}
LibraVariant::Dimensions LibraVariant::dimensions() const {return {num_detectors,impl_->hyper.size(),impl_->gedges.size(),impl_->rules};}
LibraVariant::Statistics LibraVariant::statistics() const{return impl_->stats;}
std::pair<SyndromeType,bool> LibraVariant::syndrome_of(const Ids& errors)const{return impl_->syndrome(errors);}
void LibraVariant::print_stats(std::ostream& out)const {auto s=statistics();out<<"LibraVariant shots="<<s.shots<<" ensemble_shots="<<s.ensemble_shots<<" improving_cycles="<<s.improving_cycles<<" complementary_ns="<<s.complementary_ns<<" ensemble_ns="<<s.ensemble_ns<<" synthesis_ns="<<s.synthesis_ns<<" degeneracy_ns="<<s.degeneracy_ns<<'\n';}
}
