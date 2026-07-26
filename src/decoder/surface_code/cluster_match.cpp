/*
 *  author: Suhas Vitta
 *  date:   20 May 2026
 * */

#if defined(VERILATOR_CLUSTER_MATCH)
#include <verilated.h>
#include "Vastrea.h"
#include "Vinitialize_neighbors.h"
#include "Vfilter.h"
#endif

#include "decoder/surface_code.h"

#if defined(VERILATOR_CLUSTER_MATCH)
double sc_time_stamp() { return 0; }
#endif

#include <PerfectMatching.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <limits>
#include <numeric>
#include <queue>
#include <type_traits>

namespace decoder
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

using DetIdType = ClusterMatch::det_id_type;
using AdjEntryType = ClusterMatch::adj_entry_type;
using AdjListType = ClusterMatch::adj_list_type;
using ClusterType = ClusterMatch::cluster_type;
using MwpmEdgeType = ClusterMatch::mwpm_edge_type;
using MatchingProblemType = ClusterMatch::matching_problem_type;
using QuantizationLevel = ClusterMatch::quantization_level;
using AssignmentType = MatchingData::assignment_type;

constexpr DetIdType BOUNDARY_ID{-1};

/*
 * Quantizes a weight according to the given `QuantizationLevel`
 * */
uint64_t _quantize(double, QuantizationLevel);

/*
 * Safe adjacency list update that also handles the case where the given
 * detector is already present in the adjacency list.
 * */
void _update_adjacency_list(AdjListType&, DetIdType, double p, ObsRef);

/*
 * Returns maximum of `1` and `ceil((d-1)/4)`
 * */
constexpr size_t _max_growth_steps(size_t code_distance);

/*
 * Returns number of MWPM edges that will be created for the given Hamming weight
 * */
constexpr size_t _get_mwpm_edge_count(size_t hw);

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*
 * Union find data structures used in `ClusterMatch::uf_compute_clusters()`
 * */

struct uf_type
{
public:
    uf_type* parent{nullptr};

    /*
     * `all_detectors` contains all detectors in the cluster (flipped or
     * otherwise). `flipped_detectors` only contains detectors in the
     * syndrome.
     * */
    std::vector<DetIdType> all_detectors;
    std::vector<DetIdType> flipped_detectors;

    /*
     * Boundary can have multiplicity, so we need to track it separately.
     * */
    bool has_boundary{false};
public:
    uf_type(DetIdType);

    uf_type* find();
    void merge(uf_type*);

    /*
     * Useful one-liners:
     * */
    size_t size() const { return all_detectors.size(); }
    size_t active_size() const { return flipped_detectors.size(); }
private:
    /*
     * `merge_find()` is a special find used in `merge()` that sets
     * `find()->parent = new_parent` and updates all `uf_type*` in
     * the `find()` traverse stack. This returns the old parent,
     * */
    uf_type* merge_find(uf_type* new_parent);
};

struct uf_growth_type
{
    DetIdType d;
    size_t      step{0};
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*
 * Data structure when computing distance between two detectors.
 * */

struct distance_type
{
    double w =std::numeric_limits<uint64_t>::max();
    ObsType frame_flips;
};

struct distance_queue_entry
{
    DetIdType d;
    double w;
};

struct distance_cmp
{
    bool operator()(const distance_queue_entry& a, const distance_queue_entry& b) const { return a.w > b.w; }
};

using DistanceQueueType = std::priority_queue<distance_queue_entry, 
                                                std::vector<distance_queue_entry>, 
                                                distance_cmp>;

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // anon

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

ClusterMatch::ClusterMatch(const stim::DetectorErrorModel& dem,
                                size_t _code_distance,
                                size_t _astrea_hw_max,
                                QuantizationLevel ql,
                                uint8_t _hw_emu_enable)
    :num_detectors(dem.count_detectors()),
    num_observables(dem.count_observables()),
    code_distance(_code_distance),
    astrea_hw_max(_astrea_hw_max),
    astrea_weight_quantization(ql),
    hw_emu_enable(_hw_emu_enable),
    adj_matrix_(num_detectors)
{
    // setup `adj_matrix_`
    dem.iter_flatten_error_instructions(
            [this] (const auto& inst)
            {
                const double pr = inst.arg_data[0];
                inst.for_separated_targets(
                        [this, pr] (const auto& grp)
                        {
                            std::vector<DetIdType> dets;
                            ObsType frame_flips(num_observables);
                            for (const auto& t : grp)
                            {
                                if (t.is_relative_detector_id())
                                    dets.push_back(static_cast<DetIdType>(t.val()));
                                else if (t.is_observable_id())
                                    frame_flips[t.val()] ^= 1;
                            }
                            
                            DetIdType d1 = dets[0],
                                        d2 = (dets.size() == 1) ? BOUNDARY_ID : dets[1];
                            _update_adjacency_list(adj_matrix_[d1], d2, pr, frame_flips);
                            if (d2 == BOUNDARY_ID)
                                _update_adjacency_list(boundary_adjacency_, d1, pr, frame_flips);
                            else
                                _update_adjacency_list(adj_matrix_[d2], d1, pr, frame_flips);
                        });
            });
}

const AdjListType&
ClusterMatch::adj_matrix(DetIdType d) const
{
    return (d == BOUNDARY_ID) ? boundary_adjacency_ : adj_matrix_[d];
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

result_type
ClusterMatch::decode(SyndromeRef syndrome)
{
    result_type out{.flipped_obs=ObsType(num_observables)};
    s_hamming_weight.add(syndrome.popcnt());

#if defined(VERILATOR_CLUSTER_MATCH)
    // 0. Initialize verilator context + modules:
    VerilatedContext v_ctx;
    Vastrea v_astrea(&v_ctx);
#else
    if (hw_emu_enable)
        std::cerr << "CLUSTER_MATCH::decode: `hw_emu_enable > 0` but VERILATOR_CLUSTER_MATCH macro is undefined." << _die{};
#endif

    // 1. Filter syndrome and remove isolated weight-1 errors.
    SyndromeType filtered_syndrome(syndrome);
    auto filter_out = filter_isolated_errors(filtered_syndrome);
    out.flipped_obs ^= filter_out.flipped_obs;
    out.matching_data.merge(filter_out.matching_data);

    s_filtered.add(syndrome.popcnt() - filtered_syndrome.popcnt());
    s_post_filter_hamming_weight.add(filtered_syndrome.popcnt());

    if (filtered_syndrome.popcnt() == 0)
        return out;

    // 2. Use UF algorithm to compute matching clusters
    auto clusters = uf_compute_clusters(filtered_syndrome);
    s_clusters.add(clusters.size());

    for (size_t i = 0; i < clusters.size(); i++)
    {
        auto& cl = clusters[i];

        s_cluster_size.add(cl.all.size());
        s_cluster_hamming_weight.add(cl.flipped.size());

        MatchingProblemType mp;
        result_type mp_result;

        // 3. Compute pairwise distances for all detection events in the cluster.
        mp = synthesize_matching_problem(std::move(cl));

        // 4. Run Astrea to get correction for cluster.
#if defined(VERILATOR_CLUSTER_MATCH)
        if (hw_emu_enable & hw_emu_flag::astrea)
            mp_result = v_solve_matching_problem(std::move(mp), v_astrea);
        else
            mp_result = solve_matching_problem(std::move(mp), i);
#else
        mp_result = solve_matching_problem(std::move(mp), i);
#endif
        out.flipped_obs ^= mp_result.flipped_obs;
        out.matching_data.merge(mp_result.matching_data);
    }

#if defined(VERILATOR_CLUSTER_MATCH)
    v_astrea.final();
#endif

    return out;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void
ClusterMatch::print_stats(std::ostream& ostrm) const
{
    // report both statistics and any useful RTL data:
    std::vector<size_t> degree_array(num_detectors);
    std::transform(adj_matrix_.begin(), adj_matrix_.end(), degree_array.begin(),
                    [] (const auto& al) { return al.size(); });

    const size_t max_steps = _max_growth_steps(code_distance);
    const size_t max_degree = *std::max_element(degree_array.begin(), degree_array.end());
    const size_t boundary_degree = boundary_adjacency_.size();

    ostrm << "CLUSTER_MATCH-----------------------------------\n";
    print_stat(ostrm, "DISTANCE", code_distance);
    print_stat(ostrm, "MAX_STEPS", max_steps);
    print_stat(ostrm, "MAX_DEGREE", max_degree);
    print_stat(ostrm, "BOUNDARY_DEGREE", boundary_degree);

    ostrm << "\nCLUSTER_COUNT\t" << s_clusters.to_string_some()
            << "\nFILTER_COUNT\t" << s_filtered.to_string_some()
            << "\nHAMMING_WEIGHT\t" << s_hamming_weight.to_string_some()
            << "\nPOST_FILTER_HAMMING_WEIGHT\t" << s_post_filter_hamming_weight.to_string_some()
            << "\nCLUSTER_HAMMING_WEIGHT\t" << s_cluster_hamming_weight.to_string_some()
            << "\nCLUSTER_SIZE\t" << s_cluster_size.to_string_some()
            << "\nGROWTH_TICKS\t" << s_growth_ticks.to_string_some()
            << "\nSYNTHESIS_TICKS\t" << s_synthesis_ticks.to_string_some()
            << "\nSYNTHESIS_TICKS_NORMALIZED\t" << s_synthesis_ticks_norm.to_string_some()
            << "\n";
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

result_type
ClusterMatch::filter_isolated_errors(SyndromeRef syndrome)
{
    // count active degree of all syndrome bits:
    const bool count_boundary = (syndrome.popcnt() & 1);
    std::vector<size_t> active_degree(num_detectors, 0);
    std::vector<std::pair<DetIdType, AdjEntryType>> active_companion(num_detectors,
                                                                        std::make_pair(0, AdjEntryType{.frame_flips=ObsType(1)}) );
    size_t boundary_degree{0};
    for (size_t i = 0; i < num_detectors; i++)
    {
        if (syndrome[i])
        {
            for (const auto& e : adj_matrix(i))
            {
                if (e.d == BOUNDARY_ID)
                {
                    if (count_boundary)
                    {
                        active_degree[i]++;
                        boundary_degree++;
                        active_companion[i] = std::make_pair(BOUNDARY_ID, e);
                    }
                }
                else if (syndrome[e.d])
                {
                    active_degree[i]++;
                    active_companion[i] = std::make_pair(e.d, e);
                }
            }
        }
    }

    // process the syndrome again:
    result_type out{};
    for (size_t i = 0; i < num_detectors; i++)
    {
        if (syndrome[i] && active_degree[i] == 1)
        {
            auto [j, e] = active_companion[i];
            if ((j == BOUNDARY_ID && boundary_degree == 1) || (j != BOUNDARY_ID && active_degree[j] == 1))
            {
                syndrome[i] ^= 1;
                if (j != BOUNDARY_ID)
                    syndrome[j] ^= 1;
                out.flipped_obs ^= e.frame_flips;
                // update matching data:
                auto w_qu = _quantize(-std::log(e.pr), astrea_weight_quantization);
                AssignmentType a{.d1=i, .d2=j, .pr=e.pr, .w_qu=w_qu, .frame_flips=e.frame_flips};
                out.matching_data.add(a);
            }
        }
    }
    return out;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

std::vector<ClusterType>
ClusterMatch::uf_compute_clusters(SyndromeRef syndrome)
{
    // initialize `uf_pool` (storage for UF data structures) and
    // `growth_fifo` (what detectors to traverse from) using `syndrome`
    std::vector<uf_type> uf_pool;
    std::vector<uf_type*> uf_lookup(num_detectors, nullptr);
    std::deque<uf_growth_type> growth_fifo;

    uf_pool.reserve(syndrome.popcnt());
    for (DetIdType i = 0; i < num_detectors; i++)
    {
        if (syndrome[i])
        {
            uf_pool.push_back(uf_type(i));
            uf_lookup[i] = &uf_pool.back();
            growth_fifo.push_back(uf_growth_type{.d=i});
        }
    }

    const size_t max_steps = _max_growth_steps(code_distance);
    int unique_clusters{uf_pool.size()};
    uint64_t tick{0};
    while (growth_fifo.size() > 0)
    {
        auto g = std::move(growth_fifo.front());
        growth_fifo.pop_front();

        // do not traverse too much: we can terminate once the cluster
        // has grown to an amount proportional to the code distance.
        if (g.step >= max_steps)
            continue;
        tick++;
        const DetIdType d1 = g.d;
        assert(d1 != BOUNDARY_ID);

        // run find on `uf_lookup[d1].owner` now so we have the updated owner:
        uf_lookup[d1] = uf_lookup[d1]->find();
        auto* uf1 = uf_lookup[d1];
        // no need to traverse if `uf1->active_size()` is large enough
        if (uf1->active_size() >= astrea_hw_max)
            continue;
        
        // traverse:
        for (const auto& e : adj_matrix(g.d))
        {
            const DetIdType d2 = e.d;

            // handle `d2 == BOUNDARY_ID` specially:
            if (d2 == BOUNDARY_ID)
            {
                uf1->has_boundary = true;
                continue;
            }

            // run find on `uf_lookup[d2].owner`:
            if (uf_lookup[d2] != nullptr)
                uf_lookup[d2] = uf_lookup[d2]->find();
            auto* uf2 = uf_lookup[d2];

            // compare `uf1` and `uf2` 
            if (uf1 == uf2)
                continue;

            if (uf2 == nullptr)
            {
                uf1->all_detectors.push_back(d2);
                uf_lookup[d2] = uf1;
                growth_fifo.push_back(uf_growth_type{.d=d2, .step=g.step+1});
            }
            // only merge `uf1` and `uf2` if they are beneath the HW threshold
            else
            {
                if (uf1->active_size() + uf2->active_size() <= astrea_hw_max)
                {
                    uf1->merge(uf2);
                    uf_lookup[d2] = uf1;
                    unique_clusters--;
                }
            }
        }
    }

    s_growth_ticks.add(tick);

    // form clusters:
    std::vector<ClusterType> clusters;
    clusters.reserve(unique_clusters);
    for (const auto& uf : uf_pool)
    {
        if (uf.parent != nullptr)
            continue;
        ClusterType cl{ std::move(uf.all_detectors), 
                         std::move(uf.flipped_detectors) };
        // only add boundary if `cl` is already odd:
        if (cl.flipped.size() % 2 == 1)
        {
            cl.flipped.push_back(BOUNDARY_ID);
            cl.all.push_back(BOUNDARY_ID);
        }
        else if (uf.has_boundary)
        {
            cl.all.push_back(BOUNDARY_ID);
        }
        clusters.push_back(cl);
    }
    return clusters;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

MatchingProblemType
ClusterMatch::synthesize_matching_problem(ClusterType cl)
{
    assert((cl.flipped.size() % 2) == 0);
    const size_t n = cl.all.size();
    const size_t hw = cl.flipped.size();
    
    // create a map that maps detector id to index in `cl.all`
    std::unordered_map<DetIdType, size_t> idx_map;
    idx_map.reserve(n);
    for (size_t i = 0; i < n; i++)
        idx_map.insert({cl.all[i], i});

    // run `dijkstra's` `n-1` times to create `MatchingProblemType`
    MatchingProblemType mp;
    mp.edges.reserve(hw*(hw-1)/2);

    const distance_type fill_val{.frame_flips=ObsType(num_observables)};
    std::vector<distance_type> dist(n, fill_val);
    uint64_t tick{0};
    for (size_t ii = 0; ii < hw-1; ii++)
    {
        std::fill(dist.begin(), dist.end(), fill_val);
        const DetIdType d1 = cl.flipped[ii];
        dist[idx_map[d1]].w = 0.0;
        DistanceQueueType pq;
        pq.push({d1, 0.0});
        while (pq.size() > 0)
        {
            auto e = std::move(pq.top());
            pq.pop();
            const auto z1 = e.d;
            const auto w1 = e.w;
            const auto i = idx_map.at(z1);
            if (std::abs(w1 - dist[i].w) > 1e-6)
                continue;
            tick++;
            // Accumulate the raw `-log(pr)`; quantization happens once below.
            for (const auto& x : adj_matrix(z1))
            {
                const auto z2 = x.d;
                auto idx_it = idx_map.find(z2);
                if (idx_it == idx_map.end())
                    continue;
                const auto j = idx_it->second;
                const double w2 = w1 + (-std::log(x.pr));
                if (w2 < dist[j].w)
                {
                    dist[j].w = w2;
                    dist[j].frame_flips = dist[i].frame_flips ^ x.frame_flips;
                    pq.push({z2, w2});
                }
            }
        }

        // create mwpm edges: `dist[j].w` is the accumulated `-log` path weight, so the
        // matching's error probability is `exp(-w)` and its quantized weight `_quantize(w)`.
        for (size_t jj = ii+1; jj < hw; jj++)
        {
            const DetIdType d2 = cl.flipped[jj];
            const size_t j = idx_map.at(d2);
            MwpmEdgeType e{ .d1=d1,
                                .d2=d2,
                                .pr=std::exp(-dist[j].w),
                                .w_qu=_quantize(dist[j].w, astrea_weight_quantization),
                                .frame_flips=dist[j].frame_flips };
            mp.edges.push_back(e);
        }
    }
    double tick_norm = static_cast<double>(tick) / static_cast<double>(hw);
    s_synthesis_ticks.add(tick);
    s_synthesis_ticks_norm.add(tick_norm);
    mp.detectors = std::move(cl.flipped);
    return mp;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

result_type
ClusterMatch::solve_matching_problem(MatchingProblemType mp, int cluster_id)
{
    // create index map for `mp.detectors`
    std::unordered_map<DetIdType, size_t> idx_map;
    idx_map.reserve(mp.detectors.size());
    for (size_t i = 0; i < mp.detectors.size(); i++)
        idx_map[mp.detectors[i]] = i;

    [[ maybe_unused ]]  const size_t n = mp.detectors.size(),
                                     m = mp.edges.size();
    if (n > astrea_hw_max)
    {
        std::cerr << "CLUSTER_MATCHING::solve_matching_problem: got matching problem with HW = " 
                    << n << " > HW_MAX (" << astrea_hw_max << ")" << _die{};
    }
    
    b5::PerfectMatching pm(n, m); 
    pm.options.verbose = false;
    for (size_t k = 0; k < m; k++)
    {
        const auto& e = mp.edges[k];
        const size_t i = idx_map[e.d1],
                     j = idx_map[e.d2];
        [[ maybe_unused ]] auto _k = pm.AddEdge(i, j, e.w_qu);
        assert(k == _k);
    }

    pm.Solve();

    // Retrieve the solution to the MWPM problem:
    result_type out{.flipped_obs=ObsType(num_observables)};
    for (size_t i = 0; i < m; i++)
    {
        if (pm.GetSolution(i))
        {
            const auto& e = mp.edges[i];
            out.flipped_obs ^= e.frame_flips;
            // update matching data:
            AssignmentType a{e};
            a.matching_step = 1;
            a.cluster_id = cluster_id;
            out.matching_data.add(a);
        }
    }
    
    return out;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

uint64_t
_quantize(double w, QuantizationLevel ql)
{
    uint64_t q;
    if (ql == QuantizationLevel::b4)
    {
        q = std::round(w); 
        q = std::min(q, uint64_t{15});
        return q;
    }
    else if (ql == QuantizationLevel::b8)
    {
        q = std::round(15*w);
        q = std::min(q, uint64_t{255});
    }
    else if (ql == QuantizationLevel::b16)
    {
        q = std::round(100*w);
        q = std::min(q, uint64_t{(1ull<<16)-1});
    }
    else
    {
        q = std::round(1'000'000*w);
        q = std::min(q, uint64_t{(1ull<<32)-1});
    }
    return q;
}

void
_update_adjacency_list(AdjListType& adj, DetIdType d, double p, ObsRef frame_flips)
{
    // check if `d` is in `adj`
    auto adj_it = std::find_if(adj.begin(), adj.end(),
                        [d] (const auto& e) { return e.d == d; });
    if (adj_it != adj.end())
    {
        // update weight;
        adj_it->pr = (1-adj_it->pr)*p + (1-p)*adj_it->pr;
    }
    else
    {
        AdjEntryType e{d, p, ObsType{frame_flips}};
        adj.push_back(e);
    }
}

constexpr size_t
_max_growth_steps(size_t d)
{
    double g = static_cast<double>(d-1) / 4.0 + 1.0;
    g = std::max(1.0, std::ceil(g));
    return static_cast<size_t>( std::round(g) );
}

constexpr size_t
_get_mwpm_edge_count(size_t hw)
{
    if (hw & 1)
        return _get_mwpm_edge_count(hw+1);
    else if (hw == 2)
        return 1;
    else
        return (hw-1) * _get_mwpm_edge_count(hw-2);
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

uf_type::uf_type(DetIdType d)
    :all_detectors{d},
    flipped_detectors{d}
{}

uf_type*
uf_type::find()
{
    std::vector<uf_type*> visited;
    visited.reserve(2);

    uf_type* p{this};
    while (p->parent != nullptr)
    {
        visited.push_back(p);
        p = p->parent;
    }
    for (auto* uf : visited)
        uf->parent = p;
    return p;
}

void
uf_type::merge(uf_type* other)
{
    // update `other`
    auto* p = find();
    uf_type* old_parent = other->merge_find(p);
    for (DetIdType d : old_parent->all_detectors)
        p->all_detectors.push_back(d);
    for (DetIdType d : old_parent->flipped_detectors)
        p->flipped_detectors.push_back(d);
    p->has_boundary |= old_parent->has_boundary;
}

uf_type*
uf_type::merge_find(uf_type* np)
{
    std::vector<uf_type*> visited;
    visited.reserve(2);

    uf_type* x{this};
    while (x->parent != nullptr)
    {
        visited.push_back(x);
        x = x->parent;
    }

    for (auto* uf : visited)
        uf->parent = np;
    x->parent = np;
    return x;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // anon

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace decoder

#if defined(VERILATOR_CLUSTER_MATCH)
#include "cluster_match.v.cpp"
#endif
