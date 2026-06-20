/*
 *  author: Suhas Vitta
 *  date:   20 May 2026
 * */

#include "decoder/surface_code.h"

#include <PerfectMatching.h>

#include <algorithm>
#include <cassert>
#include <cmath>
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

using det_id_type = CLUSTER_MATCH::det_id_type;
using adj_entry_type = CLUSTER_MATCH::adj_entry_type;
using adj_list_type = CLUSTER_MATCH::adj_list_type;
using cluster_type = CLUSTER_MATCH::cluster_type;
using mwpm_edge_type = CLUSTER_MATCH::mwpm_edge_type;
using matching_problem_type = CLUSTER_MATCH::matching_problem_type;
using quantization_level = CLUSTER_MATCH::quantization_level;

constexpr det_id_type BOUNDARY_ID{-1};

/*
 * Quantizes a weight according to the given `quantization_level`
 * */
uint64_t _quantize(double, quantization_level);

/*
 * Safe adjacency list update that also handles the case where the given
 * detector is already present in the adjacency list.
 * */
void _update_adjacency_list(adj_list_type&, det_id_type, double p, obs_ref);

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
 * Union find data structures used in `CLUSTER_MATCH::uf_compute_clusters()`
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
    std::vector<det_id_type> all_detectors;
    std::vector<det_id_type> flipped_detectors;

    /*
     * Boundary can have multiplicity, so we need to track it separately.
     * */
    bool has_boundary{false};
public:
    uf_type(det_id_type);

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
    det_id_type d;
    size_t      step{0};
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*
 * Data structure when computing distance between two detectors.
 * */

struct distance_type
{
    uint64_t w =std::numeric_limits<uint64_t>::max();
    obs_type frame_flips;
};

struct distance_queue_entry
{
    det_id_type d;
    uint64_t w;
};

struct distance_cmp
{
    bool operator()(const distance_queue_entry& a, const distance_queue_entry& b) const { return a.w > b.w; }
};

using distance_queue_type = std::priority_queue<distance_queue_entry, 
                                                std::vector<distance_queue_entry>, 
                                                distance_cmp>;

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // anon

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

CLUSTER_MATCH::CLUSTER_MATCH(const stim::DetectorErrorModel& dem,
                                size_t _code_distance,
                                size_t _astrea_hw_max,
                                quantization_level ql)
    :num_detectors(dem.count_detectors()),
    num_observables(dem.count_observables()),
    code_distance(_code_distance),
    astrea_hw_max(_astrea_hw_max),
    astrea_weight_quantization(ql),
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
                            std::vector<det_id_type> dets;
                            obs_type frame_flips(num_observables);
                            for (const auto& t : grp)
                            {
                                if (t.is_relative_detector_id())
                                    dets.push_back(static_cast<det_id_type>(t.val()));
                                else if (t.is_observable_id())
                                    frame_flips[t.val()] ^= 1;
                            }
                            
                            det_id_type d1 = dets[0],
                                        d2 = (dets.size() == 1) ? BOUNDARY_ID : dets[1];
                            _update_adjacency_list(adj_matrix_[d1], d2, pr, frame_flips);
                            if (d2 == BOUNDARY_ID)
                                _update_adjacency_list(boundary_adjacency_, d1, pr, frame_flips);
                            else
                                _update_adjacency_list(adj_matrix_[d2], d1, pr, frame_flips);
                        });
            });
}

const adj_list_type&
CLUSTER_MATCH::adj_matrix(det_id_type d) const
{
    return (d == BOUNDARY_ID) ? boundary_adjacency_ : adj_matrix_[d];
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

result_type
CLUSTER_MATCH::decode(syndrome_ref syndrome)
{
    result_type out{.flipped_obs=obs_type(num_observables)};

    s_hamming_weight.add(syndrome.popcnt());

    syndrome_type filtered_syndrome(syndrome);
    auto filter_out = filter_isolated_errors(filtered_syndrome);
    out.flipped_obs ^= filter_out.flipped_obs;
    s_filtered.add(syndrome.popcnt() - filtered_syndrome.popcnt());
    s_post_filter_hamming_weight.add(filtered_syndrome.popcnt());

    if (filtered_syndrome.popcnt() == 0)
        return out;

    auto clusters = uf_compute_clusters(filtered_syndrome);
    s_clusters.add(clusters.size());

    for (size_t i = 0; i < clusters.size(); i++)
    {
        auto& cl = clusters[i];

        s_cluster_size.add(cl.all.size());
        s_cluster_hamming_weight.add(cl.flipped.size());

        auto mp = synthesize_matching_problem(std::move(cl));
        auto mp_result = solve_matching_problem(std::move(mp));
        out.flipped_obs ^= mp_result.flipped_obs;
    }

    return out;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void
CLUSTER_MATCH::print_stats(std::ostream& ostrm) const
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
CLUSTER_MATCH::filter_isolated_errors(syndrome_ref syndrome)
{
    // count active degree of all syndrome bits:
    std::vector<size_t> active_degree(num_detectors, 0);
    std::vector<std::pair<det_id_type, obs_type>> active_companion(num_detectors, 
                                                                    {0,obs_type(num_observables)});
    for (size_t i = 0; i < num_detectors; i++)
    {
        if (syndrome[i])
        {
            for (const auto& e : adj_matrix(i))
            {
                if (e.d == BOUNDARY_ID)
                    continue;
                if (syndrome[e.d])
                {
                    active_degree[i]++;
                    active_companion[i] = std::make_pair(e.d, e.frame_flips);
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
            auto [j, frame_flips] = active_companion[i];
            if (active_degree[j])
            {
                syndrome[i] ^= 1;
                syndrome[j] ^= 1;
                out.flipped_obs ^= frame_flips;
            }
        }
    }
    return out;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

std::vector<cluster_type>
CLUSTER_MATCH::uf_compute_clusters(syndrome_ref syndrome)
{
    // initialize `uf_pool` (storage for UF data structures) and
    // `growth_fifo` (what detectors to traverse from) using `syndrome`
    std::vector<uf_type> uf_pool;
    std::vector<uf_type*> uf_lookup(num_detectors, nullptr);
    std::deque<uf_growth_type> growth_fifo;

    uf_pool.reserve(syndrome.popcnt());
    for (det_id_type i = 0; i < num_detectors; i++)
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
        const det_id_type d1 = g.d;
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
            const det_id_type d2 = e.d;

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
    std::vector<cluster_type> clusters;
    clusters.reserve(unique_clusters);
    for (const auto& uf : uf_pool)
    {
        if (uf.parent != nullptr)
            continue;
        cluster_type cl{ std::move(uf.all_detectors), 
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

matching_problem_type
CLUSTER_MATCH::synthesize_matching_problem(cluster_type&& cl)
{
    assert((cl.flipped.size() % 2) == 0);
    const size_t n = cl.all.size();
    const size_t hw = cl.flipped.size();
    
    // create a map that maps detector id to index in `cl.all`
    std::unordered_map<det_id_type, size_t> idx_map;
    idx_map.reserve(n);
    for (size_t i = 0; i < n; i++)
        idx_map.insert({cl.all[i], i});

    // run `dijkstra's` `n-1` times to create `matching_problem_type`
    matching_problem_type mp;
    mp.edges.reserve(hw*(hw-1)/2);

    const distance_type fill_val{.frame_flips=obs_type(num_observables)};
    std::vector<distance_type> dist(n, fill_val);
    uint64_t tick{0};
    for (size_t ii = 0; ii < hw-1; ii++)
    {
        std::fill(dist.begin(), dist.end(), fill_val);
        const det_id_type d1 = cl.flipped[ii];
        dist[idx_map[d1]].w = 0;
        distance_queue_type pq;
        pq.push({d1, 0});
        while (pq.size() > 0)
        {
            auto e = std::move(pq.top());
            pq.pop();
            const auto z1 = e.d;
            const auto w1 = e.w;
            const auto i = idx_map.at(z1);
            if (w1 != dist[i].w)
                continue;
            tick++;
            for (const auto& x : adj_matrix(z1))
            {
                const auto z2 = x.d;
                auto idx_it = idx_map.find(z2);
                if (idx_it == idx_map.end())
                    continue;
                const auto j = idx_it->second;
                const uint64_t w_qu = _quantize(-std::log(x.pr), astrea_weight_quantization);
                const auto w2 = w1 + w_qu;
                if (w2 < dist[j].w)
                {
                    dist[j].w = w2;
                    dist[j].frame_flips = dist[i].frame_flips ^ x.frame_flips;
                    pq.push({z2, w2});
                }
            }
        }

        // create mwpm edges:
        for (size_t jj = ii+1; jj < hw; jj++)
        {
            const det_id_type d2 = cl.flipped[jj];
            const size_t j = idx_map.at(d2);
            mwpm_edge_type e{ .d1=d1,
                                .d2=d2,
                                .w_qu=dist[j].w,
                                .frame_flips=std::move(dist[j].frame_flips) };
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
CLUSTER_MATCH::solve_matching_problem(matching_problem_type&& mp)
{
    // create index map for `mp.detectors`
    std::unordered_map<det_id_type, size_t> idx_map;
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
    result_type out{.flipped_obs=obs_type(num_observables)};
    for (size_t i = 0; i < m; i++)
    {
        if (pm.GetSolution(i))
        {
            const auto& e = mp.edges[i];
            out.flipped_obs ^= e.frame_flips;
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
_quantize(double w, quantization_level ql)
{
    uint64_t q;
    if (ql == quantization_level::b4)
    {
        q = std::round(w); 
        q = std::min(q, uint64_t{15});
        return q;
    }
    else if (ql == quantization_level::b8)
    {
        q = std::round(15*w);
        q = std::min(q, uint64_t{255});
    }
    else if (ql == quantization_level::b16)
    {
        q = std::round(1000*w);
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
_update_adjacency_list(adj_list_type& adj, det_id_type d, double p, obs_ref frame_flips)
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
        adj_entry_type e{d, p, obs_type{frame_flips}};
        adj.push_back(e);
    }
}

constexpr size_t
_max_growth_steps(size_t d)
{
    double g = static_cast<double>(d-1) / 4.0;
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

uf_type::uf_type(det_id_type d)
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
    for (det_id_type d : old_parent->all_detectors)
        p->all_detectors.push_back(d);
    for (det_id_type d : old_parent->flipped_detectors)
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
