/*
 *  author: Suhas Vitta
 *  date:   20 May 2026
 * */

#include "decoder/surface_code.h"

#include <PerfectMatching.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
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
int64_t _quantize(double, quantization_level);

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
 * Data structure used in Floyd-Warshall. We need to track both
 * frame flips and the distance (`w`).
 * */

struct fw_data_type
{
    int64_t w;
    obs_type frame_flips;
};

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
                            if (d2 != BOUNDARY_ID)
                                _update_adjacency_list(adj_matrix_[d2], d1, pr, frame_flips);
                        });
            });
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

result_type
CLUSTER_MATCH::decode(syndrome_ref syndrome, LOGGER& logger)
{
    // 1. compute clusters:
    logger.info(1) << "uf_compute_clusters ------------------------------\n";
    logger.tab_level++;
    auto clusters = uf_compute_clusters(syndrome, logger);
    logger.tab_level--;

    // 2. synthesize and decode matching problems (one per cluster)
    result_type out{.flipped_obs=obs_type(num_observables)};
    logger.info(1) << "performing matching on clusters (count = " << clusters.size() << ") ---------\n";
    for (size_t i = 0; i < clusters.size(); i++)
    {
        auto& cl = clusters[i];

        logger.info(1) << "cluster " << i 
                        << ", size = " << cl.all.size() 
                        << ", hw = " << cl.flipped.size() 
                        << ", detectors =";
        for (det_id_type d : cl.flipped)
            logger.info(1) << " " << d;
        logger.info(1) << "\n";
        logger.tab_level++;

        logger.info(1) << "synthesize_matching_problem:" << "\n";
        logger.tab_level++;
        auto mp = synthesize_matching_problem(std::move(cl), logger);
        logger.tab_level--;

        logger.info(1) << "solve_matching_problem:" << "\n";
        logger.tab_level++;
        auto mp_result = solve_matching_problem(std::move(mp), logger);
        logger.tab_level--;

        // merge `mp_result` with `out`
        out.flipped_obs ^= mp_result.flipped_obs;
    }
    logger.tab_level--;

    return out;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

std::vector<cluster_type>
CLUSTER_MATCH::uf_compute_clusters(syndrome_ref syndrome, LOGGER& logger)
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
    size_t unique_clusters{uf_pool.size()};
    while (growth_fifo.size() > 0 && unique_clusters > 1)
    {
        auto g = std::move(growth_fifo.front());
        growth_fifo.pop_front();

        // do not traverse too much: we can terminate once the cluster
        // has grown to an amount proportional to the code distance.
        if (g.step >= max_steps)
            continue;

        const det_id_type d1 = g.d;
        assert(d1 != BOUNDARY_ID);
        
        // run find on `uf_lookup[d1].owner` now so we have the updated
        // owner:
        uf_lookup[d1] = uf_lookup[d1]->find();
        auto* uf1 = uf_lookup[d1];

        // no need to traverse if `uf1->active_size()` is large enough
        if (uf1->active_size() >= astrea_hw_max)
            continue;
        
        // traverse:
        for (const auto& e : adj_matrix_[g.d])
        {
            const det_id_type d2 = e.d;

            // handle `d2 == BOUNDARY_ID` specially:
            if (d2 == BOUNDARY_ID)
            {
                uf1->all_detectors.push_back(d2);
                uf1->flipped_detectors.push_back(d2);
                // do not traverse for `BOUNDARY_ID`
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
                // add `d2` to `uf` and push `d2` onto the `growth_fifo`
                uf1->all_detectors.push_back(d2);
                uf_lookup[d2] = uf1;
                growth_fifo.push_back(uf_growth_type{.d=d2, .step=g.step+1});
            }
            // only merge `uf1` and `uf2` if they are beneath the HW threshold
            else if (uf1->active_size() + uf2->active_size() <= astrea_hw_max)
            {
                uf1->merge(uf2);
                uf_lookup[d2] = uf1;
                unique_clusters--;
            }
        }
    }

    // form clusters:
    std::vector<cluster_type> clusters;
    clusters.reserve(unique_clusters);
    for (const auto& uf : uf_pool)
    {
        if (uf.parent != nullptr)
            continue;
        cluster_type cl{ std::move(uf.all_detectors), 
                         std::move(uf.flipped_detectors) };
        clusters.push_back(cl);
    }
    return clusters;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

matching_problem_type
CLUSTER_MATCH::synthesize_matching_problem(cluster_type&& cl, LOGGER& logger)
{
    // first, clean up `flipped` if it is odd.
    if (cl.flipped.size() & 1)
    {
        // Check if `BOUNDARY_ID` is already in `cl.flipped`
        // If it already exists, remove it. Otherwise add it.
        auto d_it = std::find(cl.flipped.begin(), cl.flipped.end(), BOUNDARY_ID);
        if (d_it == cl.flipped.end())
            cl.flipped.push_back(BOUNDARY_ID);
        else
            cl.flipped.erase(d_it);
    }
    assert((cl.flipped.size() & 1) == 0);

    // we need to do some setup before we run floyd-warshall
    const size_t n = cl.all.size();
    const fw_data_type fw_fill_val{.frame_flips=obs_type(num_observables)};
    std::vector<fw_data_type> dist(n*n, fw_fill_val);
    
    // create a map that maps detector id to index in `cl.all`
    std::unordered_map<det_id_type, size_t> idx_map;
    idx_map.reserve(n);
    for (size_t i = 0; i < n; i++)
        idx_map.insert({cl.all[i], i});

    // initialize `dist` using the edges in `adj_matrix_`
    for (size_t i = 0; i < n; i++)
    {
        const det_id_type d1 = cl.all[i];
        if (d1 == BOUNDARY_ID)
            continue;
        for (const auto& e : adj_matrix_[d1])
        {
            const det_id_type d2{e.d};
            auto idx_it = idx_map.find(d2);
            if (idx_it == idx_map.end())
                continue;
            const size_t j = idx_it->second;
            // if the weight is already nonzero, then this entry has already
            // been allocated:
            if (dist[i*n+j].w > 0.0)
                continue;
            dist[i*n+j].w = _quantize(-std::log(e.pr), astrea_weight_quantization);
            dist[i*n+j].frame_flips ^= e.frame_flips;
            // copy the data over to `j,i`
            dist[j*n+i] = dist[i*n+j];
        }
    }

    // run `floyd_warshall`:
    for (size_t k = 0; k < n; k++)
    {
        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                if (dist[i*n+j].w > dist[i*n+k].w + dist[j*n+k].w)
                {
                    dist[i*n+j].w = dist[i*n+k].w + dist[j*n+k].w;
                    dist[i*n+j].frame_flips = dist[i*n+k].frame_flips ^ dist[j*n+k].frame_flips;
                }
            }
        }
    }

    // create matching problem:
    matching_problem_type mp{.detectors=std::move(cl.flipped)};
    mp.edges.reserve(_get_mwpm_edge_count(mp.detectors.size()));
    for (size_t i = 0; i < mp.detectors.size(); i++)
    {
        const det_id_type d1 = mp.detectors[i];
        const size_t ii = idx_map.at(d1);
        for (size_t j = i+1; j < mp.detectors.size(); j++)
        {
            const det_id_type d2 = mp.detectors[j];
            const size_t jj = idx_map.at(d2);
            mwpm_edge_type e{ .d1=d1, 
                                .d2=d2, 
                                .w_qu=dist[ii*n+jj].w, 
                                .frame_flips=std::move(dist[ii*n+jj].frame_flips) };
            mp.edges.push_back(e);
        }
    }
    return mp;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

result_type
CLUSTER_MATCH::solve_matching_problem(matching_problem_type&& mp, LOGGER& logger)
{
    // create index map for `mp.detectors`
    std::unordered_map<det_id_type, size_t> idx_map;
    idx_map.reserve(mp.detectors.size());
    for (size_t i = 0; i < mp.detectors.size(); i++)
        idx_map[mp.detectors[i]] = i;

    b5::PerfectMatching pm(mp.detectors.size(), mp.edges.size()); 
    for (size_t k = 0; k < mp.edges.size(); k++)
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
    for (size_t i = 0; i < mp.edges.size(); i++)
        if (pm.GetSolution(i))
            out.flipped_obs ^= mp.edges[i].frame_flips;
    return out;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

int64_t
_quantize(double w, quantization_level ql)
{
    if (ql == quantization_level::b4)
    {
        int64_t q = std::round(w); 
        q = std::min(q, int64_t{15});
        return q;
    }
    else if (ql == quantization_level::b8)
    {
        int64_t q = std::round(15*w);
        q = std::min(q, int64_t{255});
    }
    else if (ql == quantization_level::b16)
    {
        int64_t q = std::round(1000*w);
        q = std::min(q, int64_t{(1ll<<16)-1});
    }
    else
    {
        int64_t q = std::round(1'000'000*w);
        q = std::min(q, int64_t{(1ll<<32)-1});
    }

    return -1;
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
    parent = find();
    uf_type* old_parent = other->merge_find(parent);
    for (det_id_type d : old_parent->all_detectors)
        all_detectors.push_back(d);
    for (det_id_type d : old_parent->flipped_detectors)
        flipped_detectors.push_back(d);
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
