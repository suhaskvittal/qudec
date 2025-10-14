/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 */

#include <limits>
#include <queue>

namespace graph
{

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

template <class W, class GRAPH_TYPE, class WEIGHT_FUNCTION, class EARLY_TERM_ITER> DIJKSTRA_RESULT<W>
dijkstra(const GRAPH_TYPE& gr,
        GRAPH_COMPONENT_ID src,
        const WEIGHT_FUNCTION& wf,
        bool terminate_early,
        EARLY_TERM_ITER et_begin,
        EARLY_TERM_ITER et_end)
{
    constexpr int64_t UNDEFINED{-19243987};  // some random number -- unlikely collision

    std::vector<W> dist(gr.get_vertices().size(), std::numeric_limits<W>::max());
    std::vector<GRAPH_COMPONENT_ID> prev(gr.get_vertices().size(), UNDEFINED);

    // initialize priority queue:
    struct queue_entry
    {
        GRAPH_COMPONENT_ID id;
        W                  dist;
    };

    auto cmp = [] (const queue_entry& a, const queue_entry& b) { return a.dist > b.dist; };
    std::priority_queue<queue_entry, std::vector<queue_entry>, decltype(cmp)> pq(cmp);

    pq.push({src, 0});
    dist[src] = 0;
    prev[src] = src;

    // if `ENABLE_EARLY_TERM` is true, then we terminate once `early_term_set` is empty
    std::unordered_set<GRAPH_COMPONENT_ID> early_term_set;
    early_term_set.reserve(16);
    if (terminate_early)
        early_term_set.insert(et_begin, et_end);

    while (!pq.empty() && (!terminate_early || !early_term_set.empty()))
    {
        auto [v_id, d] = pq.top();
        pq.pop();

        // check if `v_id` is up-to-date:
        if (d > dist[v_id])
            continue;

        if (terminate_early)   
            early_term_set.erase(v_id);

        // otherwise, get the corresponding vertex:
        auto* v = gr.get_vertex(v_id);
        for (const auto& [w, e] : gr.get_adjacency_list(v))
        {
            GRAPH_COMPONENT_ID w_id = w->id;

            W new_dist = d + wf(e);
            if (new_dist < dist[w_id])
            {
                dist[w_id] = new_dist;
                prev[w_id] = v_id;
                pq.push({w_id, new_dist});
            }
        }
    }

    return DIJKSTRA_RESULT<W>{dist, prev};
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

inline std::vector<GRAPH_COMPONENT_ID>
dijkstra_path(const std::vector<GRAPH_COMPONENT_ID>& prev, 
                GRAPH_COMPONENT_ID src,
                GRAPH_COMPONENT_ID dst,
                bool reverse_ok)
{
    std::vector<GRAPH_COMPONENT_ID> path;
    path.reserve(4);

    GRAPH_COMPONENT_ID curr{dst};
    while (curr != src)
    {
        path.push_back(curr);
        curr = prev[curr];
    }
    path.push_back(src);

    if (!reverse_ok)
        std::reverse(path.begin(), path.end());

    return path;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

}  // namespace graph
