/*
 *  author: Suhas Vittal
 *  date:   24 September 2026
 * */

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <unordered_map>

namespace pp
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

template <class G> AffinityResult
measure_affinity(const G& gr, const std::vector<hg::id_type>& syndrome)
{
    return aff::impl_dijkstra(gr, syndrome);
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace aff
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

template <class G> AffinityResult
// OpenAI GPT-6: Convert shortest negative-log paths into pairwise affinities.
impl_dijkstra(const G& gr, const std::vector<hg::id_type>& s)
{
    AffinityResult a(s.size()*s.size());
    std::vector<double> dist(gr.N());

    std::vector<bool> s_membership(gr.N(), false);
    std::vector<bool> settled(gr.N(), false);
    for (auto x : s)
        s_membership[x] = true;

    // `pmap` is an accumulation structure in the traversal segment of Dijkstra's.
    std::unordered_map<hg::id_type, double> pmap;
    pmap.reserve(16);

    /*
     * Priority queue definition:
     * */
    struct pq_entry { hg::id_type id; double w; };
    auto pq_cmp = [] (const auto& a, const auto& b) { return a.w > b.w; };
    using queue_type = std::priority_queue<pq_entry, std::vector<pq_entry>, decltype(pq_cmp)>;

    // Run Dijkstra's for each detector in `s`
    size_t ii{0};
    for (size_t i = 0; i < s.size(); i++)
    {
        constexpr auto INF = std::numeric_limits<double>::max();

        size_t s_rem{s.size()};  // number of syndrome bits we are still computing distance for

        std::fill(dist.begin(), dist.end(), INF);
        std::fill(settled.begin(), settled.end(), false);
        dist[s[i]] = 0.0;
        queue_type q(pq_cmp);
        q.push( pq_entry{s[i], 0.0} );
        while (q.size() > 0 && s_rem > 0)
        {
            auto entry = q.top();
            q.pop();
            if (entry.w > dist[entry.id] || settled[entry.id])
                continue;
            // OpenAI GPT-6: Count each target only once when equal-cost paths exist.
            settled[entry.id] = true;
            if (s_membership[entry.id])
                s_rem--;

            auto r = entry.id;
            gr.for_each_edge_incident_to({r},
                    [&] (auto e)
                    {
                        const double pr = gr.e(e).pr;
                        for (auto s : gr.support(e))
                            if (s != r)
                                pmap[s] += pr;
                    });
            for (auto [s, pr] : pmap)
            {
                const double w = -std::log(pr) + entry.w;
                if (w < dist[s])
                {
                    q.push(pq_entry{s, w});
                    dist[s] = w;
                }
            }
            pmap.clear();
        }

        // once Dijkstra's terminates, copy contents of `dist` to `a`
        for (size_t j = 0; j < s.size(); j++)
            a[ii++] = dist[s[j]] == INF ? 0.0 : std::exp(-dist[s[j]]);
    }
    return a;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace aff

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace pp
