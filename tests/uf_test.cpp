/*
 *  contributor: OpenAI Codex
 * */

#include "decoder/uf.h"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

void
check_predictions(const std::string& model,
                  const std::vector<std::pair<unsigned, unsigned>>& cases)
{
    const stim::DetectorErrorModel dem(model);
    dec::UnionFind decoder(dem);
    dec::SyndromeType syndrome(dem.count_detectors());
    dec::ObsType truth(dem.count_observables());
    for (const auto& [events, prediction] : cases)
    {
        syndrome.clear();
        for (size_t d = 0; d < dem.count_detectors(); d++)
            syndrome[d] = (events >> d) & 1;
        const auto original = syndrome;
        const auto result = decoder.decode(syndrome, truth);
        if (syndrome != original)
            throw std::runtime_error("UnionFind modified the input syndrome");
        for (size_t obs = 0; obs < dem.count_observables(); obs++)
            if (bool(result.flipped_obs[obs]) != bool((prediction >> obs) & 1))
                throw std::runtime_error("UnionFind returned the wrong observable correction");
        if (result.matching_data.frame_flips != result.flipped_obs)
            throw std::runtime_error("UnionFind returned inconsistent frame flips");
    }
}

} // anon

int
main()
{
    // Equal hop counts but unequal likelihoods: the more likely boundary wins,
    // regardless of whether its observable is zero or one.
    check_predictions("error(0.01) D0\nerror(0.2) D0 L0\n", {{1, 1}, {0, 0}, {1, 1}});
    check_predictions("error(0.2) D0\nerror(0.01) D0 L0\n", {{1, 0}});

    // A longer route can grow first because its total weight is smaller.
    check_predictions("error(0.01) D0\nerror(0.2) D0 D1\nerror(0.2) D1 L0\n",
                      {{1, 1}, {0, 0}, {3, 0}});

    // Two active clusters meet twice as quickly as growth from one endpoint.
    check_predictions("error(0.2) D0 D1 L0\nerror(0.01) D0\nerror(0.01) D1\n",
                      {{3, 1}, {1, 0}, {2, 0}, {0, 0}});

    // Multiple independently absorbing boundaries and observable XOR.
    check_predictions("error(0.1) D0 L0\nerror(0.1) D1 L0\n",
                      {{3, 0}, {1, 1}, {2, 1}});
    check_predictions("error(0.1) D0 L0 ^ D1 L1\n", {{3, 3}, {1, 1}, {2, 2}});

    // Negative log-odds and deterministic channels are baseline corrections.
    check_predictions("error(0.8) D0 L0\n", {{0, 0}, {1, 1}});
    check_predictions("error(1) D0 L0\n", {{1, 1}, {1, 1}});
    check_predictions("error(0.8) L0\n", {{0, 1}});

    // Identical channels combine by odd-parity probability, not by summing p.
    // Two p=0.8 channels combine to p=0.32, below the competing p=0.4 edge.
    check_predictions("error(0.8) D0 L0\nerror(0.8) D0 L0\nerror(0.4) D0\n",
                      {{1, 0}, {0, 0}});
    check_predictions("error(0.5) D0 L0\n", {{1, 1}, {0, 0}});

    // A closed graph can decode even parity but must reject odd parity.
    const stim::DetectorErrorModel closed("error(0.1) D0 D1 L0\n");
    dec::UnionFind decoder(closed);
    dec::SyndromeType syndrome(2);
    dec::ObsType truth(1);
    syndrome[0] = 1;
    bool rejected{false};
    try
    {
        decoder.decode(syndrome, truth);
    }
    catch (const std::runtime_error&)
    {
        rejected = true;
    }
    if (!rejected)
        throw std::runtime_error("UnionFind accepted an impossible odd syndrome");
    syndrome[1] = 1;
    if (!decoder.decode(syndrome, truth).flipped_obs[0])
        throw std::runtime_error("UnionFind did not recover after rejecting a syndrome");

    // Hyperedges must be explicitly decomposed, rather than silently discarded.
    rejected = false;
    try
    {
        dec::UnionFind invalid(stim::DetectorErrorModel("error(0.1) D0 D1 D2\n"));
    }
    catch (const std::invalid_argument&)
    {
        rejected = true;
    }
    if (!rejected)
        throw std::runtime_error("UnionFind accepted an undecomposed hyperedge");

    // Observable storage is not limited to a 64-bit mask.
    dec::UnionFind wide(stim::DetectorErrorModel("error(0.1) D0 L70\n"));
    dec::SyndromeType one(1);
    one[0] = 1;
    dec::ObsType wide_truth(71);
    if (!wide.decode(one, wide_truth).flipped_obs[70])
        throw std::runtime_error("UnionFind lost an observable above index 63");

    std::cout << "UnionFind tests passed\n";
}
