/*
 *  contributor: OpenAI Codex
 * */

#ifndef QUDEC_LIBRA_LOCAL_SPLITTING_H
#define QUDEC_LIBRA_LOCAL_SPLITTING_H

#include "libra_splitting.h"
#include <cstdint>

// Adds count-preserving nearby swaps to the exact conditional-failure kernel.
// Nearby means that two raw DEM instructions share at least one detector.
// This changes mixing, not the target distribution or bridge likelihood ratio.
class LibraLocalSplittingModel : public LibraSplittingModel {
public:
    explicit LibraLocalSplittingModel(const stim::DetectorErrorModel& dem)
        : LibraSplittingModel(dem), neighbors_(probabilities.size()) {
        if (probabilities.size() >= UINT32_MAX)
            throw std::invalid_argument("Local splitting channel index overflow");
        std::vector<std::vector<uint32_t>> incidence(num_detectors);
        for (size_t i = 0; i < patterns.size(); ++i)
            for (uint32_t d : patterns[i].detectors)
                incidence[d].push_back(static_cast<uint32_t>(i));
        for (size_t i = 0; i < patterns.size(); ++i) {
            auto& nearby = neighbors_[i];
            for (uint32_t d : patterns[i].detectors)
                nearby.insert(nearby.end(), incidence[d].begin(), incidence[d].end());
            std::sort(nearby.begin(), nearby.end());
            nearby.erase(std::unique(nearby.begin(), nearby.end()), nearby.end());
            nearby.erase(std::remove(nearby.begin(), nearby.end(), i), nearby.end());
        }
    }

    template<class Decide>
    bool step(State& state, double theta, std::mt19937_64& rng, const Decide& decide) const {
        // Each component is a reversible kernel. The fixed mixture is also
        // reversible: birth/death each have probability 1/4, swap has 1/2.
        if (rng() & 1)
            return LibraSplittingModel::step(state, theta, rng, decide);
        ++state.attempts;
        if (state.active.empty())
            return false;
        std::uniform_int_distribution<size_t> choose_active(0, state.active.size() - 1);
        const size_t position = choose_active(rng);
        const size_t removed = state.active[position];
        double forward_sum = 0;
        for (uint32_t i : neighbors_[removed])
            if (state.position[i] < 0)
                forward_sum += odds[i];
        if (!(forward_sum > 0))
            return false;
        std::uniform_real_distribution<double> unit(0, 1);
        double draw = unit(rng) * forward_sum;
        size_t added = SIZE_MAX;
        for (uint32_t i : neighbors_[removed]) {
            if (state.position[i] >= 0 || odds[i] == 0)
                continue;
            added = i;
            draw -= odds[i];
            if (draw < 0)
                break;
        }
        if (added == SIZE_MAX)
            throw std::runtime_error("Local splitting proposal has no candidate");
        double reverse_sum = 0;
        for (uint32_t i : neighbors_[added])
            if (state.position[i] < 0 || i == removed)
                reverse_sum += odds[i];
        if (!(reverse_sum > 0))
            throw std::runtime_error("Local splitting reverse proposal has zero mass");

        // q(S->T)=odds_added/(K*Z_removed), and
        // q(T->S)=odds_removed/(K*Z_added). The target odds ratio cancels
        // their odds terms, leaving acceptance min(1,Z_removed/Z_added).
        if (unit(rng) >= std::min(1.0, forward_sum / reverse_sum))
            return false;
        ++state.physical_accepts;
        flip(state, removed);
        flip(state, added);
        const auto labels = decide(state.syndrome, state.truth);
        if (!labels.any()) {
            flip(state, added);
            flip(state, removed);
            return false;
        }
        ++state.event_accepts;
        state.labels = labels;
        state.active[position] = added;
        state.position[removed] = -1;
        state.position[added] = static_cast<int32_t>(position);
        state.active_odds += odds[added] - odds[removed];
        return true;
    }

private:
    std::vector<std::vector<uint32_t>> neighbors_;

    void flip(State& state, size_t id) const {
        for (uint32_t d : patterns[id].detectors)
            state.syndrome[d] ^= 1;
        state.truth[0] ^= patterns[id].observable;
    }
};

#endif
