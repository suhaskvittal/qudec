/*
 *  author: OpenAI Codex
 *  date:   22 October 2026
 *
 *  Single-pass correlation by fractional discounts of original log-odds weights.
 * */
#ifndef DECODER_CORRELATED_RESTRICTION_MATCHING_h
#define DECODER_CORRELATED_RESTRICTION_MATCHING_h

#include "decoder/restriction/lattice.h"
#include <set>
#include <string>

namespace dec::correlated_restriction
{
using namespace restriction;
using EdgeSet = std::set<edge_type>;
using State = std::array<EdgeSet,3>;

std::vector<Assignment> decompose_edges(const std::vector<Vertex>&, const syndrome_type&,
    const State&, const std::array<std::map<edge_type,edge_type>,3>&);

struct Event
{
    double probability;
    std::vector<det_id_type> checks;
    uint64_t observable;
    std::array<edge_type,3> projection;
};

struct Outcome
{
    uint64_t observable{0};
    double weight{0};
    size_t iterations{0};
    size_t period{0};
    size_t lift_failures{0};
    std::string termination;
    std::vector<std::string> lifting_errors;
    bool chose_initial{false};
    State state;
};

class Decoder
{
public:
    Decoder(uint32_t distance, uint32_t rounds, double p, bool correlated);
    Outcome decode(const syndrome_type&);
    // Initial matching and one synchronous update; no baseline safeguard.
    std::array<Outcome,2> decode_single_pass(const syndrome_type&);
    Outcome inspect_lift(const syndrome_type&, const State&) const;
    syndrome_type syndrome(const std::vector<size_t>&) const;
    uint64_t observable(const std::vector<size_t>&) const;
    std::vector<std::array<double,3>> conditional_probabilities(const State&) const;

    const Lattice lattice;
    const uint32_t rounds;
    const size_t num_detectors;
    std::vector<Event> events;
    std::vector<Vertex> vertices;
    std::array<det_id_type,3> boundaries;
    double probability_floor{1e-12};

private:
    bool correlated_;
    restriction::RestrictionDecoder lifter_;
    std::array<pm::Mwpm,3> initial_matchers_;
    std::array<std::map<edge_type,edge_type>,3> colored_edges_;
    std::array<std::map<edge_type,std::vector<size_t>>,3> incident_events_;
    State match(const syndrome_type&, const std::vector<std::array<double,3>>*, bool);
    std::vector<Assignment> assignments(const syndrome_type&, const State&) const;
    Outcome lift(const syndrome_type&, const State&) const;
    double score_history(const syndrome_type&, const restriction::Result&) const;
};
} // namespace dec::correlated_restriction
#endif
