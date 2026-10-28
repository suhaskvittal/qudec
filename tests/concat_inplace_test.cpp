/* Author: OpenAI Codex. Persistent graph updates versus fresh PyMatching graphs. */
#include "decoder/concat_mwpm/matching.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>

namespace
{
void require(bool ok,const char* message)
{
    if (!ok) throw std::runtime_error(message);
}
void tests(size_t ensemble)
{
    using namespace dec::concat;
    Model topology{{0,1,2,0,1,2},{
        {.1,{0},0,BoundaryHint::Spatial},
        {.1,{0,1,2},0,BoundaryHint::Spatial},
        {.1,{0,3},0,BoundaryHint::Temporal},
        {.1,{0,4,5},0,BoundaryHint::Spatial},
        {.1,{1,2},0,BoundaryHint::Spatial},
        {.1,{1,4},0,BoundaryHint::Temporal},
        {.1,{2,5},0,BoundaryHint::Temporal},
        {.1,{3,4,5},0,BoundaryHint::Spatial}}};
    Decoder persistent(topology,false,{ensemble,1729,1});
    std::mt19937_64 rng(1729);
    for (size_t trial=0;trial<60;trial++)
    {
        Model active{topology.colors,{},trial%2==0};
        std::vector<uint8_t> syndrome(topology.colors.size());
        for (const auto& original : topology.events)
        {
            if (trial%10==0 || rng()%3==0) continue;
            auto event=original;
            event.probability=std::pow(10.0,-1-int(rng()%8));
            if (trial%11==0) event.probability=1;
            event.observable=rng()%2;
            event.flag_supported=rng()%2;
            if (rng()%2) for(auto v:event.detectors) syndrome[v]^=1;
            active.events.push_back(std::move(event));
        }
        persistent.update(active);
        Decoder rebuilt(active,false,{ensemble,1729,1});
        require(persistent.audit_weights()==rebuilt.audit_weights(),"Active graph edge counts differ");
        const auto actual=persistent.decode_candidates(syndrome);
        const auto expected=rebuilt.decode_candidates(syndrome);
        require(actual.size()==expected.size(),"In-place candidate count changed");
        for(size_t i=0;i<actual.size();i++)
        {
            const auto& a=actual[i];const auto& b=expected[i];
            require(a.observable==b.observable && std::abs(a.weight-b.weight)<1e-10,
                    "In-place correction frame/weight differs from a fresh graph");
            require(a.events==b.events && a.restricted_edges==b.restricted_edges &&
                    a.second_stage_paths==b.second_stage_paths && a.virtual_detectors==b.virtual_detectors,
                    "In-place matching paths differ from a fresh graph");
            require(a.color==b.color && a.member==b.member,
                    "In-place ensemble provenance changed");
        }
    }
}
}
int main()
{
    try
    {
        for(size_t ensemble:{1,5}) tests(ensemble);
        std::cout<<"In-place graph updates reproduce fresh matching graphs\n";
    }
    catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
