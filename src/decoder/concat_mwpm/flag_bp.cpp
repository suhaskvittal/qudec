/*
 *  author: OpenAI Codex
 *  date:   27 October 2026
 * */
#include "decoder/concat_mwpm/flag_bp.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dec::concat
{
FlagBP::FlagBP(size_t num_detectors,std::vector<FlagBPVariable> variables)
    :variables_(std::move(variables)),checks_(num_detectors),variable_edges_(variables_.size()),
     star_posterior_(variables_.size()),prior_llr_(variables_.size())
{
    for (size_t i=0;i<variables_.size();i++)
    {
        const auto& variable=variables_[i];
        if (!std::isfinite(variable.probability) || variable.probability<0 || variable.probability>1)
            throw std::invalid_argument("Invalid auxiliary BP prior");
        stars_ &= variable.flags.size()<=1;
        const double p=std::clamp(variable.probability,1e-15,1-1e-15);
        prior_llr_[i]=std::log((1-p)/p);
        for (auto flag:variable.flags)
        {
            if (flag>=num_detectors) throw std::invalid_argument("Invalid auxiliary BP detector");
            variable_edges_[i].push_back(owners_.size());
            checks_[flag].push_back(owners_.size());
            owners_.push_back(i);
        }
        star_posterior_[i]={p,p};
    }
    if (!stars_) return;
    // Prefix/suffix products avoid division and handle a p=1/2 neighbor.
    for (const auto& edges:checks_) if (!edges.empty())
    {
        std::vector<double> prefix(edges.size()+1,1),suffix(edges.size()+1,1);
        for (size_t j=0;j<edges.size();j++)
            prefix[j+1]=prefix[j]*std::tanh(prior_llr_[owners_[edges[j]]]*.5);
        for (size_t j=edges.size();j-->0;)
            suffix[j]=suffix[j+1]*std::tanh(prior_llr_[owners_[edges[j]]]*.5);
        for (size_t j=0;j<edges.size();j++) for (size_t outcome=0;outcome<2;outcome++)
        {
            const auto id=owners_[edges[j]];
            const double product=(outcome ? -1:1)*prefix[j]*suffix[j+1];
            const double llr=prior_llr_[id]+2*std::atanh(std::clamp(product,-1+1e-14,1-1e-14));
            star_posterior_[id][outcome]=1/(1+std::exp(std::clamp(llr,-700.0,700.0)));
        }
    }
}

std::vector<double>
FlagBP::posterior(const std::vector<uint8_t>& syndrome) const
{
    if (syndrome.size()!=checks_.size()) throw std::invalid_argument("Wrong auxiliary BP syndrome length");
    std::vector<double> result(variables_.size());
    if (stars_)
    {
        for (size_t i=0;i<variables_.size();i++)
            result[i]=star_posterior_[i][variables_[i].flags.empty() ? 0:bool(syndrome[variables_[i].flags[0]])];
        return result;
    }
    // General correlated auxiliary-factor fallback. Converge marginals, rather than
    // stopping merely because a hard decision satisfies the observed parity.
    constexpr size_t max_iterations=30;
    constexpr double damping=.5,tolerance=1e-10;
    std::vector<double> vc(owners_.size()),cv(owners_.size(),0),next(owners_.size()),total(variables_.size());
    for (size_t e=0;e<owners_.size();e++) vc[e]=prior_llr_[owners_[e]];
    std::vector<double> prefix,suffix;
    for (size_t iteration=0;iteration<max_iterations;iteration++)
    {
        double change=0;
        for (size_t d=0;d<checks_.size();d++)
        {
            const auto& edges=checks_[d];
            if (edges.empty()) continue;
            prefix.assign(edges.size()+1,1);suffix.assign(edges.size()+1,1);
            for (size_t j=0;j<edges.size();j++) prefix[j+1]=prefix[j]*std::tanh(vc[edges[j]]*.5);
            for (size_t j=edges.size();j-->0;) suffix[j]=suffix[j+1]*std::tanh(vc[edges[j]]*.5);
            for (size_t j=0;j<edges.size();j++)
            {
                const double product=(syndrome[d] ? -1:1)*prefix[j]*suffix[j+1];
                next[edges[j]]=damping*cv[edges[j]]+(1-damping)*2*std::atanh(std::clamp(product,-1+1e-14,1-1e-14));
            }
        }
        cv.swap(next);
        for (size_t i=0;i<variables_.size();i++)
        {
            total[i]=prior_llr_[i];
            for (auto edge:variable_edges_[i]) total[i]+=cv[edge];
            const double probability=1/(1+std::exp(std::clamp(total[i],-700.0,700.0)));
            change=std::max(change,std::abs(probability-result[i]));result[i]=probability;
            for (auto edge:variable_edges_[i]) vc[edge]=std::clamp(total[i]-cv[edge],-50.0,50.0);
        }
        if (iteration && change<tolerance) break;
    }
    return result;
}
} // namespace dec::concat
