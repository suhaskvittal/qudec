/* Author: OpenAI Codex. Chain connectivity uses shared flipped detectors,
 * including flags. DEM separators do not split correlated fault events.
 */
#include "error_chain.h"
#include "error_chain_display.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>
#include <memory>
#include <ostream>
#include <random>
#include <set>
#include <stdexcept>
#include <vector>

namespace
{
struct Event
{
    size_t original_id;
    double probability;
    std::vector<size_t> detectors;
    std::vector<size_t> observables;
};
void parity_reduce(std::vector<size_t>& ids)
{
    std::sort(ids.begin(),ids.end());
    std::vector<size_t> out;
    for (size_t i = 0; i < ids.size();)
    {
        size_t j = i+1;
        while (j < ids.size() && ids[j] == ids[i]) j++;
        if ((j-i)%2) out.push_back(ids[i]);
        i = j;
    }
    ids = std::move(out);
}
void increment(uint64_t& count)
{
    if (count == std::numeric_limits<uint64_t>::max())
        throw std::overflow_error("Chain counter exceeds uint64_t");
    count++;
}
}

std::string error_chain_combination_count(size_t n,size_t k)
{
    if (k > n) return "0";
    k = std::min(k,n-k);
    std::string digits = "1"; // Least-significant digit first.
    for (size_t i = 1; i <= k; i++)
    {
        __uint128_t carry = 0;
        for (char& digit : digits)
        {
            carry += static_cast<unsigned>(digit-'0')*static_cast<__uint128_t>(n-k+i);
            digit = '0'+carry%10;
            carry /= 10;
        }
        while (carry) {digits.push_back('0'+carry%10); carry /= 10;}
        __uint128_t remainder = 0;
        for (size_t j = digits.size(); j-- > 0;)
        {
            remainder = remainder*10+static_cast<unsigned>(digits[j]-'0');
            digits[j] = '0'+remainder/i;
            remainder %= i;
        }
        while (digits.size()>1 && digits.back() == '0') digits.pop_back();
    }
    std::reverse(digits.begin(),digits.end());
    return digits;
}

ErrorChainSummary evaluate_error_chains(const stim::DetectorErrorModel& dem,
    const ErrorChainConfig& config,const ErrorChainDecode& decode,std::ostream& out)
{
    if (!config.length || (!config.exhaustive && !config.samples))
        throw std::invalid_argument("Chain length and sample count must be positive");
    std::vector<Event> events;
    size_t original_id = 0;
    dem.iter_flatten_error_instructions([&](const auto& inst)
    {
        const auto id = original_id++;
        const double p = inst.arg_data[0];
        if (!std::isfinite(p) || p < 0 || p > 1) throw std::invalid_argument("Invalid DEM probability");
        if (!p) return;
        Event event{.original_id=id,.probability=p};
        for (const auto& target : inst.target_data)
        {
            if (target.is_relative_detector_id()) event.detectors.push_back(target.val());
            else if (target.is_observable_id()) event.observables.push_back(target.val());
        }
        parity_reduce(event.detectors);
        parity_reduce(event.observables);
        events.push_back(std::move(event));
    });
    if (config.length > events.size()) throw std::invalid_argument("Chain length exceeds positive-probability DEM events");
    std::vector<std::vector<size_t>> incidence(dem.count_detectors()), adjacency(events.size());
    for (size_t i = 0; i < events.size(); i++)
        for (auto detector : events[i].detectors) incidence.at(detector).push_back(i);
    for (size_t i = 0; i < events.size(); i++)
    {
        auto& neighbors = adjacency[i];
        for (auto detector : events[i].detectors)
            for (auto j : incidence[detector]) if (i != j) neighbors.push_back(j);
        std::sort(neighbors.begin(),neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(),neighbors.end()),neighbors.end());
    }
    std::vector<uint8_t> visited(events.size());
    std::vector<size_t> eligible;
    for (size_t root = 0; root < events.size(); root++)
    {
        if (visited[root]) continue;
        std::vector<size_t> component{root};
        visited[root] = 1;
        for (size_t i = 0; i < component.size(); i++)
            for (auto next : adjacency[component[i]])
                if (!visited[next]) {visited[next] = 1; component.push_back(next);}
        if (component.size() >= config.length) eligible.insert(eligible.end(),component.begin(),component.end());
    }
    out << "DEM events: " << events.size() << " (positive probability; IDs are zero-based original flattened instructions)\n";
    if (config.exhaustive)
        out << "Total combinations before connectivity filtering: "
            << error_chain_combination_count(events.size(),config.length) << std::endl;
    else
        out << "Random chains: " << config.samples << " (seed/frontier growth, with replacement; not uniform over connected subsets)\n";
    if (eligible.empty()) throw std::invalid_argument("No connected chain of the requested length exists");

    std::vector<uint8_t> selected(events.size());
    auto connected = [&](const std::vector<size_t>& pattern)
    {
        for (auto id : pattern) selected[id] = 1;
        std::vector<size_t> reached{pattern[0]};
        selected[pattern[0]] = 2;
        for (size_t j = 0; j < reached.size(); j++)
            for (auto id : adjacency[reached[j]])
                if (selected[id] == 1) {selected[id] = 2; reached.push_back(id);}
        for (auto id : pattern) selected[id] = 0;
        return reached.size() == pattern.size();
    };
    ErrorChainSummary summary;
    std::unique_ptr<ErrorChainDisplay> display;
    // Author: OpenAI Codex. Preserve the requested OR convention for flags
    // in chain stress tests, with XOR for checks and logical effects.
    std::set<uint64_t> detector_ids;
    for (size_t d = 0; d < dem.count_detectors(); d++) detector_ids.insert(d);
    std::vector<bool> is_flag(dem.count_detectors(),false);
    for (const auto& [id,c] : dem.get_detector_coordinates(detector_ids))
        is_flag[id] = c.size() == 6 && c[5] == 1;
    dec::SyndromeType syndrome(dem.count_detectors()), truth(dem.count_observables());
    auto test = [&](const std::vector<size_t>& pattern)
    {
        syndrome.clear(); truth.clear();
        for (auto id : pattern)
        {
            for (auto d : events[id].detectors)
                if (is_flag[d]) syndrome[d] = 1;
                else syndrome[d] ^= 1;
            for (auto obs : events[id].observables) truth[obs] ^= 1;
        }
        const auto result = decode(syndrome,truth);
        bool failed = false;
        for (size_t obs = 0; obs < dem.count_observables(); obs++)
            failed |= result.flipped_obs[obs] != truth[obs];
        increment(summary.tested);
        if (failed) increment(summary.failures);
        if ((failed && summary.failures <= config.print_failures) || summary.tested <= config.print_chains)
        {
            if (!display) display = std::make_unique<ErrorChainDisplay>(dem,config.code,config.distance);
            if (failed) out << "Failure " << summary.failures << " at chain " << summary.tested;
            else out << "Chain " << summary.tested << " (decoded correctly)";
            out << ": events=";
            for (auto id : pattern) out << events[id].original_id << ',';
            out << " detectors=";
            for (size_t d = 0; d < dem.count_detectors(); d++) if (syndrome[d]) out << 'D' << d << ',';
            out << " truth=";
            for (size_t obs = 0; obs < dem.count_observables(); obs++) out << truth[obs];
            out << " predicted=";
            for (size_t obs = 0; obs < dem.count_observables(); obs++) out << result.flipped_obs[obs];
            out << '\n';
            for (auto id : pattern)
            {
                out << "  Event " << events[id].original_id << " (p=" << events[id].probability << ") flips detectors {";
                bool first = true;
                for (auto detector : events[id].detectors)
                {
                    if (!first) out << ", ";
                    out << display->detector_label(detector); first = false;
                }
                out << "}, observables {"; first = true;
                for (auto obs : events[id].observables)
                {
                    if (!first) out << ", ";
                    out << 'L' << obs; first = false;
                }
                out << "}\n";
            }
            display->draw(syndrome,out);
            if (!result.debug_info.empty()) out << "Decoder decisions:\n" << result.debug_info;
        }
        if (config.print_progress && summary.tested%10000 == 0)
            out << "Progress: tested=" << summary.tested << " failures=" << summary.failures << std::endl;
    };
    if (config.exhaustive)
    {
        std::vector<size_t> pattern(config.length);
        std::iota(pattern.begin(),pattern.end(),0);
        while (true)
        {
            increment(summary.combinations_visited);
            if (connected(pattern)) test(pattern);
            size_t j = pattern.size();
            while (j && pattern[j-1] == events.size()-pattern.size()+j-1) j--;
            if (!j) break;
            pattern[j-1]++;
            for (size_t i = j; i < pattern.size(); i++) pattern[i] = pattern[i-1]+1;
        }
    }
    else
    {
        std::mt19937_64 rng(config.seed);
        for (uint64_t trial = 0; trial < config.samples; trial++)
        {
            std::vector<size_t> pattern{eligible[std::uniform_int_distribution<size_t>(0,eligible.size()-1)(rng)]};
            std::vector<size_t> frontier;
            selected[pattern[0]] = 1;
            auto add_frontier = [&](size_t id)
            {
                for (auto next : adjacency[id])
                    if (!selected[next]) {selected[next] = 2; frontier.push_back(next);}
            };
            add_frontier(pattern[0]);
            while (pattern.size() < config.length)
            {
                if (frontier.empty()) throw std::logic_error("Connected component exhausted");
                const size_t index = std::uniform_int_distribution<size_t>(0,frontier.size()-1)(rng);
                const size_t next = frontier[index];
                frontier[index] = frontier.back(); frontier.pop_back();
                selected[next] = 1; pattern.push_back(next); add_frontier(next);
            }
            for (auto id : frontier) selected[id] = 0;
            for (auto id : pattern) selected[id] = 0;
            std::sort(pattern.begin(),pattern.end());
            test(pattern);
        }
    }
    out << "Chains tested: " << summary.tested << "\nLogical failures: " << summary.failures
        << "\nChain failure fraction: " << std::setprecision(10)
        << static_cast<double>(summary.failures)/summary.tested
        << " (connected-chain stress test, not IID logical error rate)\n";
    if (config.exhaustive) out << "Combinations visited: " << summary.combinations_visited << '\n';
    return summary;
}
