/*
 *  contributor: OpenAI Codex
 *  date:        13 October 2026
 * */

#include "decoder/libra/assignment.h"
#include <bit>
#include <iostream>
#include <random>

using dec::libra_assignment::Pair;

int64_t exhaustive_cost(const std::vector<double>& singles, const std::vector<Pair>& pairs) {
    const size_t n = singles.size();
    std::vector<int64_t> memo(size_t{1} << n, -1);
    memo[0] = 0;
    auto visit = [&](auto&& self, unsigned mask) -> int64_t {
        if (memo[mask] >= 0)
            return memo[mask];
        const unsigned a = std::countr_zero(mask);
        const unsigned remaining = mask ^ (1U << a);
        int64_t best = dec::libra_assignment::quantize(singles[a]) + self(self, remaining);
        for (const auto& pair : pairs) {
            const unsigned b = pair.a == a ? pair.b : pair.b == a ? pair.a : unsigned(n);
            if (b < n && (remaining & (1U << b)))
                best = std::min(best, dec::libra_assignment::quantize(pair.weight) +
                                      self(self, remaining ^ (1U << b)));
        }
        return memo[mask] = best;
    };
    return visit(visit, (1U << n) - 1);
}

void check(const std::vector<double>& singles, const std::vector<Pair>& pairs) {
    const auto result = dec::libra_assignment::solve(singles, pairs);
    if (result.quantized_cost != exhaustive_cost(singles, pairs))
        throw std::runtime_error("PyMatching assignment differs from exhaustive optimum");
    std::vector<unsigned> covered(singles.size(), 0);
    int64_t cost = 0;
    for (size_t i = 0; i < singles.size(); ++i) {
        if (result.mate[i] < 0) {
            ++covered[i];
            cost += dec::libra_assignment::quantize(singles[i]);
        } else if (size_t(result.mate[i]) >= singles.size() ||
                   result.mate[result.mate[i]] != int32_t(i)) {
            throw std::runtime_error("PyMatching assignment mate invariant failed");
        }
    }
    for (size_t index : result.pair_indices) {
        const auto& pair = pairs.at(index);
        ++covered[pair.a];
        ++covered[pair.b];
        cost += dec::libra_assignment::quantize(pair.weight);
        if (result.mate[pair.a] != int32_t(pair.b))
            throw std::runtime_error("PyMatching assignment pair identity failed");
    }
    if (cost != result.quantized_cost ||
        std::find_if(covered.begin(), covered.end(), [](unsigned c) {return c != 1;}) != covered.end())
        throw std::runtime_error("PyMatching assignment cover invariant failed");
}

int main() {
    // Unpenalized parity decoding would incorrectly choose all three star
    // edges, cost 3. The exact-once assignment costs 21.
    check({10, 10, 10, 10}, {{0,1,1}, {0,2,1}, {0,3,1}});
    check({0, 0, 0, 0}, {{0,1,0}, {1,2,0}, {2,3,0}, {3,0,0}});
    check({4, 5, 6}, {{0,1,9}, {1,2,1}, {1,2,2}, {0,2,0.5}});
    check({}, {});
    std::mt19937_64 rng(20261007);
    std::uniform_real_distribution<double> weight(0.01, 20.0);
    for (unsigned trial = 0; trial < 3000; ++trial) {
        const unsigned n = 1 + rng() % 12;
        std::vector<double> singles(n);
        for (double& w : singles)
            w = trial % 3 ? weight(rng) : double(rng() % 8);
        std::vector<Pair> pairs;
        for (unsigned a = 0; a < n; ++a)
            for (unsigned b = a + 1; b < n; ++b)
                if (rng() % 3)
                    pairs.push_back({a, b, trial % 3 ? weight(rng) : double(rng() % 8)});
        check(singles, pairs);
    }
    bool overflow_rejected = false;
    try {
        dec::libra_assignment::solve(std::vector<double>(9, 1000), {{0,1,1900}});
    } catch (const std::overflow_error&) {
        overflow_rejected = true;
    }
    if (!overflow_rejected)
        throw std::runtime_error("PyMatching assignment overflow guard failed");
    std::cout << "3000 exhaustive assignment oracles, star reduction, ties, duplicate pairs, and overflow passed\n";
}
