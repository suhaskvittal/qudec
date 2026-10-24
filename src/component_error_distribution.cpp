// Implementation authored by OpenAI Codex.
#include "component_error_distribution.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
uint64_t mix(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
void hash_combine(uint64_t& h, uint64_t x) { h = mix(h ^ mix(x)); }

stim::Circuit transform(const stim::Circuit& input, const ComponentErrorDistribution& distribution) {
    stim::Circuit out;
    for (const auto& op : input.operations) {
        if (op.gate_type == stim::GateType::REPEAT) {
            out.append_repeat_block(op.repeat_block_rep_count(),
                transform(op.repeat_block_body(input), distribution), op.tag);
            continue;
        }
        // Contributor: OpenAI Codex. Each Pauli product is one fixed check-
        // readout component, so phenomenological MPP probabilities vary by
        // check/basis and remain fixed when that check is repeated in time.
        if (op.gate_type == stim::GateType::MPP && !op.args.empty()) {
            for (size_t start = 0; start < op.targets.size();) {
                size_t end = start + 1;
                while (end < op.targets.size() && op.targets[end].is_combiner()) end += 2;
                uint64_t id = mix(static_cast<uint64_t>(op.gate_type));
                for (unsigned char ch : op.tag) hash_combine(id,ch);
                std::vector<uint32_t> targets;
                for (size_t k = start; k < end; k++) {
                    targets.push_back(op.targets[k].data);
                    hash_combine(id,op.targets[k].data);
                }
                out.safe_append_ua("MPP",targets,distribution.probability(op.args[0],id),op.tag);
                start = end;
            }
            continue;
        }
        bool pair = op.gate_type == stim::GateType::DEPOLARIZE2;
        bool supported = pair || op.gate_type == stim::GateType::DEPOLARIZE1 ||
            op.gate_type == stim::GateType::X_ERROR || op.gate_type == stim::GateType::Y_ERROR ||
            op.gate_type == stim::GateType::Z_ERROR;
        if (!supported) {
            if ((stim::GATE_DATA[op.gate_type].flags & stim::GATE_IS_NOISY) && !op.args.empty())
                throw std::invalid_argument("Nonuniform noise: unsupported channel " +
                    std::string(stim::GATE_DATA[op.gate_type].name));
            out.safe_append(op);
            continue;
        }
        const size_t stride = pair ? 2 : 1;
        for (size_t k = 0; k < op.targets.size(); k += stride) {
            uint64_t id = mix(static_cast<uint64_t>(op.gate_type));
            for (unsigned char ch : op.tag) hash_combine(id, ch);
            std::vector<uint32_t> targets;
            for (size_t j = 0; j < stride; ++j) {
                targets.push_back(op.targets[k+j].data);
                hash_combine(id, targets.back());
            }
            // Stim's DEM representation cannot analyze over-mixing
            // depolarization, although such probabilities are valid circuits.
            const double cap = pair ? 15.0/16 :
                op.gate_type == stim::GateType::DEPOLARIZE1 ? 3.0/4 : 1;
            out.safe_append_ua(stim::GATE_DATA[op.gate_type].name, targets,
                std::min(cap, distribution.probability(op.args[0], id)), op.tag);
        }
    }
    return out;
}
}

void ComponentErrorDistribution::validate() const {
    if (!std::isfinite(log_stddev) || log_stddev < 0 ||
        log_stddev > std::sqrt(std::numeric_limits<double>::max()))
        throw std::invalid_argument("Noise log standard deviation must be finite and nonnegative");
}

double ComponentErrorDistribution::probability(double nominal, uint64_t component_id) const {
    validate();
    if (!std::isfinite(nominal) || nominal < 0 || nominal > 1)
        throw std::invalid_argument("Nominal component probability must be in [0, 1]");
    if (log_stddev == 0 || nominal == 0) return nominal;
    // Box-Muller from stateless SplitMix64 bits: independent of traversal,
    // circuit batching, REPEAT expansion, MPI rank, or standard-library RNG.
    const uint64_t bits1 = mix(component_id ^ mix(seed));
    const uint64_t bits2 = mix(bits1);
    const double u1 = ((bits1 >> 12) + 0.5) * 0x1p-52;
    const double u2 = ((bits2 >> 12) + 0.5) * 0x1p-52;
    const double z = std::sqrt(-2 * std::log(u1)) * std::cos(6.2831853071795864769 * u2);
    const double log_rate = std::log(nominal) + log_stddev*z - 0.5*log_stddev*log_stddev;
    return log_rate >= 0 ? 1 : std::exp(log_rate);
}

stim::Circuit apply_component_error_distribution(
    const stim::Circuit& circuit, const ComponentErrorDistribution& distribution) {
    distribution.validate();
    return distribution.log_stddev == 0 ? circuit : transform(circuit, distribution);
}
