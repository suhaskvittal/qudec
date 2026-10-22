/* Export shared-shot SI1000 syndrome, decoder labels, and detector geometry. */
#include "circuit_generator.h"
#include "decoder/pymatching.h"
#include "preprocessor/affinity.h"

#include <stim.h>
#include <stim/simulators/dem_sampler.h>
#include <stim/util_top/circuit_to_dem.h>

#include <charconv>
#include <bit>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct Options {
    uint32_t distance = 7;
    uint32_t rounds = 0;
    bool rounds_were_set = false;
    double p = 1e-3;
    uint64_t shots = 10000;
    uint64_t seed = 42;
    fs::path output_dir;
};

template <typename T>
T parse_unsigned(const std::string &value, const std::string &option) {
    if (value.empty() || value.front() == '-')
        throw std::invalid_argument(option + " must be a nonnegative integer");
    T result{};
    const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), result);
    if (error == std::errc::result_out_of_range)
        throw std::invalid_argument(option + " is out of range");
    if (error != std::errc{} || end != value.data() + value.size())
        throw std::invalid_argument(option + " must be an integer: " + value);
    return result;
}

double parse_double(const std::string &value, const std::string &option) {
    double result{};
    const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), result);
    if (error == std::errc::result_out_of_range || !std::isfinite(result))
        throw std::invalid_argument(option + " must be finite and in range");
    if (error != std::errc{} || end != value.data() + value.size())
        throw std::invalid_argument(option + " must be a number: " + value);
    return result;
}

Options parse_args(int argc, char **argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        if (i + 1 >= argc) throw std::invalid_argument("missing value for " + key);
        const std::string value = argv[++i];
        if (key == "--distance") o.distance = parse_unsigned<uint32_t>(value, key);
        else if (key == "--rounds") {
            o.rounds = parse_unsigned<uint32_t>(value, key);
            o.rounds_were_set = true;
        }
        else if (key == "--p") o.p = parse_double(value, key);
        else if (key == "--shots") o.shots = parse_unsigned<uint64_t>(value, key);
        else if (key == "--seed") o.seed = parse_unsigned<uint64_t>(value, key);
        else if (key == "--output-dir") o.output_dir = value;
        else throw std::invalid_argument("unknown option: " + key);
    }
    if (o.output_dir.empty()) throw std::invalid_argument("--output-dir is required");
    if (!o.rounds_were_set) o.rounds = o.distance;
    if (o.distance < 3 || o.distance % 2 == 0 || o.rounds == 0 || o.shots == 0 ||
        !std::isfinite(o.p) || !(o.p > 0 && o.p < 0.5))
        throw std::invalid_argument("distance must be odd >= 3; rounds/shots positive; 0 < p < 0.5");
    if (fs::exists(o.output_dir)) throw std::invalid_argument("output directory already exists: " + o.output_dir.string());
    return o;
}

void finish_output(std::ofstream &out, const fs::path &path) {
    out.flush();
    if (!out) throw std::runtime_error("failed flushing " + path.string());
    out.close();
    if (!out) throw std::runtime_error("failed closing " + path.string());
}

void write_u64_le(std::ostream &out, uint64_t value) {
    for (unsigned i = 0; i < 8; ++i) out.put(static_cast<char>(value >> (8 * i)));
}

void write_f64(std::ostream &out, double value) {
    static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559);
    write_u64_le(out, std::bit_cast<uint64_t>(value));
}

std::ofstream binary_file(const fs::path &path) {
    std::ofstream f(path, std::ios::binary | std::ios::trunc);
    if (!f) throw std::runtime_error("failed to create " + path.string());
    return f;
}

void write_graph_edges(const stim::DetectorErrorModel &dem, const fs::path &path) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("failed to create " + path.string());
    out << "detector_u,detector_v,weight_natural_log,boundary_class\n" << std::setprecision(17);
    dem.iter_flatten_error_instructions([&](const stim::DemInstruction &instruction) {
        const double pr = instruction.arg_data[0];
        if (pr <= 0 || pr >= 1) return;
        instruction.for_separated_targets([&](const auto &group) {
            std::vector<uint64_t> detectors;
            bool logical = false;
            for (const auto &target : group) {
                if (target.is_relative_detector_id()) detectors.push_back(target.val());
                if (target.is_observable_id()) logical = !logical;
            }
            if (detectors.empty() || detectors.size() > 2) return;
            const double weight = std::log((1.0 - pr) / pr);
            const int boundary_class = detectors.size() == 1 ? static_cast<int>(logical) : -1;
            out << detectors[0] << ',';
            if (detectors.size() == 2) out << detectors[1];
            else out << -1;
            out << ',' << weight << ',' << boundary_class << '\n';
        });
    });
    finish_output(out, path);
}

void write_coordinates(const stim::Circuit &circuit, size_t n, const fs::path &path) {
    std::set<uint64_t> ids;
    for (size_t i = 0; i < n; ++i) ids.insert(i);
    const auto coords = circuit.get_detector_coordinates(ids);
    std::ofstream out(path);
    if (!out) throw std::runtime_error("failed to create " + path.string());
    out << "detector_id,x,y,t\n" << std::setprecision(17);
    for (size_t i = 0; i < n; ++i) {
        out << i;
        const auto it = coords.find(i);
        for (size_t k = 0; k < 3; ++k) {
            out << ',';
            if (it != coords.end() && k < it->second.size()) out << it->second[k];
        }
        out << '\n';
    }
    finish_output(out, path);
}

void write_metadata(const Options &o, size_t detectors, const fs::path &path) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("failed to create " + path.string());
    out << std::setprecision(17)
        << "{\n  \"schema\": 1,\n  \"distance\": " << o.distance
        << ",\n  \"rounds\": " << o.rounds << ",\n  \"p\": " << o.p
        << ",\n  \"shots\": " << o.shots << ",\n  \"seed\": " << o.seed
        << ",\n  \"memory_basis\": \"Z\",\n  \"include_opposite_basis_detectors\": true"
        << ",\n  \"num_detectors\": " << detectors
        << ",\n  \"sampling\": \"Stim DemSampler from undecomposed circuit DEM; decoders use decomposed DEM\""
        << ",\n  \"binary_format\": \"little-endian float64 where applicable\""
        << ",\n  \"syndromes.u8\": {\"shape\": [" << o.shots << ", " << detectors << "], \"dtype\": \"uint8\"}"
        << ",\n  \"labels.f64\": {\"shape\": [" << o.shots << ", 6], \"columns\": [\"unsigned_complementary_gap_db\", \"cluster_gap_db\", \"logical_truth\", \"ordinary_prediction\", \"ordinary_decode_ns\", \"cluster_decode_ns\"]}"
        << ",\n  \"affinity.f64\": {\"shape\": [" << detectors + 1 << ", " << detectors + 1 << "], \"last_index\": \"shared boundary\"}"
        << ",\n  \"hops.u8\": {\"shape\": [" << detectors + 1 << ", " << detectors + 1 << "], \"last_index\": \"shared boundary\", \"saturated_value\": 31}"
        << ",\n  \"graph_edges.csv\": {\"boundary_class\": \"0 or 1 for v=-1; -1 for detector-detector edges\", \"weight\": \"natural-log odds\", \"rows\": \"unmerged separated DEM components; repeated pairs are not XOR-merged\"}\n}\n";
    finish_output(out, path);
}

} // namespace

int main(int argc, char **argv) {
    try {
        const Options o = parse_args(argc, argv);
        std::error_code dir_error;
        if (!fs::create_directory(o.output_dir, dir_error)) {
            if (dir_error) throw std::runtime_error("could not create output directory: " + dir_error.message());
            throw std::runtime_error("output directory already exists: " + o.output_dir.string());
        }
        const auto circuit = sc_si1000(o.distance, o.rounds, o.p, false, true);
        const auto raw_dem = stim::circuit_to_dem(circuit, {.decompose_errors = false});
        const auto matching_dem = stim::circuit_to_dem(circuit, {.decompose_errors = true});
        const size_t n = raw_dem.count_detectors();
        if (matching_dem.count_detectors() != n || raw_dem.count_observables() != 1 || matching_dem.count_observables() != 1)
            throw std::runtime_error("unexpected detector or observable layout");

        { const auto path = o.output_dir / "circuit.stim"; std::ofstream f(path); if (!f) throw std::runtime_error("failed to create " + path.string()); f << circuit; finish_output(f, path); }
        { const auto path = o.output_dir / "raw.dem"; std::ofstream f(path); if (!f) throw std::runtime_error("failed to create " + path.string()); f << raw_dem; finish_output(f, path); }
        { const auto path = o.output_dir / "decomposed.dem"; std::ofstream f(path); if (!f) throw std::runtime_error("failed to create " + path.string()); f << matching_dem; finish_output(f, path); }
        write_coordinates(circuit, n, o.output_dir / "coords.csv");
        write_graph_edges(matching_dem, o.output_dir / "graph_edges.csv");

        auto affinity_graph = pp::aff::build_hypergraph_from_dem(raw_dem);
        std::vector<hg::id_type> all_ids;
        all_ids.reserve(n + 1);
        for (size_t i = 0; i <= n; ++i) all_ids.push_back(static_cast<hg::id_type>(i));
        const auto affinity = pp::measure_affinity(affinity_graph, all_ids);
        auto affinity_out = binary_file(o.output_dir / "affinity.f64");
        for (double x : affinity) write_f64(affinity_out, x);
        finish_output(affinity_out, o.output_dir / "affinity.f64");

        const pp::AffinityTable hops(raw_dem);
        auto hops_out = binary_file(o.output_dir / "hops.u8");
        for (size_t i = 0; i <= n; ++i)
            for (size_t j = 0; j <= n; ++j)
                hops_out.put(static_cast<char>(static_cast<uint8_t>(hops(i, j))));
        finish_output(hops_out, o.output_dir / "hops.u8");

        auto syndromes_out = binary_file(o.output_dir / "syndromes.u8");
        auto labels_out = binary_file(o.output_dir / "labels.f64");
        dec::PyMatching ordinary(matching_dem);
        dec::PyMatching complementary(matching_dem, {.enable_gap_estimation = true});
        dec::PyMatching cluster(matching_dem, {.enable_cluster_gap = true});
        constexpr size_t batch_size = 512;
        stim::DemSampler<64> sampler(raw_dem, std::mt19937_64(o.seed), batch_size);
        uint64_t completed = 0;
        while (completed < o.shots) {
            sampler.resample(false);
            const auto detectors = sampler.det_buffer.transposed();
            const auto observables = sampler.obs_buffer.transposed();
            const size_t rows = std::min<uint64_t>(batch_size, o.shots - completed);
            for (size_t row = 0; row < rows; ++row) {
                for (size_t det = 0; det < n; ++det) syndromes_out.put(detectors[row][det] ? '\1' : '\0');
                const auto baseline = ordinary.decode(detectors[row], observables[row]);
                const auto comp = complementary.decode(detectors[row], observables[row]);
                const auto cl = cluster.decode(detectors[row], observables[row]);
                write_f64(labels_out, std::abs(comp.matching_data.gap));
                write_f64(labels_out, cl.cluster_gap_db);
                write_f64(labels_out, observables[row][0] ? 1.0 : 0.0);
                write_f64(labels_out, baseline.flipped_obs[0] ? 1.0 : 0.0);
                write_f64(labels_out, static_cast<double>(baseline.decode_latency_ns));
                write_f64(labels_out, static_cast<double>(cl.decode_latency_ns));
                ++completed;
            }
            if (completed % 102400 == 0 || completed == o.shots)
                std::cerr << "exported " << completed << " / " << o.shots << " shots\n";
        }
        finish_output(syndromes_out, o.output_dir / "syndromes.u8");
        finish_output(labels_out, o.output_dir / "labels.f64");
        write_metadata(o, n, o.output_dir / "metadata.json");
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "export_affinity_dataset: " << e.what() << '\n';
        return 1;
    }
}
