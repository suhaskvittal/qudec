/* Author: OpenAI Codex. Geometry is read from the existing circuit generator;
 * no independent tiling or inferred circuit-fault/data-qubit correspondence.
 */
#include "error_chain_display.h"
#include "circuit_generator.h"
#include <stim/util_top/circuit_to_dem.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <ostream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace
{
std::map<uint64_t,std::vector<double>> detector_coordinates(const stim::DetectorErrorModel& dem)
{
    std::set<uint64_t> ids;
    for (size_t id = 0; id < dem.count_detectors(); id++) ids.insert(id);
    return dem.get_detector_coordinates(ids);
}
}

ErrorChainDisplay::ErrorChainDisplay(const stim::DetectorErrorModel& dem,
    const std::string& code,uint32_t distance)
    :code_(code), coordinates_(detector_coordinates(dem))
{
    if (code.empty()) return; // Generic toy DEMs need no lattice metadata.
    if (code != "cc666" && code != "cc488") throw std::invalid_argument("ASCII display requires a color code");
    const CodeCapacityNoisePolicy noise(0);
    const auto reference = code == "cc666" ? cc666_circuit(distance,1,noise,false)
                                           : cc488_circuit(distance,1,noise,false);
    qubits_ = reference.get_final_qubit_coords();
    // One Z-product measurement per spatial check, in generator check-ID order.
    for (const auto& op : reference.operations) if (op.gate_type == stim::GateType::MPP)
    {
        // Stim may fuse adjacent MPP instructions into a single operation.
        // Split its products at adjacent Pauli targets without a combiner.
        std::vector<uint32_t> product;
        bool z_product = false, joined = false;
        auto finish_product = [&]()
        {
            if (z_product && !product.empty())
            {
                Check check{.qubits=product};
                for (auto q : check.qubits)
                {
                    check.x += qubits_.at(q)[0]/check.qubits.size();
                    check.y += qubits_.at(q)[1]/check.qubits.size();
                }
                checks_.push_back(std::move(check));
            }
            product.clear();
        };
        for (const auto& target : op.targets)
        {
            if (target.is_combiner()) {joined = true; continue;}
            if (!joined && !product.empty()) finish_product();
            if (product.empty()) z_product = target.is_z_target();
            product.push_back(target.qubit_value());
            joined = false;
        }
        finish_product();
    }
    if (checks_.empty()) throw std::runtime_error("No color-code check supports for ASCII display");
    const auto reference_coordinates = detector_coordinates(stim::circuit_to_dem(reference,{}));
    for (size_t id = 0; id < checks_.size(); id++)
    {
        const auto& c = reference_coordinates.at(id);
        checks_[id].color = static_cast<int>(c[3]);
        if (std::abs(checks_[id].x-c[0])>1e-8 || std::abs(checks_[id].y-c[1])>1e-8)
            throw std::runtime_error("Check geometry and detector IDs disagree");
    }
    for (const auto& [detector,c] : coordinates_)
    {
        if (c.size() != 6) throw std::invalid_argument("Missing color-code detector coordinates");
        auto match = checks_.size();
        for (size_t id = 0; id < checks_.size(); id++)
            if (std::abs(checks_[id].x-c[0])<1e-8 && std::abs(checks_[id].y-c[1])<1e-8)
            {
                if (match != checks_.size()) throw std::runtime_error("Ambiguous spatial check coordinates");
                match = id;
            }
        if (match == checks_.size()) throw std::runtime_error("Detector does not map to a spatial check");
        detector_checks_[detector] = match;
    }
}

std::string ErrorChainDisplay::detector_label(size_t id) const
{
    std::ostringstream out;
    out << 'D' << id;
    const auto found = detector_checks_.find(id);
    if (found == detector_checks_.end()) return out.str();
    const auto check = found->second;
    const auto& c = coordinates_.at(id);
    out << '=' << (c[5] == 1 ? 'F' : "RGB"[checks_[check].color])
        << '(' << check << ", " << static_cast<int64_t>(std::llround(c[2])) << "):"
        << (c[4] == 0 ? 'X' : 'Z');
    return out.str();
}

void ErrorChainDisplay::draw(dec::SyndromeRef syndrome,std::ostream& out) const
{
    if (checks_.empty()) return;
    std::vector<bool> fired(checks_.size()), flagged(checks_.size());
    for (const auto& [detector,check] : detector_checks_) if (syndrome[detector])
    {
        if (coordinates_.at(detector)[5] == 1) flagged[check] = true;
        else fired[check] = true;
    }
    std::vector<std::string> labels;
    for (size_t id = 0; id < checks_.size(); id++)
    {
        std::string label = std::string(1,"RGB"[checks_[id].color])+std::to_string(id);
        if (fired[id]) label = '['+label+']';
        if (flagged[id]) label += '!';
        labels.push_back(std::move(label));
    }
    // Scale until labels fit without overwriting lattice edges or data vertices.
    std::vector<std::string> canvas;
    for (int scale = 3; scale <= 24; scale++)
    {
        const double sx = 2*scale, sy = code_ == "cc666" ? 2*scale/std::sqrt(3.0) : scale;
        const int margin = 12;
        auto point = [&](double x,double y)
        {return std::pair<int,int>{margin+static_cast<int>(std::lround(sx*x)),margin+static_cast<int>(std::lround(sy*y))};};
        int width = 0,height = 0;
        for (const auto& [q,c] : qubits_)
        {
            auto [x,y] = point(c[0],c[1]);
            width = std::max(width,x+margin); height = std::max(height,y+margin);
        }
        canvas.assign(height+1,std::string(width+1,' '));
        for (const auto& check : checks_)
        {
            auto vertices = check.qubits;
            std::sort(vertices.begin(),vertices.end(),[&](auto a,auto b)
            {
                const auto &ca = qubits_.at(a), &cb = qubits_.at(b);
                return std::atan2(ca[1]-check.y,ca[0]-check.x) < std::atan2(cb[1]-check.y,cb[0]-check.x);
            });
            for (size_t j = 0; j < vertices.size(); j++)
            {
                const auto &a = qubits_.at(vertices[j]), &b = qubits_.at(vertices[(j+1)%vertices.size()]);
                const auto [x,y] = point(a[0],a[1]);
                const auto [xx,yy] = point(b[0],b[1]);
                const int dx = xx-x,dy = yy-y, steps = std::max(std::abs(dx),std::abs(dy));
                const char edge = dy == 0 ? '-' : dx == 0 ? '|' : dx*dy>0 ? '\\' : '/';
                for (int k = 1; k < steps; k++)
                    canvas[y+std::lround(static_cast<double>(dy)*k/steps)][x+std::lround(static_cast<double>(dx)*k/steps)] = edge;
            }
        }
        for (const auto& [q,c] : qubits_)
        {
            auto [x,y] = point(c[0],c[1]); canvas[y][x] = '.';
        }
        bool collision = false;
        for (size_t id = 0; id < checks_.size(); id++)
        {
            auto [x,y] = point(checks_[id].x,checks_[id].y);
            x -= labels[id].size()/2;
            for (size_t j = 0; j < labels[id].size(); j++)
            {
                if (canvas.at(y).at(x+j) != ' ') collision = true;
                canvas[y][x+j] = labels[id][j];
            }
        }
        if (!collision) break;
        if (scale == 24) throw std::runtime_error("Could not fit ASCII check labels");
    }
    out << "Color-code lattice (" << code_ << "), combined syndrome across all rounds (checks XOR, flags OR):\n"
        << "  . = data qubit; [R/G/B id] = check fired in at least one round; ! = colocated colorless flag fired\n"
        << "  Unbracketed checks are quiet. Data-fault locations are not inferred from DEM events.\n";
    // Trim the extra canvas margin, preserving the relative position of every edge.
    size_t left = std::numeric_limits<size_t>::max();
    for (const auto& row : canvas)
    {
        const auto first = row.find_first_not_of(' ');
        if (first != std::string::npos) left = std::min(left,first);
    }
    for (const auto& row : canvas)
    {
        const auto last = row.find_last_not_of(' ');
        if (last != std::string::npos) out << "  " << row.substr(left,last-left+1) << '\n';
    }
    out << "Fired detector identities: (check-id, round), zero-based rounds; final closure is round=extraction-rounds.\n";
    bool any = false;
    for (const auto& [detector,check] : detector_checks_) if (syndrome[detector])
    {
        any = true;
        out << "  " << detector_label(detector) << (coordinates_.at(detector)[5] == 1 ? " (colorless flag)" : " (check)") << '\n';
    }
    if (!any) out << "  none\n";
}
