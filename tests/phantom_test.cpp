#include "decoder/phantom.h"

#include <stdexcept>
#include <string>

namespace
{

void check(const stim::DetectorErrorModel& dem, size_t distance,
           const std::vector<size_t>& fired, bool expected)
{
    dec::Phantom decoder(dem, distance);
    dec::SyndromeType syndrome(dem.count_detectors());
    dec::ObsType observations(dem.count_observables());
    for (size_t detector : fired)
        syndrome[detector] = 1;
    const auto result = decoder.decode(syndrome, observations);
    if (bool(result.flipped_obs[0]) != expected)
        throw std::runtime_error("Phantom returned the wrong correction");
}

} // namespace

int main()
{
    // A fired detector beyond the first popcnt() positions must be collected.
    const stim::DetectorErrorModel single("error(0.2) D2 L0\n");
    check(single, 5, {2}, true);

    // An even-sized pair connected by one edge must be decoded together.
    const stim::DetectorErrorModel pair(
        "error(0.2) D0 D2 L0\n"
        "error(0.001) D0\n"
        "error(0.001) D2\n");
    check(pair, 5, {0, 2}, true);

    // Transitive edges must merge at their roots instead of overwriting the
    // connection made by an earlier edge.
    const stim::DetectorErrorModel transitive(
        "error(0.2) D0 D2\n"
        "error(0.3) D1 D2 L0\n"
        "error(0.1) D0\n"
        "error(0.001) D1\n"
        "error(0.001) D2\n");
    check(transitive, 5, {0, 1, 2}, true);

    // A connected cluster above the hardware limit receives no correction.
    std::string chain = "error(0.1) D0 L0\n";
    for (size_t i = 0; i < 10; i++)
        chain += "error(0.1) D" + std::to_string(i) + " D" + std::to_string(i + 1) + "\n";
    const stim::DetectorErrorModel large(chain);
    check(large, 5, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, false);

    // An isolated small cluster is still decoded beside a large component.
    std::string wide_chain = "error(0.1) D0 L0\n";
    std::vector<size_t> wide_fired;
    for (size_t i = 0; i < 65; i++)
    {
        wide_fired.push_back(i);
        if (i < 64)
            wide_chain += "error(0.1) D" + std::to_string(i) + " D" + std::to_string(i + 1) + "\n";
    }
    wide_chain += "error(0.2) D65 L0\n";
    wide_fired.push_back(65);
    check(stim::DetectorErrorModel(wide_chain), 5, wide_fired, true);

    bool rejected = false;
    try { dec::Phantom invalid(single, 0); }
    catch (const std::invalid_argument&) { rejected = true; }
    if (!rejected)
        throw std::runtime_error("Phantom accepted zero code distance");
}
