#ifndef QUDEC_CODE_CATALOG_H
#define QUDEC_CODE_CATALOG_H

#include <array>
#include <cstdint>
#include <string_view>

// Named non-planar codes. For a distance interval, default rounds use its
// proven lower bound; the interval is included in the code name.
struct CatalogCode {
    std::string_view name;
    std::string_view directory;
    uint32_t default_rounds;
};

inline constexpr std::array<CatalogCode, 15> external_codes{{
    {"hycc_888_64_8_6", "8.8.8_L2_D444-32_n64_k8", 6},
    {"hycc_888_128_12_8", "8.8.8_L2_D444-40_n128_k12", 8},
    {"hycc_888_512_36_12", "8.8.8_L2_D444-58_n512_k36", 12},
    {"hycc_888_1600_104_14to16", "8.8.8_L2_D444-109_n1600_k104", 14},
    {"hycc_888_1344_88_14", "8.8.8_L2_D444-96_n1344_k88", 14},
    {"hycc_888_5376_340_18to20", "8.8.8_L2_D444-271_n5376_k340", 18},
    {"hycc_888_24288_1522_22", "8.8.8_L2_D444-PSL2-23_n24288_k1522", 22},
    {"hycc_888_31200_1954_22to24", "8.8.8_L2_D444-PSL2-25_n31200_k1954", 22},
    {"hycc_688_384_20_12", "6.8.8_L2_D344-18_n384_k20", 12},
    {"hycc_688_768_36_14", "6.8.8_L2_D344-22_n768_k36", 14},
    {"hycc_688_1536_68_16", "6.8.8_L2_D344-32_n1536_k68", 16},
    {"hycc_688_3840_164_20", "6.8.8_L2_D344-45_n3840_k164", 20},
    {"hycc_688_10368_436_22", "6.8.8_L2_D344-79_n10368_k436", 22},
    {"hycc_688_19584_820_24", "6.8.8_L2_D344-PGL2-17_n19584_k820", 24},
    {"hycc_688_48576_2028_28", "6.8.8_L2_D344-PGL2-23_n48576_k2028", 28},
}};

inline const CatalogCode* find_external_code(std::string_view name) {
    for (const auto& code : external_codes)
        if (code.name == name) return &code;
    return nullptr;
}

#endif
