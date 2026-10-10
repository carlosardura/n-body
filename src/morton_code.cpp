#include "morton_code.h"
#include <libmorton/morton.h>
#include <vector>
#include <algorithm>
#include <cmath>

namespace MortonCode {

    struct ParticleCode {
        uint64_t code;
        uint32_t original_index;
    };

    void sort_system(BodySystem& system, const Cell& bbox) {
        size_t n = system.size();
        if (n <= 1 || bbox.max_span == 0.0) return;

        double inv_span = 1.0 / bbox.max_span;
        const uint32_t max_bins = (1 << 21) - 1;

        std::vector<ParticleCode> codes(n);
        for (size_t i = 0; i < n; ++i) {
            uint32_t bx = static_cast<uint32_t>((system.pos_mass[i].x - bbox.xmin) * inv_span * max_bins);
            uint32_t by = static_cast<uint32_t>((system.pos_mass[i].y - bbox.ymin) * inv_span * max_bins);
            uint32_t bz = static_cast<uint32_t>((system.pos_mass[i].z - bbox.zmin) * inv_span * max_bins);
            
            codes[i].code = libmorton::morton3D_64_encode(bx, by, bz);
            codes[i].original_index = static_cast<uint32_t>(i);
        }

        std::sort(codes.begin(), codes.end(), [](const ParticleCode& a, const ParticleCode& b) {
            return a.code < b.code;
        });

        std::vector<PosMass> temp_pos_mass(n);
        std::vector<Velocity> temp_velocities(n);
        std::vector<uint32_t> temp_ids(n);

        for (size_t i = 0; i < n; ++i) {
            uint32_t orig_idx = codes[i].original_index;
            temp_pos_mass[i] = system.pos_mass[orig_idx];
            temp_velocities[i] = system.velocities[orig_idx];
            temp_ids[i] = system.ids[orig_idx];
        }

        system.pos_mass = std::move(temp_pos_mass);
        system.velocities = std::move(temp_velocities);
        system.ids = std::move(temp_ids);
    }
}