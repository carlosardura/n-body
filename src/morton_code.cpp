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

    void sort_system(BodySystem& system, BodySystem& morton_system, const OctantDomain& bbox) {
        size_t n = system.size();
        if (n <= 1 || bbox.max_span == 0.0) return;

        double inv_span = 1.0 / bbox.max_span;
        const uint32_t max_bins = (1 << 21) - 1;

        std::vector<ParticleCode> codes(n);
        for (size_t i = 0; i < n; ++i) {
            uint32_t bx = static_cast<uint32_t>((system.x[i] - bbox.xmin) * inv_span * max_bins);
            uint32_t by = static_cast<uint32_t>((system.y[i] - bbox.ymin) * inv_span * max_bins);
            uint32_t bz = static_cast<uint32_t>((system.z[i] - bbox.zmin) * inv_span * max_bins);
            
            codes[i].code = libmorton::morton3D_64_encode(bx, by, bz);
            codes[i].original_index = static_cast<uint32_t>(i);
        }

        std::sort(codes.begin(), codes.end(), [](const ParticleCode& a, const ParticleCode& b) {
            return a.code < b.code;
        });

        for (size_t i = 0; i < n; ++i) {
            uint32_t orig_idx = codes[i].original_index;
            morton_system.x[i] = system.x[orig_idx];
            morton_system.y[i] = system.y[orig_idx];
            morton_system.z[i] = system.z[orig_idx];
            morton_system.mass[i] = system.mass[orig_idx];

            morton_system.vx[i] = system.vx[orig_idx];
            morton_system.vy[i] = system.vy[orig_idx];
            morton_system.vz[i] = system.vz[orig_idx];

            morton_system.ids[i] = system.ids[orig_idx];
        }

        std::swap(system.x, morton_system.x);
        std::swap(system.y, morton_system.y);
        std::swap(system.z, morton_system.z);
        std::swap(system.mass, morton_system.mass);
        std::swap(system.vx, morton_system.vx);
        std::swap(system.vy, morton_system.vy);
        std::swap(system.vz, morton_system.vz);
        std::swap(system.ids, morton_system.ids);
    }
}