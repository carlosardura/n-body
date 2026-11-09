#include "integrator.h"
#include <omp.h>

namespace Integrator {

    void velverlet_step1(BodySystem& system, const std::vector<double>& ax, const std::vector<double>& ay, const std::vector<double>& az, double dt) {
        double dt_half = 0.5 * dt;
        size_t n = system.size();
        
        #pragma omp parallel for simd schedule(static)
        for (size_t i = 0; i < n; ++i) {
            system.vx[i] += ax[i] * dt_half;
            system.vy[i] += ay[i] * dt_half;
            system.vz[i] += az[i] * dt_half;

            system.x[i] += system.vx[i] * dt;
            system.y[i] += system.vy[i] * dt;
            system.z[i] += system.vz[i] * dt;
        }
    }

    void velverlet_step2(BodySystem& system, const std::vector<double>& ax, const std::vector<double>& ay, const std::vector<double>& az, double dt) {
        double dt_half = 0.5 * dt;
        size_t n = system.size();

        #pragma omp parallel for simd schedule(dynamic, 64)
        for (size_t i = 0; i < n; ++i) {
            system.vx[i] += ax[i] * dt_half;
            system.vy[i] += ay[i] * dt_half;
            system.vz[i] += az[i] * dt_half;
        }
    }

}