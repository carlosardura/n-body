#include "integrator.h"

namespace Integrator {
    void velverlet_step1(BodySystem& system, const std::vector<PosMass>& accelerations, double dt) {
        double dt_half = 0.5 * dt;
        
        for (size_t i = 0; i < system.size(); ++i) {
            system.velocities[i].vx += accelerations[i].x * dt_half;
            system.velocities[i].vy += accelerations[i].y * dt_half;
            system.velocities[i].vz += accelerations[i].z * dt_half;

            system.pos_mass[i].x += system.velocities[i].vx * dt;
            system.pos_mass[i].y += system.velocities[i].vy * dt;
            system.pos_mass[i].z += system.velocities[i].vz * dt;
        }
    }

    void velverlet_step2(BodySystem& system, const std::vector<PosMass>& accelerations, double dt) {
        double dt_half = 0.5 * dt;
        
        for (size_t i = 0; i < system.size(); ++i) {
            system.velocities[i].vx += accelerations[i].x * dt_half;
            system.velocities[i].vy += accelerations[i].y * dt_half;
            system.velocities[i].vz += accelerations[i].z * dt_half;
        }
    }

}