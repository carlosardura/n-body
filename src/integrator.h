#pragma once
#include "types.h"
#include <vector>

namespace Integrator {
    void velverlet_step1(BodySystem& system, const std::vector<PosMass>& accelerations, double dt);
    void velverlet_step2(BodySystem& system, const std::vector<PosMass>& accelerations, double dt);
}