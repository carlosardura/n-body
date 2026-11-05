#pragma once
#include "types.h"
#include <vector>

namespace Integrator {
    void velverlet_step1(BodySystem& system, const std::vector<double>& ax, const std::vector<double>& ay, const std::vector<double>& az, double dt);
    void velverlet_step2(BodySystem& system, const std::vector<double>& ax, const std::vector<double>& ay, const std::vector<double>& az, double dt);
}