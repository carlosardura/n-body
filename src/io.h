#pragma once
#include "types.h"
#include <string>

namespace IO {
    void dump_frame(const std::string& directory, const BodySystem& system, int step, int total_steps);
}