#pragma once
#include "types.h"

namespace MortonCode {
    void sort_system(BodySystem& system, BodySystem& morton_system, const OctantDomain& bbox);
}