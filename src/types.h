#pragma once
#include <vector>
#include <cstdint>

struct alignas(32) PosMass {
    double x = 0.0, y = 0.0, z = 0.0;
    double mass = 0.0;
};

struct alignas(32) Velocity {
    double vx = 0.0, vy = 0.0, vz = 0.0;
    double padding = 0.0; 
};

struct Cell {
    double xmin = 0.0, ymin = 0.0, zmin = 0.0;
    double mid_x = 0.0, mid_y = 0.0, mid_z = 0.0;
    double size = 0.0; 

    double x_cm = 0.0, y_cm = 0.0, z_cm = 0.0; 
    double M = 0.0;

    int32_t children[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
    int32_t body_idx = -1;
    bool external = true;
    int32_t depth = 0;
};

struct BodySystem {
    std::vector<PosMass> pos_mass;
    std::vector<Velocity> velocities;
    std::vector<uint32_t> ids; 
 
    size_t size() const { return pos_mass.size(); }
    
    void reserve(size_t n) {
        pos_mass.reserve(n);
        velocities.reserve(n);
        ids.reserve(n);
    }
};