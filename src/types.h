#pragma once
#include <vector>
#include <cstdint>

struct alignas(64) OctantDomain {
    double xmin = 0.0, ymin = 0.0, zmin = 0.0;
    double mid_x = 0.0, mid_y = 0.0, mid_z = 0.0;
    double size = 0.0; 
    double max_span = 0.0;
};

struct alignas(64) MultipoleNode {
    double x_cm = 0.0, y_cm = 0.0, z_cm = 0.0;
    double M = 0.0;
    int32_t first_child = -1;
    int32_t body_idx = -1;
    int32_t depth = 0;
    int32_t padding_node[5] = {0, 0, 0, 0, 0};
};
    

struct BodySystem {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<double> mass;

    std::vector<double> vx;
    std::vector<double> vy;
    std::vector<double> vz;

    std::vector<uint32_t> ids;

    size_t num_bodies = 0;

    size_t size() const { return num_bodies; }

    void reserve(size_t n) {
        x.reserve(n); y.reserve(n); z.reserve(n); mass.reserve(n);
        vx.reserve(n); vy.reserve(n); vz.reserve(n);
        ids.reserve(n);
    }
    
    void resize(size_t n) {
        x.resize(n, 0.0); y.resize(n, 0.0); z.resize(n, 0.0); mass.resize(n, 0.0);
        vx.resize(n, 0.0); vy.resize(n, 0.0); vz.resize(n, 0.0);
        ids.resize(n, 0);
        num_bodies = n;
    }
};