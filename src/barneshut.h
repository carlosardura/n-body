#pragma once
#include "types.h"
#include <vector>

class BarnesHut {
public:
    BarnesHut();
    void partition_space(const BodySystem& system, const OctantDomain& bbox, double theta);
    void compute_force(const BodySystem& system, std::vector<double>& ax_out, std::vector<double>& ay_out, std::vector<double>& az_out, double epsilon, double G);

    private:
    std::vector<MultipoleNode> nodes;
    std::vector<OctantDomain> octants;
    int32_t node_count;
    double inv_theta2 = 0.0;
    
    int32_t construct_cell(double xmin, double ymin, double zmin, double size, int32_t depth);
    void new_octant(int32_t cell_idx);
    void insert_body(int32_t cell_idx, size_t start, size_t end, const BodySystem& system);
    void validation_tree(int32_t node_idx, int32_t next_ptr);
};