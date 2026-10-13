#pragma once
#include "types.h"
#include <vector>

class BarnesHut {
public:
    BarnesHut();
    void partition_space(const std::vector<PosMass>& pos_mass, const OctantDomain& bbox);
    void compute_force(const std::vector<PosMass>& pos_mass, std::vector<PosMass>& accelerations, double theta, double epsilon, double G);

private:
    std::vector<MultipoleNode> nodes;
    std::vector<OctantDomain> octants;
    int32_t node_count;
    
    int32_t construct_cell(double xmin, double ymin, double zmin, double size, int32_t depth);
    int32_t get_octant(int32_t cell_idx, const PosMass& p) const;
    void new_octant(int32_t cell_idx);
    void insert_body(int32_t cell_idx, int32_t body_idx, const std::vector<PosMass>& pos_mass);
};