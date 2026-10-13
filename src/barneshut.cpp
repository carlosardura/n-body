#include "barneshut.h"
#include <cmath>
#include <algorithm>

BarnesHut::BarnesHut() {
    nodes.resize(2500000); 
    octants.resize(2500000); 
    node_count = 0;
}


int32_t BarnesHut::construct_cell(double xmin, double ymin, double zmin, double size, int32_t depth) {
    int32_t  idx = node_count++;

    MultipoleNode& node = nodes[idx];
    OctantDomain& octant = octants[idx];

    octant.xmin = xmin;
    octant.ymin = ymin;
    octant.zmin = zmin;
    octant.mid_x = xmin + size / 2.0;
    octant.mid_y = ymin + size / 2.0;
    octant.mid_z = zmin + size / 2.0;
    octant.size = size;
    octant.max_span = size;

    node.depth = depth;
    node.first_child = -1;
    node.body_idx = -1;
    node.M = 0.0;
    
    return idx;
}


int32_t BarnesHut::get_octant(int32_t idx, const PosMass& p) const {
    const OctantDomain& octant = octants[idx];
    return (p.x >= octant.mid_x) | ((p.y >= octant.mid_y) << 1) | ((p.z >= octant.mid_z) << 2);
}


inline void BarnesHut::new_octant(int32_t cell_idx) {
    if (nodes[cell_idx].first_child != -1) return;

    nodes[cell_idx].first_child = node_count;
    node_count += 8;

    const OctantDomain& parent = octants[cell_idx];
    double mid_x = parent.mid_x;
    double mid_y = parent.mid_y;
    double mid_z = parent.mid_z;
    double half_size = parent.size / 2.0;
    int32_t next_depth = nodes[cell_idx].depth + 1;

    for (int i = 0; i < 8; ++i) {
        int32_t child_idx = nodes[cell_idx].first_child + i;
        
        MultipoleNode& child_node = nodes[child_idx];
        OctantDomain& child_octant = octants[child_idx];

        child_node.depth = next_depth;
        child_node.first_child = -1;
        child_node.body_idx = -1;
        child_node.M = 0.0;

        double nx = (i & 1) ? mid_x : parent.xmin;
        double ny = (i & 2) ? mid_y : parent.ymin;
        double nz = (i & 4) ? mid_z : parent.zmin;

        child_octant.xmin = nx;
        child_octant.ymin = ny;
        child_octant.zmin = nz;
        child_octant.mid_x = nx + half_size / 2.0;
        child_octant.mid_y = ny + half_size / 2.0;
        child_octant.mid_z = nz + half_size / 2.0;
        child_octant.size = half_size;
        child_octant.max_span = half_size;
    }
}


void BarnesHut::insert_body(int32_t cell_idx, int32_t body_idx, const std::vector<PosMass>& pos_mass) {
    const PosMass& p = pos_mass[body_idx];
    MultipoleNode& node = nodes[cell_idx];

    if (node.M > 0) { 
        if (node.depth < 20) { 
            if (node.first_child == -1) {
                new_octant(cell_idx);
                
                int32_t old_body_idx = node.body_idx;
                int32_t old_octant = get_octant(cell_idx, pos_mass[old_body_idx]);
                insert_body(node.first_child + old_octant, old_body_idx, pos_mass);
                
                node.body_idx = -1;
            }

            int32_t octant = get_octant(cell_idx, p);
            insert_body(node.first_child + octant, body_idx, pos_mass);
        }

        // updates the center of mass
        double m_total = node.M + p.mass;
        node.x_cm = (node.M * node.x_cm + p.mass * p.x) / m_total;
        node.y_cm = (node.M * node.y_cm + p.mass * p.y) / m_total;
        node.z_cm = (node.M * node.z_cm + p.mass * p.z) / m_total;
        node.M += p.mass;

    } else {   // empty nodes
        node.body_idx = body_idx;
        node.x_cm = p.x;
        node.y_cm = p.y;
        node.z_cm = p.z;
        node.M = p.mass;
    }
}


void BarnesHut::partition_space(const std::vector<PosMass>& pos_mass, const OctantDomain& bbox) {
    node_count = 0; 
    
    if (pos_mass.empty()) return;

    int32_t root = construct_cell(bbox.xmin, bbox.ymin, bbox.zmin, bbox.max_span, 0);

    for (size_t i = 0; i < pos_mass.size(); ++i) {
        insert_body(root, static_cast<int32_t>(i), pos_mass);
    }
}


void BarnesHut::compute_force(const std::vector<PosMass>& pos_mass, std::vector<PosMass>& accelerations, double theta, double epsilon, double G) {
    if (nodes.empty()) return;

    for (size_t i = 0; i < pos_mass.size(); ++i) {
        double ax = 0.0, ay = 0.0, az = 0.0;

        int32_t stack[512]; 
        int stack_ptr = 0;
        stack[stack_ptr++] = 0;

        while (stack_ptr > 0) {
            int32_t current_idx = stack[--stack_ptr];
            const MultipoleNode& node = nodes[current_idx];
            const OctantDomain& octant = octants[current_idx];

            if (node.M == 0.0 || node.body_idx == static_cast<int32_t>(i)) {
                continue;
            }

            double dx = node.x_cm - pos_mass[i].x;
            double dy = node.y_cm - pos_mass[i].y;
            double dz = node.z_cm - pos_mass[i].z;
            double d2 = dx*dx + dy*dy + dz*dz + epsilon*epsilon;

            if (((octant.size * octant.size) < (theta * theta * d2)) || (node.first_child == -1)) {
                // multipole acceptance criterion
                double inv_d = 1.0 / std::sqrt(d2);
                double a = G * node.M * inv_d * inv_d * inv_d; 

                ax += a * dx; 
                ay += a * dy; 
                az += a * dz; 
            } else {
                for (int j = 0; j < 8; ++j) {
                    stack[stack_ptr++] = node.first_child + j;
                }
            }
        }
        accelerations[i].x = ax;
        accelerations[i].y = ay;
        accelerations[i].z = az;
    }
}