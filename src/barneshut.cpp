#include "barneshut.h"
#include <cmath>
#include <algorithm>

BarnesHut::BarnesHut() {
    nodes.resize(2500000); 
    octants.resize(2500000); 
    node_count = 0;
}


int32_t BarnesHut::construct_cell(double xmin, double ymin, double zmin, double size, int32_t depth) {
    int32_t idx = node_count++;

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


void BarnesHut::insert_body(int32_t cell_idx, size_t start, size_t end, const BodySystem& system) {
    if (start >= end) return;

    MultipoleNode& node = nodes[cell_idx];
    const OctantDomain& octant = octants[cell_idx];

    if (end - start == 1) {
        size_t body_idx = start;
        node.body_idx = static_cast<int32_t>(body_idx);
        node.x_cm = system.x[body_idx];
        node.y_cm = system.y[body_idx];
        node.z_cm = system.z[body_idx];
        node.M = system.mass[body_idx];
        return;
    }

    double total_m = 0.0, x_cm = 0.0, y_cm = 0.0, z_cm = 0.0;
    for (size_t i = start; i < end; ++i) {
        double m = system.mass[i];
        total_m += m;
        x_cm += m * system.x[i];
        y_cm += m * system.y[i];
        z_cm += m * system.z[i];
    }
    node.M = total_m;
    node.x_cm = x_cm / total_m;
    node.y_cm = y_cm / total_m;
    node.z_cm = z_cm / total_m;

    if (node.depth >= 20) {
        return; 
    }

    new_octant(cell_idx);
    double mid_x = octant.mid_x;
    double mid_y = octant.mid_y;
    double mid_z = octant.mid_z;

    size_t octant_starts[9];
    octant_starts[0] = start;
    octant_starts[8] = end;

    size_t current = start;
    for (int oct = 0; oct < 7; ++oct) {
        size_t count = end - current;
        size_t first = current;
        while (count > 0) {
            size_t step = count / 2;
            size_t it = first + step;
            int o = (system.x[it] >= mid_x) | ((system.y[it] >= mid_y) << 1) | ((system.z[it] >= mid_z) << 2);
            if (o <= oct) {
                first = it + 1;
                count -= step + 1;
            } else {
                count = step;
            }
        }
        current = first;
        octant_starts[oct + 1] = current;
    }

    for (int i = 0; i < 8; ++i) {
        size_t o_start = octant_starts[i];
        size_t o_end = octant_starts[i + 1];
        if (o_start < o_end) {
            insert_body(node.first_child + i, o_start, o_end, system);
        }
    }
}


void BarnesHut::partition_space(const BodySystem& system, const OctantDomain& bbox) {
    node_count = 0; 
    if (system.size() == 0) return;

    int32_t root_idx = construct_cell(bbox.xmin, bbox.ymin, bbox.zmin, bbox.max_span, 0);
    insert_body(root_idx, 0, system.size(), system);
}


void BarnesHut::compute_force(const BodySystem& system, std::vector<double>& ax_out, std::vector<double>& ay_out, std::vector<double>& az_out, double theta, double epsilon, double G) {
    if (nodes.empty()) return;
    size_t n = system.size();

    for (size_t i = 0; i < n; i += 4) {
        int block_size = (n - i < 4) ? (n - i) : 4;
        double ax[4] = {0.0}, ay[4] = {0.0}, az[4] = {0.0};
        double px[4], py[4], pz[4];
        
        for (int k = 0; k < block_size; ++k) {
            px[k] = system.x[i + k];
            py[k] = system.y[i + k];
            pz[k] = system.z[i + k];
        }

        int32_t stack[512]; 
        int stack_ptr = 0;
        stack[stack_ptr++] = 0;

        while (stack_ptr > 0) {
            int32_t current_idx = stack[--stack_ptr];
            const MultipoleNode& node = nodes[current_idx];
            const OctantDomain& octant = octants[current_idx];

            if (node.M == 0.0) continue;

            bool self_interaction = false;
            for(int k = 0; k < block_size; ++k) {
                if (node.body_idx == static_cast<int32_t>(i + k)) self_interaction = true;
            }
            if (self_interaction) continue;
            
            double d2[4];
            bool open_node = false;

            for (int k = 0; k < block_size; ++k) {
                double dx = node.x_cm - px[k];
                double dy = node.y_cm - py[k];
                double dz = node.z_cm - pz[k];
                d2[k] = dx*dx + dy*dy + dz*dz + epsilon*epsilon;
                
                if ((octant.size * octant.size) >= (theta * theta * d2[k])) {
                    // multipole acceptance criterion
                    open_node = true; 
                }
            }

            if (open_node && node.first_child != -1) {
                for (int j = 0; j < 8; ++j) {
                    stack[stack_ptr++] = node.first_child + j;
                }
            } else {
                for (int k = 0; k < block_size; ++k) {
                    double dx = node.x_cm - px[k];
                    double dy = node.y_cm - py[k];
                    double dz = node.z_cm - pz[k];
                    
                    double inv_d = 1.0 / std::sqrt(d2[k]);
                    double a = G * node.M * inv_d * inv_d * inv_d;
                    
                    ax[k] += a * dx;
                    ay[k] += a * dy;
                    az[k] += a * dz;
                }
            }
        }

        for (int k = 0; k < block_size; ++k) {
            ax_out[i + k] = ax[k];
            ay_out[i + k] = ay[k];
            az_out[i + k] = az[k];
        }
    }
}