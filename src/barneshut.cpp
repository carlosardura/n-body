#include "barneshut.h"
#include <cmath>
#include <algorithm>

BarnesHut::BarnesHut() {
    cells.reserve(2000000); 
}


int32_t BarnesHut::cell(double xmin, double ymin, double zmin, double size, int32_t depth) {
    Cell cell;
    cell.xmin = xmin; 
    cell.ymin = ymin; 
    cell.zmin = zmin;

    cell.mid_x = xmin + size / 2.0; 
    cell.mid_y = ymin + size / 2.0; 
    cell.mid_z = zmin + size / 2.0;

    cell.size = size; 
    cell.depth = depth;
    
    cells.push_back(cell);
    return static_cast<int32_t>(cells.size() - 1);
}


int32_t BarnesHut::get_octant(int32_t cell_idx, const PosMass& p) const {
    const Cell& cell = cells[cell_idx];
    return (p.x >= cell.mid_x) | ((p.y >= cell.mid_y) << 1) | ((p.z >= cell.mid_z) << 2); 
}


inline void BarnesHut::new_octant(int32_t cell_idx, int32_t octant) {
    if (cells[cell_idx].children[octant] == -1) {
        double nx = (octant & 1) ? cells[cell_idx].mid_x : cells[cell_idx].xmin;
        double ny = (octant & 2) ? cells[cell_idx].mid_y : cells[cell_idx].ymin;
        double nz = (octant & 4) ? cells[cell_idx].mid_z : cells[cell_idx].zmin;
        cells[cell_idx].children[octant] = cell(nx, ny, nz, cells[cell_idx].size / 2.0, cells[cell_idx].depth + 1);
    }
}


void BarnesHut::insert_body(int32_t cell_idx, int32_t body_idx, const std::vector<PosMass>& pos_mass) {
    const PosMass& p = pos_mass[body_idx];

    if (cells[cell_idx].M > 0) { 
        if (cells[cell_idx].depth < 20) { 
            if (cells[cell_idx].external) {
                cells[cell_idx].external = false; 
                
                int32_t old_body_idx = cells[cell_idx].body_idx;
                int32_t old_octant = get_octant(cell_idx, pos_mass[old_body_idx]);
                new_octant(cell_idx, old_octant);
                insert_body(cells[cell_idx].children[old_octant], old_body_idx, pos_mass);
                
                cells[cell_idx].body_idx = -1;
            }

            int32_t octant = get_octant(cell_idx, p);
            new_octant(cell_idx, octant);
            insert_body(cells[cell_idx].children[octant], body_idx, pos_mass);
        }

        // updates the center of mass
        double m_total = cells[cell_idx].M + p.mass;
        cells[cell_idx].x_cm = (cells[cell_idx].M * cells[cell_idx].x_cm + p.mass * p.x) / m_total;
        cells[cell_idx].y_cm = (cells[cell_idx].M * cells[cell_idx].y_cm + p.mass * p.y) / m_total;
        cells[cell_idx].z_cm = (cells[cell_idx].M * cells[cell_idx].z_cm + p.mass * p.z) / m_total;
        cells[cell_idx].M += p.mass;

    } else {   // empty nodes
        cells[cell_idx].body_idx = body_idx;
        cells[cell_idx].x_cm = p.x;
        cells[cell_idx].y_cm = p.y;
        cells[cell_idx].z_cm = p.z;
        cells[cell_idx].M = p.mass;
        cells[cell_idx].external = true;
    }
}


void BarnesHut::partition_space(const std::vector<PosMass>& pos_mass) {
    cells.clear();
    if (pos_mass.empty()) return;

    double min_x = pos_mass[0].x, max_x = pos_mass[0].x;
    double min_y = pos_mass[0].y, max_y = pos_mass[0].y;
    double min_z = pos_mass[0].z, max_z = pos_mass[0].z;

    for (const auto& p : pos_mass) {
        if (p.x < min_x) min_x = p.x; if (p.x > max_x) max_x = p.x;
        if (p.y < min_y) min_y = p.y; if (p.y > max_y) max_y = p.y;
        if (p.z < min_z) min_z = p.z; if (p.z > max_z) max_z = p.z;
    }

    double max_span = std::max({max_x - min_x, max_y - min_y, max_z - min_z});
    int32_t root = cell(min_x, min_y, min_z, max_span, 0);

    for (size_t i = 0; i < pos_mass.size(); ++i) {
        insert_body(root, static_cast<int32_t>(i), pos_mass);
    }
}


void BarnesHut::update_accelerations(const std::vector<PosMass>& pos_mass, std::vector<PosMass>& accelerations, double theta, double epsilon, double G) {
    if (cells.empty()) return;
    
    std::vector<int32_t> stack;
    stack.reserve(256);

    for (size_t i = 0; i < pos_mass.size(); ++i) {
        double ax = 0.0, ay = 0.0, az = 0.0;
        stack.clear();
        stack.push_back(0);

        while (!stack.empty()) {
            int32_t current_idx = stack.back();
            stack.pop_back(); 

            const Cell& current = cells[current_idx];

            if (current.M == 0.0 || current.body_idx == static_cast<int32_t>(i)) {
                continue;
            }

            double dx = current.x_cm - pos_mass[i].x;
            double dy = current.y_cm - pos_mass[i].y;
            double dz = current.z_cm - pos_mass[i].z;
            double d2 = dx*dx + dy*dy + dz*dz + epsilon*epsilon;

            if (((current.size * current.size) < (theta * theta * d2)) || current.external) {
                // multipole acceptance criterion
                double a = G * current.M / (d2 * std::sqrt(d2)); 
                ax += a * dx; 
                ay += a * dy; 
                az += a * dz; 
            } else {
                for (int j = 0; j < 8; ++j) {
                    if (current.children[j] != -1) {
                        stack.push_back(current.children[j]);
                    }
                }
            }
        }
        accelerations[i].x = ax;
        accelerations[i].y = ay;
        accelerations[i].z = az;
    }
}