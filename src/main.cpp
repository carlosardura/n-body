#include <H5Cpp.h>
#include <vector>
#include <iostream>
#include <string>
#include <stdexcept>
#include <thread>
#include <cmath>
#include <algorithm>
#include "types.h"
#include "integrator.h"
#include "barneshut.h"
#include "morton_code.h"
#include "io.h"

OctantDomain compute_bounds(const BodySystem& system) {
    if (system.size() == 0) return {};

    double min_x = system.x[0], max_x = system.x[0];
    double min_y = system.y[0], max_y = system.y[0];
    double min_z = system.z[0], max_z = system.z[0];

    for (size_t i = 1; i < system.size(); ++i) {
        if (system.x[i] < min_x) min_x = system.x[i]; if (system.x[i] > max_x) max_x = system.x[i];
        if (system.y[i] < min_y) min_y = system.y[i]; if (system.y[i] > max_y) max_y = system.y[i];
        if (system.z[i] < min_z) min_z = system.z[i]; if (system.z[i] > max_z) max_z = system.z[i];
    }

    double span_x = max_x - min_x;
    double span_y = max_y - min_y;
    double span_z = max_z - min_z;

    OctantDomain bbox;
    bbox.xmin = min_x;
    bbox.ymin = min_y;
    bbox.zmin = min_z;
    bbox.max_span = std::max({span_x, span_y, span_z});
    bbox.size = bbox.max_span;
    bbox.mid_x = min_x + bbox.size / 2.0;
    bbox.mid_y = min_y + bbox.size / 2.0;
    bbox.mid_z = min_z + bbox.size / 2.0;

    return bbox;
}


int main(int argc, char* argv[]) {
    if (argc != 10) {
        return 1; 
    }

    std::string init_file = argv[1];
    std::string output_dir = argv[2];
    size_t n = std::stoul(argv[3]);
    double G = std::stod(argv[4]);
    int total_steps = std::stoi(argv[5]);
    double dt = std::stod(argv[6]);
    int dump_freq = std::stoi(argv[7]);
    double theta = std::stod(argv[8]);
    double epsilon = std::stod(argv[9]);

    BodySystem system;
    system.resize(n);

    BodySystem morton_system;
    morton_system.resize(n);

    std::vector<double> ax(n, 0.0);
    std::vector<double> ay(n, 0.0);    
    std::vector<double> az(n, 0.0);

    try {
        H5::H5File file(init_file, H5F_ACC_RDONLY);

        H5::DataSet dataset_pos = file.openDataSet("positions");
        H5::DataSet dataset_vel = file.openDataSet("velocities");
        H5::DataSet dataset_mass = file.openDataSet("masses");

        std::vector<double> rbuffer(n * 3);
        std::vector<double> vbuffer(n * 3);
        std::vector<double> mbuffer(n);

        dataset_pos.read(rbuffer.data(), H5::PredType::NATIVE_DOUBLE);
        dataset_vel.read(vbuffer.data(), H5::PredType::NATIVE_DOUBLE);
        dataset_mass.read(mbuffer.data(), H5::PredType::NATIVE_DOUBLE);

        for (size_t i = 0; i < n; ++i) {
            system.x[i] = rbuffer[i * 3];
            system.y[i] = rbuffer[i * 3 + 1];
            system.z[i] = rbuffer[i * 3 + 2];
            system.mass[i] = mbuffer[i];

            system.vx[i] = vbuffer[i * 3];
            system.vy[i] = vbuffer[i * 3 + 1];
            system.vz[i] = vbuffer[i * 3 + 2];

            system.ids[i] = static_cast<uint32_t>(i);
        } 
    
    } catch (const H5::Exception& e) {
        throw std::runtime_error("Unable to read initial conditions in'" + init_file + "'.");
    }

    BarnesHut barnes_hut;
    std::thread io_thread;

    OctantDomain bbox = compute_bounds(system);
    MortonCode::sort_system(system, morton_system, bbox);
    barnes_hut.partition_space(system, bbox);
    barnes_hut.compute_force(system, ax, ay, az, theta, epsilon, G);


    for (int step = 0; step < total_steps; ++step) {
        Integrator::velverlet_step1(system, ax, ay, az, dt);
        OctantDomain bbox = compute_bounds(system);
        if (bbox.max_span == 0.0) break;

        MortonCode::sort_system(system, morton_system, bbox);

        barnes_hut.partition_space(system, bbox);
        barnes_hut.compute_force(system, ax, ay, az, theta, epsilon, G);

        Integrator::velverlet_step2(system, ax, ay, az, dt);
    
        if (step % dump_freq == 0) {
            if (io_thread.joinable()) {
                io_thread.join();
            }
            
            auto system_snapshot = std::make_shared<BodySystem>(system);
            
            io_thread = std::thread([output_dir, system_snapshot, step, total_steps]() {
                IO::dump_frame(output_dir, *system_snapshot, step, total_steps);
            });
        }
    }
    if(io_thread.joinable()) {
        io_thread.join();
    }

    return 0;
}