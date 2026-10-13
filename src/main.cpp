#include <H5Cpp.h>
#include <vector>
#include <iostream>
#include <string>
#include <stdexcept>
#include "types.h"
#include "integrator.h"
#include "barneshut.h"
#include "morton_code.h"
#include "io.h"

Cell compute_bounds(const std::vector<PosMass>& pos_mass) {
    if (pos_mass.empty()) return {};

    double min_x = pos_mass[0].x, max_x = pos_mass[0].x;
    double min_y = pos_mass[0].y, max_y = pos_mass[0].y;
    double min_z = pos_mass[0].z, max_z = pos_mass[0].z;

    for (const auto& p : pos_mass) {
        if (p.x < min_x) min_x = p.x; if (p.x > max_x) max_x = p.x;
        if (p.y < min_y) min_y = p.y; if (p.y > max_y) max_y = p.y;
        if (p.z < min_z) min_z = p.z; if (p.z > max_z) max_z = p.z;
    }

    double span_x = max_x - min_x;
    double span_y = max_y - min_y;
    double span_z = max_z - min_z;

    Cell bbox;
    bbox.xmin = min_x;
    bbox.ymin = min_y;
    bbox.zmin = min_z;
    bbox.max_span = std::max({span_x, span_y, span_z});

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
    system.reserve(n);
    std::vector<PosMass> accelerations(n, {0.0, 0.0, 0.0, 0.0});


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
            PosMass pm;
            pm.x = rbuffer[i * 3];
            pm.y = rbuffer[i * 3 + 1];
            pm.z = rbuffer[i * 3 + 2];
            pm.mass = mbuffer[i];

            Velocity v;
            v.vx = vbuffer[i * 3];
            v.vy = vbuffer[i * 3 + 1];
            v.vz = vbuffer[i * 3 + 2];
            v.padding = 0.0;

            system.pos_mass.push_back(pm);
            system.velocities.push_back(v);
            system.ids.push_back(static_cast<uint32_t>(i));

        } 
    
    } catch (const H5::Exception& e) {
        throw std::runtime_error("Unable to read initial conditions in'" + init_file + "'.");
    }
}