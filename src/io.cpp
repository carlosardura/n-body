#include <H5Cpp.h>
#include <vector>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include "io.h"

namespace IO {
    void dump_frame(const std::string& directory, const BodySystem& system, int step, int total_steps) {
        size_t n = system.size();
        
        std::vector<double> out_r(n * 3);
        std::vector<double> out_v(n * 3);

        for (size_t i = 0; i < n; ++i) {
            out_r[i * 3] = system.pos_mass[i].x;
            out_r[i * 3 + 1] = system.pos_mass[i].y;
            out_r[i * 3 + 2] = system.pos_mass[i].z;

            out_v[i * 3] = system.velocities[i].vx;
            out_v[i * 3 + 1] = system.velocities[i].vy;
            out_v[i * 3 + 2] = system.velocities[i].vz;
        }

        int width = std::to_string(total_steps).length();
        std::ostringstream filename;
        filename << directory << "/output_" << std::setfill('0') << std::setw(width) << step << ".h5";

        try {
            H5::H5File file(filename.str(), H5F_ACC_TRUNC);
            
            hsize_t dims[2] = {n, 3};
            H5::DataSpace dataspace(2, dims);
            
            H5::DataSet dataset_pos = file.createDataSet("positions", H5::PredType::NATIVE_DOUBLE, dataspace);
            dataset_pos.write(out_r.data(), H5::PredType::NATIVE_DOUBLE);
            
            H5::DataSet dataset_vel = file.createDataSet("velocities", H5::PredType::NATIVE_DOUBLE, dataspace);
            dataset_vel.write(out_v.data(), H5::PredType::NATIVE_DOUBLE);
            
        } catch (const H5::Exception& e) {
            throw std::runtime_error("Unable to write'" + filename.str() + "'data to file.");
        }
    }
}