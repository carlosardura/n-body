#include <H5Cpp.h>
#include <vector>
#include <string>
#include <stdexcept>
#include "types.h"

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
            v.v_padding = 0.0;

            system.pos_mass.push_back(pm);
            system.velocities.push_back(v);
            system.ids.push_back(static_cast<uint32_t>(i));
        } 
    
    } catch (const H5::Exception& e) {
        throw std::runtime_error("Unable to read initial conditions in'" + init_file + "'.");
    }
}