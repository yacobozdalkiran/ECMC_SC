#include <omp.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "../gauge/GaugeField.h"
#include "../heatbath/heatbath.h"
#include "../io/ildg.h"
#include "../io/io.h"
#include "../observables/observables.h"

namespace fs = std::filesystem;

void generate_hb(const RunParamsHbCB& rp, bool existing) {
    //========================Objects initialization====================
    // Lattice creation + RNG
    int L = rp.L;
    GeometryCB geo(L);
    GaugeField field(geo);
    std::mt19937_64 rng(rp.seed);
    if (!rp.cold_start) {
        field.hot_start(rng);
    }

    if (existing) {
        read_ildg(field, geo, rp.run_name, rp.run_dir);

        fs::path state_path = fs::path(rp.run_dir) / rp.run_name / (rp.run_name + "_seed") /
                              (rp.run_name + "_seed.txt");
        std::ifstream ifs(state_path);
        if (ifs.is_open()) {
            ifs >> rng;
        } else {
            // Optionnel : Alerte si on attend un checkpoint mais qu'il manque un morceau
            std::cerr << "Seed file missing, using initial seed." << std::endl;
        }
    }

    // Params Hb
    HbParams h = rp.hp;
    int N_unit = 1000;

    // Measure vectors
    std::vector<double> plaquette;
    plaquette.reserve(rp.save_each / rp.N_plaquette + 2);

    io::save_params(rp, rp.run_name, rp.run_dir);

    // Print params
    print_parameters(rp);

    //==============================ECMC Checkboard===========================
    // Thermalisation
    // Skipped if existing (for successive jobs in slurm)

    if (!existing) {
        std::cout << "\n\n===========================================\n";
        std::cout << "Thermalisation : " << rp.N_therm << " configurations \n";
        std::cout << "===========================================\n";

        for (int i = 0; i < rp.N_therm; i++) {
            for (int j = 0; j < h.N_sweeps; j++) heatbath::sweep(field, geo, h.beta, h.N_hits, rng);

            // Plaquette measure (not saved for thermalization)
            if (i % rp.N_plaquette == 0) {
                double p = mean_plaquette(field, geo);
                std::cout << "\n====== Configuration " << i << " ======\n";
                std::cout << "(Therm) Sample " << i / rp.N_plaquette << ", <P> = " << p << "\n";
            }
            if (i % N_unit == 0 and i>0) {
                field.project_field_su3();
            }
        }
    }

    // Sampling
    std::cout << "\n\n===========================================\n";
    std::cout << "Sampling : " << rp.N_samples / rp.N_plaquette << " <P> samples\n "
              << "===========================================\n";

    for (int i = 0; i < rp.N_samples; i++) {
        for (int j = 0; j < h.N_sweeps; j++) heatbath::sweep(field, geo, h.beta, h.N_hits, rng);

        if (i % N_unit == 0 and i>0) {
            field.project_field_su3();
        }

        // Plaquette measure
        if ((i % rp.N_plaquette == 0) and (i > 0 or !existing)) {
            double p = mean_plaquette(field, geo);
            std::cout << "\n====== Configuration " << i << " ======\n";
            std::cout << "Sample " << i / rp.N_plaquette << ", <P> = " << p << "\n";
            plaquette.emplace_back(p);
        }

        // Save conf/seed/chain state/obs
        if (i > 0 and i % rp.save_each == 0) {
            std::cout << "\n\n==========================================\n";
            // Write the output
            int precision = 10;
            io::save_plaquette(plaquette, rp.run_name, rp.run_dir, precision);
            io::add_shift(i, rp.run_name, rp.run_dir);
            // Save conf
            save_ildg(field, geo, rp.run_name, rp.run_dir);
            // Save seeds
            io::save_seed(rng, rp.run_name, rp.run_dir);
            // Save chain state
            std::cout << "==========================================\n";
            // Clear the measures
            plaquette.clear();
            plaquette.reserve(rp.save_each / rp.N_plaquette + 2);
        }
    }

    //===========================Output======================================

    // Save conf/seed/obs
    std::cout << "\n\n==========================================\n";
    // Write the output
    int precision = 10;
    io::save_plaquette(plaquette, rp.run_name, rp.run_dir, precision);
    io::add_shift(rp.N_samples, rp.run_name, rp.run_dir);
    io::add_finished(rp.run_name, rp.run_dir);
    // Save conf
    save_ildg(field, geo, rp.run_name, rp.run_dir);
    // Save seeds
    io::save_seed(rng, rp.run_name, rp.run_dir);
    std::cout << "==========================================\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.txt>" << std::endl;
        return 1;
    }

    // Charging the parameters of the run
    RunParamsHbCB params;
    bool existing = io::read_params(params, argv[1]);

    // Measuring time
    auto start_time = std::chrono::high_resolution_clock::now();

    generate_hb(params, existing);

    auto end_time = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> total_time = end_time - start_time;
    std::cout << std::fixed << std::setprecision(4);
    print_time(total_time.count());
    // End MPI
}
