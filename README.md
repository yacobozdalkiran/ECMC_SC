## Event-Chain Monte Carlo implementation for quenched QCD (No parallelization)

This repo contains an implementation of Event-Chain Monte Carlo, Metropolis and Heatbath for pure gauge 4d lattice QCD.

Written in C++ using Eigen, compiled using CMake.

### Compilation

Once in the root folder of the repo, execute the following commands to compile with CMake.

```bash
mkdir build
cd build
cmake ..
make -j 4
```

### Execution

The exectables are 

* `gauge_ecmc_norev` for Forward ECMC
* `gauge_ecmc` for Reflective ECMC
* `gauge_heatbath` for Heatbath
* `gauge_metropolis` for Metropolis

For each of these, you need to pass an input .txt file to run the simulation. Two example input files `etest.txt` and `htest.txt` are in the folder `inputs/`, resp. for ECMC and Heatbath executables.

Example :

```bash
./build/gauge_ecmc inputs/etest.txt
```

The output files are located in the `run_dir/run_name` folder following the parameters of the input file.
