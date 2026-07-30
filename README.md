## Event-Chain Monte Carlo implementation for quenched QCD (No parallelization)

This repo contains an implementation of Event-Chain Monte Carlo and Heatbath for pure gauge 4d lattice QCD.

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

The exectables are `gauge_ecmc_norev`, `gauge_ecmc` and `gauge_heatbath`

