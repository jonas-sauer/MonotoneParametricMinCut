# ParametricMaxFlow
This repository contains algorithms for the source-sink monotone parametric minimum s-t-cut problem.

# Setup
Install the googletest submodule by running
```bash
git submodule init
git submodule update
```
To compile all executables, run
```bash
mkdir -p cmake-build-release
cd cmake-build-release
cmake .. -DCMAKE_BUILD_TYPE=Release && cmake --build . --config Release
```

Use ```InstanceLoader``` to load max-flow instances and ```Benchmark``` to run the algorithms. Additional experiments to evaluate the precision of the computed solutions are contained in ```PrecisionExperiment```.

# Benchmarking
A guide for running the benchmark experiments can be found in `Benchmark/README.md`.