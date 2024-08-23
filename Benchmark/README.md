# Setup
Our benchmarking toolchain uses simexpal and pandas. If necessary, install them with
```bash
pip3 install simexpal
pip3 install pandas
```

# Preparing the Instances
Before you run the loading scripts, you must compile `Runnables/InstanceLoader.cpp`. Run

```bash
cd ..
mkdir -p cmake-build-release
cd cmake-build-release
cmake .. && cmake --build . --target InstanceLoader --config Release
```

To download and prepare the instances used in our experiments, run
  ```bash
  python3 downloadAggregationInstances.py
  python3 loadParametricInstances.py
  python3 downloadVisionInstances.py
  python3 unzip.py
  python3 loadStaticInstances.py
  ```

# Using Your Own Instances
This framework supports two types of instances:
* Instances for the static max-flow-min-cut problem in DIMACS format (`.max` file extension). These are made parametric by randomly assigning parametric weights to the source-incident edges. See `Runnables/InstanceLoader.cpp` for more details on how this is done.
* Instances for the (source-sink-monotone) parametric min-cut problem in our own DIMACS-like format (`.pmax` file extension).

All inputs must be placed inside the `Data/RawInstances/` directory or a subdirectory of it. If you have instances that come in compressed `.tbz2` archives, you can unzip them by running
```bash
python3 unzip.py
```

To prepare parametric instances, run

```bash
python3 loadParametricInstances.py
```

To prepare static instances, run

```bash
python3 loadStaticInstances.py
```

A few things to note about static instances:
* The DIMACS format has no direct way of specifying infinite edge weights. They are represented by choosing a sufficiently high value such that no feasible flow can saturate them. In the vision instances used in our experiments, this value is 999999. This approach works for the static max-flow problem, but in the parametric problem it is usually possible to choose the parameter high enough that they do become saturated. To prevent this, our loading scripts allow the user to specify a threshold value for infinite weights. All weights with this value or higher will be treated as infinite. If your instances contain actual edge weights above 999999, or if they use a smaller value to represent infinite weights, you need to change the value of the `infinity` variable in `loadVisionInstances.py`.
* By default, the script assigns a random parametric weight to all source-incident edges and retains the original, static weight for all sink-incident edges. An exception are the instances that are listed in `detailed_analysis_instances` in `loadVisionInstances.py`. For these, the script will generate multiple configurations that vary the amount of parametric source- and sink-incident edges.

# Running the Experiments
Set up the experimental environment with
```bash
simex b remake
```

If you are replicating our experiments, the instance list will already be configured correctly. Otherwise, you need to modify it in `experiments.yml`. You can check that all instances are available with

```bash
simex i list
```

Launch the experiments with
```bash
simex e launch
```

Aggregate the results with
```bash
python3 eval.py
```

You can then find CSV files containing the aggregated results in the `results/` folder.