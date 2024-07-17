# How to run the benchmark for parametric Max Flow
Build the executable with
```bash
cd .. # This is because this ReadMe is not in the main directory.
mkdir -p build
cd build
cmake .. && cmake --build . --target ParametricMaxFlowBenchmark
```

Run the BashScriptGenerator
```bash
python3 PMFmakeBashScript.py
```

This generates a bash script that will automatically execute the benchmarking
```bash
/bin/bash PMFrunBenchmark.sh
```


# How to modify/expand the benchmark
Read the PMFrunBenchmarkConfig.txt. It contains an explanation and 