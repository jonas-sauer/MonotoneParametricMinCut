# How to run the benchmark experiments
Build the executable with
```bash
cd .. # This is because this readme is not in the main directory.
mkdir -p build
cd build
cmake .. && cmake --build . --target Benchmark --config Release
```

Run the bash script generator
```bash
python3 PMFmakeBashScript.py
```

This generates a bash script that will automatically run the benchmarks
```bash
/bin/bash PMFrunBenchmark.sh
```


# How to modify/expand the benchmark
Read the PMFrunBenchmarkConfig.txt. It contains an explanation and the benchmark settings.