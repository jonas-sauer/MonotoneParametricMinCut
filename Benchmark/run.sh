cd ..
mkdir -p cmake-build-release
cd cmake-build-release
cmake .. -DCMAKE_BUILD_TYPE=Release && cmake --build . --target InstanceLoader --config Release
cd ../Benchmark/
python3 downloadAggregationInstances.py
python3 loadParametricInstances.py
python3 downloadVisionInstances.py
python3 unzipInstances.py
python3 loadStaticInstances.py
simex b remake
simex e launch
python3 eval.py
