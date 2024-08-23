# MonotoneParametricMinCut
This repository contains algorithms for the parametric minimum s-t-cut problem with source-sink-monotone edge capacities. It is developed by the [Computational Analytics](https://ca.cs.uni-bonn.de) group at the University of Bonn. Auxiliary code, including the graph data structures, is taken from the [ULTRA](https://github.com/kit-algo/ULTRA) framework developed by the KIT Algorithmics group.

## What is this for?
The (static) minimum s-t-cut problem asks for a cut with minimal capacity between two vertices s and t in a flow network. This is equivalent to finding a maximum s-t-flow. In the parametric version of the problem, the edge capacities are functions of a single parameter $\lambda$. The objective is to find a minimum s-t-cut for every value of $\lambda$ within some interval $[\lambda_\text{min}, \lambda_\text{max}]$.

This framework provides algorithms for the special case in which the functions are *source-sink-monotone*. This means they are non-decreasing for edges incident to s, non-increasing for edges incident to t, and constant otherwise. Here, the cuts are *nested*: the source component of a minimum cut for $\lambda_1$ is fully contained in the source component of a minimum cut for $\lambda_2 > \lambda_1$. We therefore represent the solution as a *breakpoint function*: for each vertex, it provides the value of $\lambda$ for which the vertex switches from the sink component to the source component.

This framework was originally developed to solve [polygon aggregation for map simplification](https://doi.org/10.4230/LIPIcs.GIScience.2021.II.6), which can be formulated as a monotone parametric min-cut problem. Other applications include [various tasks in computer vision](https://doi.org/10.1109/ICCV.2007.4408910), such as binary image segmentation, multiview reconstruction, and surface fitting.

## Algorithms
The folder ``Algorithms/StaticMinCut/`` contains static max-flow-min-cut algorithms:
* `IBFS.h`: The [incremental breadth-first search](https://doi.org/10.1007/978-3-642-23719-5\_39) algorithm by Goldberg, Hed, Kaplan, Tarjan and Werneck.
* `EIBFS.h`: The [IBFS variant with vertex excesses](https://doi.org/10.1007/978-3-662-48350-3\_52) by the same authors plus Kohli. Our implementations of both IBFS variants are based on the implementation [provided by the authors](https://www.cs.tau.ac.il/~sagihed/ibfs/). Note that the original link no longer works. A more recent fork can be found [here](https://github.com/PolarNick239/IBFS/tree/master), for example.
* `PushRelabel.h`: The classic [push-relabel](https://doi.org/10.1145/48014.61051) algorithm by Goldberg and Tarjan. Our implementation is based on the one from the [Boost graph library](https://www.boost.org/doc/libs/1_77_0/libs/graph/doc/push_relabel_max_flow.html), which in turn is based on the implementation by [Cherkassky and Goldberg](https://doi.org/10.1007/PL00009180).

Algorithms for the monotone parametric min-cut problem can be found under ``Algorithms/ParametricMinCut``:
* `DichotomicScheme.h`: An algorithm based on the classic dichotomic search technique for bicriteria optimization problems. It was first described by [Eisner and Severance](https://doi.org/10.1145/321978.321982) for this problem. The dichotomic scheme makes $O(n)$ calls to a static max-flow algorithm, which must be provided as the template parameter `SEARCH_ALGORITHM`. The graph is successively shrunk by contracting vertices that are known to lie in the source or sink component.
* `DichotomicSchemeNoContraction.h`: A variant of the dichotomic scheme that does not shrink the graph. This is the algorithm proposed for the polygon aggregation problem by [Rottmann et al](https://doi.org/10.4230/LIPIcs.GIScience.2021.II.6). Note that this variant is prohibitively slow on large instances and is only provided for comparison.
* `ParametricIBFS.h`: Our new IBFS-based algorithm.

All parametric algorithms are supplied with a template parameter `FLOW_FUNCTION`. This should be a class derived from `FlowFunction` (see `DataStructures/MaxFlow/FlowFunction.h`) that describes the edge capacity functions and how they are evaluated. Our implementation uses the class `LinearFlowFunction`, which describes functions of the form $f(x) = a \cdot x + b$.

## Setup
Install the googletest submodule by running
```bash
git submodule init
git submodule update
```
To compile all executables, run
```bash
mkdir -p cmake-build-release
cd cmake-build-release
cmake .. && cmake --build . --config Release
```

## Usage
The algorithms can be run with the executable ```Benchmark```. It expects the input networks to be in a special binary format. To convert max-flow instances in [DIMACS format](https://lpsolve.sourceforge.net/5.5/DIMACS_maxf.htm) to this binary format, use the executable ```InstanceLoader```. Launch the executables without arguments to get specific usage instructions.

## Benchmarking
A guide for running the benchmark experiments can be found in `Benchmark/README.md`.