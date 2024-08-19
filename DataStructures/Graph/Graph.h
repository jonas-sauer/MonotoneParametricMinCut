#pragma once

#include "Classes/GraphInterface.h"

#include "Utils/Utils.h"

#include "Classes/DynamicGraph.h"
#include "Classes/StaticGraph.h"
#include "Classes/EdgeList.h"

using NoVertexAttributes = List<>;
using NoEdgeAttributes = List<>;
template<typename CAP>
using WithReverseEdgesAndCapacity = List<Attribute<ReverseEdge, Edge>, Attribute<Capacity, CAP>>;
template<typename CAP>
using WithParameterizedCapacity = List<Attribute<Capacity, CAP>>;

template<typename CAP>
using ParametricFlowGraph = StaticGraph<NoVertexAttributes, WithReverseEdgesAndCapacity<CAP>>;
template<typename CAP>
using ParametricFlowGraphEdgeList = EdgeList<NoVertexAttributes, WithParameterizedCapacity<CAP>>;
template<typename CAP>
using DynamicParametricFlowGraph = DynamicGraph<NoVertexAttributes, WithReverseEdgesAndCapacity<CAP>>;

#include "Utils/Conversion.h"
#include "Utils/IO.h"
