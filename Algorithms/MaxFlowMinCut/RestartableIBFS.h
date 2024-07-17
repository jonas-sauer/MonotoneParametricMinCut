#pragma once

#include <queue>
#include <vector>

#include "../../DataStructures/Graph/Graph.h"
#include "../../DataStructures/MaxFlowMinCut/MaxFlowInstance.h"

#include "../../Helpers/Assert.h"
#include "../../Helpers/Types.h"
#include "../../Helpers/Vector/Vector.h"

template<typename MAX_FLOW_INSTANCE, bool MEASUREMENTS = false>
class RestartableIBFS {

public:
    using MaxFlowInstance = MAX_FLOW_INSTANCE;
    using FlowType = MaxFlowInstance::FlowType;
    using GraphType = MaxFlowInstance::GraphType;

private:
    struct ExcessBuckets {
        ExcessBuckets(const int n) :
            buckets(n), positionOfVertex(n, -1), maxBucket(-1) {
        }

        inline void addVertex(const Vertex vertex, const int dist) noexcept {
            if (positionOfVertex[vertex] != -1) return;
            positionOfVertex[vertex] = buckets[dist].size();
            buckets[dist].emplace_back(vertex);
            maxBucket = std::max(maxBucket, dist);
        }

        inline void removeVertex(const Vertex vertex, const int dist) noexcept {
            if (positionOfVertex[vertex] == -1) return;
            const Vertex other = buckets[dist].back();
            positionOfVertex[other] = positionOfVertex[vertex];
            buckets[dist][positionOfVertex[vertex]] = other;
            buckets[dist].pop_back();
            positionOfVertex[vertex] = -1;
            if (maxBucket == dist && buckets[dist].empty()) {
                while (maxBucket >= 0 && buckets[maxBucket].empty()) maxBucket--;
            }
        }

        inline void increaseBucket(const Vertex vertex, const int oldDist, const int newDist) noexcept {
            Assert(newDist > oldDist, "Distance has not increased!");
            Assert(static_cast<size_t>(positionOfVertex[vertex]) < buckets[oldDist].size(), "Vertex is not in bucket!");
            Assert(buckets[oldDist][positionOfVertex[vertex]] == vertex, "Vertex is not in bucket!");
            const Vertex other = buckets[oldDist].back();
            positionOfVertex[other] = positionOfVertex[vertex];
            buckets[oldDist][positionOfVertex[vertex]] = other;
            buckets[oldDist].pop_back();
            positionOfVertex[vertex] = buckets[newDist].size();
            buckets[newDist].emplace_back(vertex);
            maxBucket = std::max(maxBucket, newDist);
        }

        inline bool empty() const noexcept {
            return maxBucket < 0;
        }

        inline Vertex front() const noexcept {
            Assert(!empty(), "Buckets are empty!");
            return buckets[maxBucket].back();
        }

        std::vector<std::vector<Vertex>> buckets;
        std::vector<int> positionOfVertex;
        int maxBucket;
    };

    struct OrphanBuckets {
        OrphanBuckets(const int n) :
            positionOfVertex_(n, -1), minBucket_(INFTY) {
        }

        inline void assertVertexInBucket(const Vertex vertex, const int dist) const noexcept {
            Assert(positionOfVertex_[vertex] != -1, "Vertex is not in bucket!");
            Assert(static_cast<size_t>(dist) < buckets_.size(), "Vertex is not in bucket!");
            Assert(static_cast<size_t>(positionOfVertex_[vertex]) < buckets_[dist].size(), "Vertex is not in bucket!");
            Assert(buckets_[dist][positionOfVertex_[vertex]] == vertex, "Vertex is not in bucket!");
        }

        inline bool contains(const Vertex vertex) const noexcept {
            return positionOfVertex_[vertex] != -1;
        }

        inline void addVertex(const Vertex vertex, const int dist) noexcept {
            if (positionOfVertex_[vertex] != -1) {
                assertVertexInBucket(vertex, dist);
                return;
            }
            if (static_cast<size_t>(dist) >= buckets_.size()) buckets_.resize(dist + 1);
            positionOfVertex_[vertex] = buckets_[dist].size();
            buckets_[dist].emplace_back(vertex);
            minBucket_ = std::min(minBucket_, dist);
        }

        inline void decreaseBucket(const Vertex vertex, const int oldDist, const int newDist) noexcept {
            Assert(newDist < oldDist, "Distance has not decreased!");
            assertVertexInBucket(vertex, oldDist);
            const Vertex other = buckets_[oldDist].back();
            positionOfVertex_[other] = positionOfVertex_[vertex];
            buckets_[oldDist][positionOfVertex_[vertex]] = other;
            buckets_[oldDist].pop_back();
            positionOfVertex_[vertex] = buckets_[newDist].size();
            buckets_[newDist].emplace_back(vertex);
            if (static_cast<size_t>(oldDist) == buckets_.size() - 1) {
                while (!buckets_.empty() && buckets_.back().empty()) buckets_.pop_back();
            }
            minBucket_ = std::min(minBucket_, newDist);
        }

        inline bool empty() const noexcept {
            return buckets_.empty();
        }

        inline Vertex front() noexcept {
            Assert(!empty(), "Buckets are empty!");
            return buckets_[minBucket_].back();
        }

        inline Vertex pop() noexcept {
            Assert(!empty(), "Buckets are empty!");
            const Vertex vertex = buckets_[minBucket_].back();
            buckets_[minBucket_].pop_back();
            positionOfVertex_[vertex] = -1;
            while (static_cast<size_t>(minBucket_) < buckets_.size() && buckets_[minBucket_].empty()) minBucket_++;
            if (static_cast<size_t>(minBucket_) == buckets_.size()) {
                buckets_.clear();
                minBucket_ = INFTY;
            }
            return vertex;
        }

        std::vector<std::vector<Vertex>> buckets_;
        std::vector<int> positionOfVertex_;
        int minBucket_;
    };

    struct TreeData {
        TreeData(const size_t n) :
            parentEdge(n, noEdge),
            parentVertex(n , noVertex),
            firstChild(n, noVertex),
            nextSibling(n, noVertex),
            prevSibling(n, noVertex) {
        }

        inline void addVertex(const Vertex parent, const Vertex child, const Edge edge) noexcept {
            parentEdge[child] = edge;
            parentVertex[child] = parent;
            prevSibling[child] = noVertex;
            nextSibling[child] = firstChild[parent];
            if (nextSibling[child] != noVertex)
                prevSibling[nextSibling[child]] = child;
            firstChild[parent] = child;
        }

        inline void removeVertex(const Vertex vertex) noexcept {
            if (parentVertex[vertex] == noVertex) return;
            if (nextSibling[vertex] != noVertex)
                prevSibling[nextSibling[vertex]] = prevSibling[vertex];
            if (prevSibling[vertex] == noVertex)
                firstChild[parentVertex[vertex]] = nextSibling[vertex];
            else
                nextSibling[prevSibling[vertex]] = nextSibling[vertex];
            parentEdge[vertex] = noEdge;
            parentVertex[vertex] = noVertex;
        }

        template<typename FUNCTION>
        inline void removeChildren(const Vertex parent, const FUNCTION& callback) noexcept {
            Vertex child = firstChild[parent];
            while (child != noVertex) {
                const Vertex next = nextSibling[child];
                parentEdge[child] = noEdge;
                parentVertex[child] = noVertex;
                callback(child);
                child = next;
            }
            firstChild[parent] = noVertex;
        }

        // Always directed towards sink: (parent, child) for source tree, (child, parent) for sink tree
        std::vector<Edge> parentEdge;
        std::vector<Vertex> parentVertex;
        std::vector<Vertex> firstChild;
        std::vector<Vertex> nextSibling;
        std::vector<Vertex> prevSibling;
    };

    struct Cut {
        Cut(const int n) : inSinkComponent(n, false) {}

        inline void compute(const std::vector<int>& dist) {
            for (size_t i = 0; i < dist.size(); i++) {
                inSinkComponent[i] = (dist[i] < 0);
            }
        }

        inline std::vector<Vertex> getSourceComponent() const noexcept {
            std::vector<Vertex> component;
            for (size_t i = 0; i < inSinkComponent.size(); i++) {
                if (!inSinkComponent[i]) component.emplace_back(Vertex(i));
            }
            return component;
        }

        inline std::vector<Vertex> getSinkComponent() const noexcept {
            std::vector<Vertex> component;
            for (size_t i = 0; i < inSinkComponent.size(); i++) {
                if (inSinkComponent[i]) component.emplace_back(Vertex(i));
            }
            return component;
        }

        std::vector<bool> inSinkComponent;
    };

public:
    explicit RestartableIBFS(const MaxFlowInstance& instance) :
        instance(instance),
        graph(instance.graph),
        n(graph.numVertices()),
        terminal{instance.source, instance.sink},
        residualCapacity(instance.getCurrentCapacities()),
        distance(n, 0),
        excess(n, 0),
        maxDistance{0, 0},
        currentEdge(n, noEdge),
        treeData(n),
        excessVertices{ExcessBuckets(n), ExcessBuckets(n)},
        orphans{n, n},
        threePassOrphans{n, n},
        processedOrphans(0),
        processedUniqueOrphans(0),
        orphanTimestamp(n, 0),
        currentTimestamp(0),
        threePass(false),
        cut(n) {
    }

public:
    inline void run() noexcept {
        if (MEASUREMENTS) timer.restart();
        initialize<FORWARD>();
        initialize<BACKWARD>();
        runAfterInitialize();
        if (MEASUREMENTS) flowTime += timer.elapsedMicroseconds();
    }

    inline void continueAfterUpdate() noexcept {
        if (MEASUREMENTS) timer.restart();
        std::swap(Q[FORWARD], nextQ[FORWARD]);
        std::swap(Q[BACKWARD], nextQ[BACKWARD]);
        updateCapacities();
        if (MEASUREMENTS) {
            updateTime += timer.elapsedMicroseconds();
            timer.restart();
        }
        runAfterInitialize();
        if (MEASUREMENTS) flowTime += timer.elapsedMicroseconds();
    }

    inline std::vector<Vertex> getSourceComponent() const noexcept {
        return cut.getSourceComponent();
    }

    inline std::vector<Vertex> getSinkComponent() const noexcept {
        return cut.getSinkComponent();
    }

    //TODO: Maintain the flow value throughout the algorithm.
    inline FlowType getFlowValue() const noexcept {
        FlowType flow = 0;
        for (const Vertex from : cut.getSourceComponent()) {
            for (const Edge edge : graph.edgesFrom(from)) {
                const Vertex to = graph.get(ToVertex, edge);
                if (!cut.inSinkComponent[to]) continue;
                flow += instance.getCapacity(edge);
            }
        }
        return flow;
    }

    inline double getUpdateTime() const noexcept {
        return updateTime;
    }

    inline double getFlowTime() const noexcept {
        return flowTime;
    }

private:
    template<int DIRECTION>
    inline void initialize() noexcept {
        setDistance<DIRECTION>(terminal[DIRECTION], 1);
        addExcess<DIRECTION>(terminal[DIRECTION], -INFTY);
        maxDistance[DIRECTION] = 1;
        currentEdge[terminal[DIRECTION]] = graph.beginEdgeFrom(terminal[DIRECTION]);
        nextQ[DIRECTION].emplace_back(terminal[DIRECTION]);
    }

    inline void updateCapacities() noexcept {
        Edge edgeFromSource = graph.beginEdgeFrom(terminal[FORWARD]);
        const std::vector<FlowType>& sourceDiff = instance.getSourceDiff();
        for (size_t i = 0; i < sourceDiff.size(); i++, edgeFromSource++) {
            Assert(sourceDiff[i] >= 0, "Capacity of source-incident edge has decreased!");
            if (sourceDiff[i] == 0) continue;
            residualCapacity[edgeFromSource] += sourceDiff[i];
            const Vertex to = graph.get(ToVertex, edgeFromSource);
            const FlowType add = residualCapacity[edgeFromSource];
            const Edge edgeToSource = graph.get(ReverseEdge, edgeFromSource);
            residualCapacity[edgeFromSource] = 0;
            residualCapacity[edgeToSource] += add;
            excess[to] += add;
            handleUpdateExcess(to);
        }

        Edge edgeFromSink = graph.beginEdgeFrom(terminal[BACKWARD]);
        const std::vector<FlowType>& sinkDiff = instance.getSinkDiff();
        for (size_t i = 0; i < sinkDiff.size(); i++, edgeFromSink++) {
            Assert(sinkDiff[i] <= 0, "Capacity of sink-incident edge has increased!");
            const Edge edgeToSink = graph.get(ReverseEdge, edgeFromSink);
            const Vertex from = graph.get(ToVertex, edgeFromSink);
            residualCapacity[edgeToSink] += sinkDiff[i];
            if (residualCapacity[edgeToSink] < 0) {
                const FlowType add = -residualCapacity[edgeToSink];
                residualCapacity[edgeFromSink] -= add;
                residualCapacity[edgeToSink] = 0;
                excess[from] += add;
                handleUpdateExcess(from);
            }
            if (!pmf::isNumberPositive(residualCapacity[edgeToSink]) && treeData.parentEdge[from] == edgeToSink) {
                makeOrphan<BACKWARD>(from);
            }
        }
        resetAdoptionCounter();
        adoptOrphans<BACKWARD>();
        drainExcesses<BACKWARD>();
    }

    inline void handleUpdateExcess(const Vertex vertex) noexcept {
        if (distance[vertex] == 0) {
            // Add vertex as a root to the source tree and continue the BFS from there
            setDistance<FORWARD>(vertex, maxDistance[FORWARD]);
            nextQ[FORWARD].emplace_back(vertex);
        } else if (isVertexInTree<BACKWARD>(vertex)) {
            // Register the excess for draining.
            Assert(getExcess<BACKWARD>(vertex) > 0, "Vertex does not have excess!");
            Assert(treeData.parentVertex[vertex] != noVertex, "Sink component is not a tree!");
            excessVertices[BACKWARD].addVertex(vertex, getDistance<BACKWARD>(vertex));
        } else {
            // Turn vertex into a root
            treeData.removeVertex(vertex);
        }
    }

    inline void runAfterInitialize() noexcept {
        while (true) {
            if (!grow<FORWARD>()) {
                // Continue the backward search to ensure we get the minimal sink component
                while (grow<BACKWARD>());
                break;
            }
            if (!grow<BACKWARD>()) break;
        }
        cut.compute(distance);
    }

    template<int DIRECTION>
    inline bool grow() noexcept {
        Q[DIRECTION].clear();
        std::swap(Q[DIRECTION], nextQ[DIRECTION]);
        maxDistance[DIRECTION]++;
        for (const Vertex from : Q[DIRECTION]) {
            if (getDistance<DIRECTION>(from) != maxDistance[DIRECTION] - 1) continue;

            Edge edge = graph.beginEdgeFrom(from);
            while (edge < graph.endEdgeFrom(from)) {
                const Edge edgeTowardsSink = getForwardEdge<DIRECTION>(edge);
                if (!isEdgeResidual(edgeTowardsSink)) {
                    edge++;
                    continue;
                }
                const Vertex to = graph.get(ToVertex, edge);
                if (distance[to] == 0) {
                    setDistance<DIRECTION>(to, maxDistance[DIRECTION]);
                    currentEdge[to] = graph.beginEdgeFrom(to);
                    treeData.addVertex(from, to, edgeTowardsSink);
                    nextQ[DIRECTION].emplace_back(to);
                } else if (!isVertexInTree<DIRECTION>(to)) {
                    const Vertex sourceEndpoint = (DIRECTION == BACKWARD) ? to : from;
                    const Vertex sinkEndpoint = (DIRECTION == BACKWARD) ? from : to;
                    augment(sourceEndpoint, sinkEndpoint, edgeTowardsSink);
                    if (getDistance<DIRECTION>(from) != maxDistance[DIRECTION] - 1) break;
                    if (isEdgeResidual(edgeTowardsSink)) continue;
                }
                edge++;
            }
        }
        if (nextQ[DIRECTION].empty()) {
            maxDistance[DIRECTION]--;
            return false;
        }
        return true;
    }

    inline void augment(const Vertex sourceEndpoint, const Vertex sinkEndpoint, const Edge edgeTowardsSink) noexcept {
        const Edge edgeTowardsSource = graph.get(ReverseEdge, edgeTowardsSink);
        const FlowType flow = findBottleneckCapacity(sourceEndpoint, sinkEndpoint, edgeTowardsSink);
        pushFlow<BACKWARD>(sourceEndpoint, sinkEndpoint, edgeTowardsSink, edgeTowardsSource, flow);
        resetAdoptionCounter();
        registerAndDrainExcess<FORWARD>(sourceEndpoint);
        registerAndDrainExcess<BACKWARD>(sinkEndpoint);
    }

    inline FlowType findBottleneckCapacity(const Vertex sourceEndpoint, const Vertex sinkEndpoint, const Edge edgeTowardsSink) const noexcept {
        auto[sourceBottleneck, sourceRoot] = findBottleneckCapacity(sourceEndpoint);
        FlowType bottleneck = std::min({excess[sourceRoot], sourceBottleneck, residualCapacity[edgeTowardsSink]});
        if (sourceRoot != terminal[FORWARD]) {
            auto[sinkBottleneck, sinkRoot] = findBottleneckCapacity(sinkEndpoint);
            bottleneck = std::min(bottleneck, sinkBottleneck);
        }
        return bottleneck;
    }

    inline std::pair<FlowType, Vertex> findBottleneckCapacity(const Vertex start) const noexcept {
        FlowType bottleneck = INFTY;
        Vertex vertex = start;
        while (treeData.parentVertex[vertex] != noVertex) {
            const Edge edge = treeData.parentEdge[vertex];
            bottleneck = std::min(bottleneck, residualCapacity[edge]);
            vertex = treeData.parentVertex[vertex];
        }
        return {bottleneck, vertex};
    }

    template<int DIRECTION>
    inline void registerAndDrainExcess(const Vertex start) noexcept {
        const FlowType exc = getExcess<DIRECTION>(start);
        if (exc < 0) return;
        else if (exc > 0) {
            excessVertices[DIRECTION].addVertex(start, getDistance<DIRECTION>(start));
            drainExcesses<DIRECTION>();
        }
        else {
            makeOrphan<DIRECTION>(start);
            adoptOrphans<DIRECTION>();
        }
    }

    template<int DIRECTION>
    inline void drainExcesses() noexcept {
        while (!excessVertices[DIRECTION].empty()) {
            const Vertex vertex = excessVertices[DIRECTION].front();
            Assert(hasNonNegativeExcess<DIRECTION>(vertex), "Trying to drain zero excess!");
            drainExcess<DIRECTION>(vertex);
            adoptOrphans<DIRECTION>();
        }
    }

    template<int DIRECTION>
    inline void drainExcess(Vertex vertex) noexcept {
        while (treeData.parentVertex[vertex] != noVertex) {
            const Vertex parentVertex = treeData.parentVertex[vertex];
            const Edge edgeTowardsSink = treeData.parentEdge[vertex];
            Assert(isEdgeResidualLax(edgeTowardsSink), "Tree edge is not residual!");
            const Edge edgeTowardsSource = graph.get(ReverseEdge, edgeTowardsSink);
            const FlowType exc = getExcess<DIRECTION>(vertex);
            const FlowType res = residualCapacity[edgeTowardsSink];
            const FlowType flow = std::min(res, exc);
            pushFlow<DIRECTION>(vertex, parentVertex, edgeTowardsSink, edgeTowardsSource, flow);
            if (pmf::areNumbersEqual(flow, res)) {
                makeOrphan<DIRECTION>(vertex);
            }
            if (flow == exc) {
                excessVertices[DIRECTION].removeVertex(vertex, getDistance<DIRECTION>(vertex));
            }
            vertex = parentVertex;
            if (getExcess<DIRECTION>(vertex) > 0) {
                excessVertices[DIRECTION].addVertex(vertex, getDistance<DIRECTION>(vertex));
            }
            else Assert(treeData.parentVertex[vertex] == noVertex, "Non-root vertex has zero excess!");
        }

        if (hasNonNegativeExcess<DIRECTION>(vertex))
            makeOrphan<DIRECTION>(vertex);
    }

    template<int DIRECTION>
    inline void makeOrphan(const Vertex vertex) noexcept {
        orphans[DIRECTION].addVertex(vertex, getDistance<DIRECTION>(vertex));
        treeData.removeVertex(vertex);
    }

    inline void resetAdoptionCounter() noexcept {
        processedOrphans = 0;
        processedUniqueOrphans = 0;
        currentTimestamp++;
        threePass = false;
    }

    template<int DIRECTION>
    inline void adoptOrphans() noexcept {
        if (threePass) {
            adoptOrphansThreePass<DIRECTION>();
            return;
        }
        while (!orphans[DIRECTION].empty()) {
            const Vertex orphan = orphans[DIRECTION].front();
            processedOrphans++;
            if (orphanTimestamp[orphan] != currentTimestamp) {
                processedUniqueOrphans++;
                orphanTimestamp[orphan] = currentTimestamp;
            }
            if (processedOrphans >= 3 * processedUniqueOrphans) {
                threePass = true;
                adoptOrphansThreePass<DIRECTION>();
                return;
            }
            orphans[DIRECTION].pop();
            if (adoptWithSameDistance<DIRECTION>(orphan)) continue;
            if (getDistance<DIRECTION>(orphan) == maxDistance[DIRECTION]) {
                removeOrphan<DIRECTION>(orphan);
                continue;
            }
            treeData.removeChildren(orphan, [&](const Vertex child) {
                orphans[DIRECTION].addVertex(child, getDistance<DIRECTION>(child));
            });
            if (!adoptWithNewDistance<DIRECTION>(orphan)) {
                removeOrphan<DIRECTION>(orphan);
            }
        }
    }

    template<int DIRECTION>
    inline void adoptOrphansThreePass() noexcept {
        adoptOrphansFirstPass<DIRECTION>();
        std::vector<Vertex> moved;
        while (!threePassOrphans[DIRECTION].empty()) {
            const Vertex vertex = threePassOrphans[DIRECTION].pop();
            if (treeData.parentEdge[vertex] == noEdge) {
                //std::cout << "Second pass of " << vertex << " with distance " << distance[vertex] << std::endl;
                if (adoptOrphansSecondPass<DIRECTION>(vertex)) {
                    //std::cout << "Distance increased to " << distance[vertex] << std::endl;
                    continue;
                }
            }
            if (treeData.parentEdge[vertex] != noEdge) {
                //std::cout << "Third pass of " << vertex << " with distance " << distance[vertex] << " and max distance " << maxDistance[DIRECTION] << std::endl;
                adoptOrphansThirdPass<DIRECTION>(vertex);
            } else {
                moved.emplace_back(vertex);
            }
        }
        for (const Vertex vertex : moved) {
            if (treeData.parentEdge[vertex] != noEdge) continue;
            //std::cout << "Remove "<< vertex << " with excess " << excess[vertex] << std::endl;
            removeOrphan<DIRECTION>(vertex);
        }
    }

    template<int DIRECTION>
    inline void adoptOrphansFirstPass() noexcept {
        while (!orphans[DIRECTION].empty()) {
            const Vertex vertex = orphans[DIRECTION].pop();
            excessVertices[DIRECTION].removeVertex(vertex, getDistance<DIRECTION>(vertex));
            if (adoptWithSameDistance<DIRECTION>(vertex)) continue;
            //std::cout << "First pass of " << vertex << " with distance " << distance[vertex] << " failed!" << std::endl;
            if (getDistance<DIRECTION>(vertex) == maxDistance[DIRECTION]) {
                removeOrphan<DIRECTION>(vertex);
                continue;
            }
            treeData.removeChildren(vertex, [&](const Vertex child) {
                orphans[DIRECTION].addVertex(child, getDistance<DIRECTION>(child));
            });
            const int newDistance = getDistance<DIRECTION>(vertex) + 1;
            setDistance<DIRECTION>(vertex, newDistance);
            threePassOrphans[DIRECTION].addVertex(vertex, newDistance);
        }
    }

    template<int DIRECTION>
    inline bool adoptOrphansSecondPass(const Vertex orphan) noexcept {
        int newDistance = maxDistance[DIRECTION];
        Edge newEdge = noEdge;
        Edge newEdgeTowardsSink = noEdge;
        Vertex newParent = noVertex;
        for (const Edge edge : graph.edgesFrom(orphan)) {
            const Edge edgeTowardsSink = getBackwardEdge<DIRECTION>(edge);
            if (!isEdgeResidual(edgeTowardsSink)) continue;
            const Vertex from = graph.get(ToVertex, edge);
            if (!isVertexInTree<DIRECTION>(from)) continue;
            if (isVertexFree<DIRECTION>(from)) continue;
            const int fromDistance = getDistance<DIRECTION>(from);
            if (fromDistance < newDistance) {
                newDistance = fromDistance;
                newEdge = edge;
                newEdgeTowardsSink = edgeTowardsSink;
                newParent = from;
            }
        }
        if (newEdge == noEdge) {
            //std::cout << "Adoption failed" << std::endl;
            return false;
        }
        treeData.parentEdge[orphan] = newEdgeTowardsSink;
        if (newDistance + 1 > getDistance<DIRECTION>(orphan)) {
            setDistance<DIRECTION>(orphan, newDistance + 1);
            threePassOrphans[DIRECTION].addVertex(orphan, newDistance + 1);
            return true;
        }
        return false;
    }

    template<int DIRECTION>
    inline void adoptOrphansThirdPass(const Vertex vertex) noexcept {
        for (const Edge edge : graph.edgesFrom(vertex)) {
            const Edge edgeTowardsSink = getBackwardEdge<DIRECTION>(edge);
            if (!isEdgeResidual(edgeTowardsSink)) continue;
            const Vertex from = graph.get(ToVertex, edge);
            if (!isEdgeAdmissible<DIRECTION>(vertex, from)) continue;
            if (isVertexFree<DIRECTION>(from)) continue;
            const int dist = getDistance<DIRECTION>(vertex);
            if (getExcess<DIRECTION>(vertex) > 0)
                excessVertices[DIRECTION].addVertex(vertex, dist);
            treeData.addVertex(from, vertex, edgeTowardsSink);
            currentEdge[vertex] = edge;
            if (dist == maxDistance[DIRECTION]) {
                nextQ[DIRECTION].emplace_back(vertex);
            }
            break;
        }

        //Don't try to adopt children beyond maxDistance. This will be done in the growth steps.
        if (getDistance<DIRECTION>(vertex) == maxDistance[DIRECTION]) return;
        for (const Edge edge : graph.edgesFrom(vertex)) {
            const Vertex from = graph.get(ToVertex, edge);
            if (from == terminal[DIRECTION] || isVertexInTree<!DIRECTION>(from)) continue;
            const bool isFree = isVertexFree<DIRECTION>(from);
            const int fromDist = getDistance<DIRECTION>(from);
            const int vDist = getDistance<DIRECTION>(vertex);
            const bool isDistanceGreater = fromDist > vDist + 1;
            if (!isFree && !isDistanceGreater) continue;
            const Edge edgeTowardsSink = getForwardEdge<DIRECTION>(edge);
            if (!isEdgeResidual(edgeTowardsSink)) continue;
            if (isDistanceGreater) {
                //std::cout << "Decrease distance of " << from << " from " << fromDist << " to " << vDist + 1 << std::endl;
                threePassOrphans[DIRECTION].decreaseBucket(from, fromDist, vDist + 1);
            } else {
                assert(isFree);
                //std::cout << "Re-add orphan " << from << " from " << fromDist << " to " << vDist + 1 << std::endl;
                if (fromDist == vDist) {
                    if (threePassOrphans[DIRECTION].contains(from)) continue;
                    threePassOrphans[DIRECTION].addVertex(from, vDist + 1);
                } else if (fromDist < vDist) {
                    threePassOrphans[DIRECTION].addVertex(from, vDist + 1);
                }
            }
            treeData.parentEdge[from] = edgeTowardsSink;
            setDistance<DIRECTION>(from, vDist + 1);
        }
    }

    template<int DIRECTION>
    inline bool adoptWithSameDistance(const Vertex orphan) noexcept {
        for (Edge edge = currentEdge[orphan]; edge < graph.endEdgeFrom(orphan); edge++) {
            const Edge edgeTowardsSink = getBackwardEdge<DIRECTION>(edge);
            if (!isEdgeResidual(edgeTowardsSink)) continue;
            const Vertex from = graph.get(ToVertex, edge);
            if (distance[from] == 0 || !isEdgeAdmissible<DIRECTION>(orphan, from)) continue;
            treeData.addVertex(from, orphan, edgeTowardsSink);
            currentEdge[orphan] = edge;
            return true;
        }
        return false;
    }

    template<int DIRECTION>
    inline bool adoptWithNewDistance(const Vertex orphan) noexcept {
        const int oldDistance = getDistance<DIRECTION>(orphan);
        int newDistance = maxDistance[DIRECTION];
        Edge newEdge = noEdge;
        Edge newEdgeTowardsSink = noEdge;
        Vertex newParent = noVertex;
        for (const Edge edge : graph.edgesFrom(orphan)) {
            const Edge edgeTowardsSink = getBackwardEdge<DIRECTION>(edge);
            if (!isEdgeResidual(edgeTowardsSink)) continue;
            const Vertex from = graph.get(ToVertex, edge);
            if (!isVertexInTree<DIRECTION>(from)) continue;
            const int fromDistance = getDistance<DIRECTION>(from);
            if (fromDistance < newDistance) {
                newDistance = fromDistance;
                newEdge = edge;
                newEdgeTowardsSink = edgeTowardsSink;
                newParent = from;
                if (newDistance == oldDistance) break;
            }
        }
        if (newEdge == noEdge) return false;
        setDistance<DIRECTION>(orphan, newDistance + 1);
        if (getExcess<DIRECTION>(orphan) > 0)
            excessVertices[DIRECTION].increaseBucket(orphan, oldDistance, newDistance + 1);
        currentEdge[orphan] = newEdge;
        treeData.addVertex(newParent, orphan, newEdgeTowardsSink);
        if (newDistance + 1 == maxDistance[DIRECTION])
            nextQ[DIRECTION].emplace_back(orphan);
        return true;
    }

    template<int DIRECTION>
    inline void removeOrphan(const Vertex orphan) noexcept {
        if (getExcess<DIRECTION>(orphan) > 0) {
            excessVertices[DIRECTION].removeVertex(orphan, getDistance<DIRECTION>(orphan));
            setDistance<!DIRECTION>(orphan, maxDistance[!DIRECTION]);
            nextQ[!DIRECTION].emplace_back(orphan);
            currentEdge[orphan] = graph.beginEdgeFrom(orphan);
        } else {
            distance[orphan] = 0;
        }
    }

    template<int DIRECTION>
    inline void pushFlow(const Vertex from, const Vertex to, const Edge edge, const Edge reverseEdge, const FlowType flow) noexcept {
        addExcess<DIRECTION>(from, -flow);
        addExcess<DIRECTION>(to, flow);
        residualCapacity[edge] -= flow;
        residualCapacity[reverseEdge] += flow;
    }

    template<int DIRECTION>
    inline Edge getForwardEdge(const Edge edge) noexcept {
        if (DIRECTION == BACKWARD)
            return graph.get(ReverseEdge, edge);
        else
            return edge;
    }

    template<int DIRECTION>
    inline Edge getBackwardEdge(const Edge edge) const noexcept {
        if (DIRECTION == FORWARD)
            return graph.get(ReverseEdge, edge);
        else
            return edge;
    }

    template<int DIRECTION>
    inline int getDistance(const Vertex vertex) const noexcept {
        return (DIRECTION == BACKWARD) ? -distance[vertex] : distance[vertex];
    }

    template<int DIRECTION>
    inline void setDistance(const Vertex vertex, const int dist) noexcept {
        distance[vertex] = (DIRECTION == BACKWARD) ? -dist : dist;
    }

    template<int DIRECTION>
    inline FlowType getExcess(const Vertex vertex) const noexcept {
        return (DIRECTION == FORWARD) ? -excess[vertex] : excess[vertex];
    }

    template<int DIRECTION>
    inline void addExcess(const Vertex vertex, const FlowType add) noexcept {
        if (DIRECTION == FORWARD)
            excess[vertex] -= add;
        else
            excess[vertex] += add;
    }

    template<int DIRECTION>
    inline bool isVertexInTree(const Vertex vertex) const noexcept {
        return getDistance<DIRECTION>(vertex) > 0;
    }

    template<int DIRECTION>
    inline bool isVertexFree(const Vertex vertex) const noexcept {
        // Vertices without a parent are either free or roots. Roots have negative excess, free vertices do not.
        return treeData.parentEdge[vertex] == noEdge && getExcess<DIRECTION>(vertex) >= 0;
    }

    template<int DIRECTION>
    inline bool isEdgeAdmissible(const Vertex from, const Vertex to) const noexcept {
        return getDistance<DIRECTION>(from) == getDistance<DIRECTION>(to) + 1;
    }

    inline bool isEdgeResidual(const Edge edge) const noexcept {
        return pmf::isNumberPositive(residualCapacity[edge]);
    }

    inline bool isEdgeResidualLax(const Edge edge) const noexcept {
        return !pmf::isNumberNegative(residualCapacity[edge]);
    }

    template<int DIRECTION>
    inline bool hasNonNegativeExcess(const Vertex vertex) const noexcept {
        return !pmf::isNumberNegative(getExcess<DIRECTION>(vertex));
    }

    template<int DIRECTION>
    inline void checkBuckets() const noexcept {
        for (size_t i = 0; static_cast<int>(i) <= excessVertices[DIRECTION].maxBucket; i++) {
            for (const Vertex vertex : excessVertices[DIRECTION].buckets[i]) {
                Assert(getExcess<DIRECTION>(vertex) > 0, "Vertex in bucket has no excess!");
                Assert(getDistance<DIRECTION>(vertex) == int(i), "Vertex is in wrong bucket!");
            }
        }
    }

    template<int DIRECTION>
    inline void checkTree() const noexcept {
        for (const Vertex v : graph.vertices()) {
            checkTree<DIRECTION>(v);
        }
    }

    template<int DIRECTION>
    inline void checkTree(const Vertex v) const noexcept {
        if (!isVertexInTree<DIRECTION>(v)) return;
        const Edge edge = treeData.parentEdge[v];
        if (edge == noEdge) return;
        const Vertex parent = treeData.parentVertex[v];
        assert(getDistance<DIRECTION>(v) == getDistance<DIRECTION>(parent) + 1);
        for (const Edge e : graph.edgesFrom(v)) {
            const Edge edgeTowardsSink = getBackwardEdge<DIRECTION>(e);
            const Vertex nonParent = graph.get(ToVertex, e);
            if (!isEdgeResidual(edgeTowardsSink)) continue;
            if (!isVertexInTree<DIRECTION>(nonParent)) continue;
            assert(getDistance<DIRECTION>(v) <= getDistance<DIRECTION>(nonParent) + 1);
        }
    }

    inline void checkDistanceInvariants(const bool allowOrphans = false) const noexcept {
        for (const Vertex vertex : graph.vertices()) {
            checkDistanceInvariants(vertex, allowOrphans);
        }
    }

    inline void checkDistanceInvariants(const Vertex vertex, const bool allowOrphans = false) const noexcept {
        const Vertex parent = treeData.parentVertex[vertex];
        if (abs(distance[vertex]) <= 1) {
            Assert(parent == noVertex, "Vertex " << vertex << " with distance <= 1 has parent " << parent << "!");
        } else {
            if (!graph.isVertex(parent)) {
                Assert(allowOrphans, "Vertex " << vertex << " has an invalid parent!");
            } else {
                Assert(abs(distance[vertex]) >= abs(distance[parent]) + 1, "Distance of " << vertex << " is " << distance[vertex] << ", but distance of " << parent << " is " << distance[parent] << "!");
            }
        }
    }

    inline void checkChildrenRelation() const noexcept {
        std::vector<std::vector<Vertex>> childrenByParent(n);
        for (const Vertex vertex : graph.vertices()) {
            if (treeData.parentVertex[vertex] != noVertex) {
                childrenByParent[treeData.parentVertex[vertex]].emplace_back(vertex);
            }
        }
        for (const Vertex vertex : graph.vertices()) {
            std::vector<Vertex> children;
            Vertex child = treeData.firstChild[vertex];
            while (child != noVertex) {
                children.emplace_back(child);
                child = treeData.nextSibling[child];
            }
            std::sort(children.begin(), children.end());
            std::sort(childrenByParent[vertex].begin(), childrenByParent[vertex].end());
            Ensure(Vector::equals(children, childrenByParent[vertex]), "Child and parent relations are inconsistent!");
        }
    }

    inline void checkChildrenRelation(const Vertex vertex) const noexcept {
        std::vector<bool> isChild(n, false);
        Vertex child = treeData.firstChild[vertex];
        while (child != noVertex) {
            isChild[child] = true;
            child = treeData.nextSibling[child];
        }
        for (const Vertex c : graph.vertices()) {
            const Vertex parent = treeData.parentVertex[c];
            if (isChild[c]) {
                Ensure(parent == vertex, "Child " << c << " of << " << vertex << " has the wrong parent!");
            } else {
                Ensure(parent != vertex, "Vertex " << c << " has parent << " << vertex << " but is not a child!");
            }
        }
    }

    inline void checkTreeResidual() const noexcept {
        for (const Vertex vertex : graph.vertices()) {
            if (treeData.parentEdge[vertex] == noEdge) continue;
            Assert(isEdgeResidual(treeData.parentEdge[vertex]), "Tree edge is not residual!");
        }
    }

    inline void checkFlowConservation() const noexcept {
        for (const Vertex vertex : graph.vertices()) {
            checkFlowConservation(vertex);
        }
    }

    inline FlowType getInflow(const Vertex vertex) const noexcept {
        FlowType inflow = 0;
        for (const Edge edge : graph.edgesFrom(vertex)) {
            const Edge reverseEdge = graph.get(ReverseEdge, edge);
            inflow += instance.getCapacity(reverseEdge) - residualCapacity[reverseEdge];
        }
        return inflow;
    }

    inline void checkFlowConservation(const Vertex vertex) const noexcept {
        if (vertex == terminal[FORWARD] || vertex == terminal[BACKWARD]) return;
        if (distance[vertex] < 0) Assert(excess[vertex] >= 0, "Vertex in sink component has a deficit!");
        const FlowType inflow = getInflow(vertex);
        Assert(pmf::areNumbersEqual(inflow, excess[vertex]), "Flow conservation not fulfilled!");
    }

    inline void checkCapacityConstraints() const noexcept {
        for (const Edge edge : graph.edges()) {
            Assert(residualCapacity[edge] >= 0, "Capacity constraint violated!");
        }
    }

    template<int DIRECTION>
    inline void checkQueue() noexcept {
        for (const Vertex vertex : graph.vertices()) {
            if (getDistance<DIRECTION>(vertex) != maxDistance[DIRECTION]) continue;
            Assert(Vector::contains(Q[DIRECTION], vertex), "Vertex missing from queue");
        }
    }

    inline void checkTreesDisconnected() noexcept {
        for (const Vertex vertex : graph.vertices()) {
            if (distance[vertex] <= 0) continue;
            for (const Edge edge : graph.edgesFrom(vertex)) {
                const Vertex to = graph.get(ToVertex, edge);
                if (distance[to] >= 0) continue;
                if (!isEdgeResidual(edge)) continue;
                Assert(false, "Found connecting edge " << vertex << " -> " << to << " with residual capacity " << residualCapacity[edge] << " and distances " << distance[vertex] << "/" << maxDistance[FORWARD] << " and " << -distance[to]  << "/" << maxDistance[BACKWARD]);
            }
        }
    }

private:
    const MaxFlowInstance& instance;
    const GraphType& graph;
    const int n;
    const Vertex terminal[2];
    std::vector<FlowType> residualCapacity;
    // positive for s-vertices, negative for t-vertices, 0 for n-vertices
    std::vector<int> distance;
    // positive for source roots, negative for sink roots
    // non-positive for other s-vertices, non-negative for other t-vertices
    std::vector<FlowType> excess;
    int maxDistance[2];
    std::vector<Edge> currentEdge;
    TreeData treeData;
    std::vector<Vertex> Q[2];
    std::vector<Vertex> nextQ[2];
    ExcessBuckets excessVertices[2];
    OrphanBuckets orphans[2];
    OrphanBuckets threePassOrphans[2];
    size_t processedOrphans;
    size_t processedUniqueOrphans;
    std::vector<int> orphanTimestamp;
    int currentTimestamp;
    bool threePass;
    Cut cut;

    double updateTime = 0;
    double flowTime = 0;
    Timer timer;
};
