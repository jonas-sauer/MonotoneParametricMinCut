#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <limits>

#include "../../../DataStructures/CSA/Data.h"
#include "../../../DataStructures/Assignment/Settings.h"

#include "../../../Helpers/BinarySearch.h"
#include "../../../Helpers/Vector/Vector.h"
#include "Beta.h"
#include "Kumaraswamy.h"
#include "LogitNormal.h"
#include "PERT.h"

namespace Assignment::MRUM {

template<uint32_t NUMBER_OF_SAMPLES, typename PROFILER>
class SimpleSamplePATComputation {

public:
    inline static constexpr uint32_t NumberOfSamples = NUMBER_OF_SAMPLES;
    using Profiler = PROFILER;
    using Type = SimpleSamplePATComputation<NumberOfSamples, Profiler>;

    using WithOriginalEdgeAndProfileIndex = List<Attribute<OriginalEdge, uint32_t>, Attribute<ProfileIndex, uint32_t>>;
    using DynamicConnectionGraph = DynamicGraph<NoVertexAttributes, WithOriginalEdgeAndProfileIndex>;
    using ConnectionGraph = StaticGraph<NoVertexAttributes, WithOriginalEdgeAndProfileIndex>;

    using TimeType = double;
    using FactorType = double;
    using TimeSample = std::array<TimeType, NumberOfSamples>;
    using FactorSample = std::array<FactorType, NumberOfSamples>;
    using ParentSample = std::array<ConnectionId, NumberOfSamples>;

    inline static constexpr TimeType NotReachable = 100000000;

    template<typename NEW_TIME, typename NEW_PARENT>
    inline static void updateLabel(TimeSample& times, ParentSample& parents, const NEW_TIME& newTime, const NEW_PARENT& newParent) noexcept {
        for (uint32_t sample = 0; sample < NumberOfSamples; sample++) {
            const bool update = times[sample] > newTime(sample);
            parents[sample] = update ? newParent(sample) : parents[sample];
            times[sample] = update ? newTime(sample) : times[sample];
        }
    }

    inline static double delayProbability(const double time, const double maxDelay) noexcept {
        if (time < 0) return 0.0;
        if (time >= maxDelay) return 1.0;
        return (31.0/30.0) - ((11.0/30.0) * (maxDelay / ((10.0 * time) + maxDelay)));
    }

    struct ProfileEntry {
        ProfileEntry(const TimeType perceivedArrivalTime = NotReachable) :
            departureTime(NotReachable) {
            perceivedArrivalTimes.fill(perceivedArrivalTime);
            connections.fill(noConnection);
        }
        int32_t departureTime;
        TimeSample perceivedArrivalTimes;
        ParentSample connections;
    };

    struct VertexLabel {
        VertexLabel() :
            index(0) {
        }
        inline void setIndexFront() noexcept {
            index = 0;
        }
        inline void setIndexBack() noexcept {
            index = profile.size() - 1;
        }
        inline void findIndex(int32_t time) noexcept {
// warning("findIndex 1 ", index, ", ", profile.size());
            while ((index > 0) && (profile[index - 1].departureTime >= time)) {
// warning("findIndex 2");
                updateLabel(profile[index - 1].perceivedArrivalTimes, profile[index - 1].connections, [&](const uint32_t sample){
                    return profile[index].perceivedArrivalTimes[sample];
                }, [&](const uint32_t sample){
                    return profile[index].connections[sample];
                });
                index--;
// warning("findIndex 3");
            }
// warning("findIndex 4");
        }
        inline const ProfileEntry& entry() const noexcept {
            return profile[index];
        }
        inline TimeSample evaluatePAT(const int32_t time, const int32_t maxDelay, const TimeSample& transferCosts) const noexcept {
            TimeSample result;
            if (maxDelay <= 0) {
                for (uint32_t sample = 0; sample < NumberOfSamples; sample++) {
                    result[sample] = profile[index].perceivedArrivalTimes[sample] + transferCosts[sample] - (time * perceivedWaitingFactor[sample]);
                }
            } else {
                result.fill(0);
                double probability = 0.0;
                for (size_t i = index; i < profile.size(); i++) {
                    const double newProbability = delayProbability(profile[i].departureTime - time, maxDelay);
                    for (uint32_t sample = 0; sample < NumberOfSamples; sample++) {
                        result[sample] = result[sample] + ((newProbability - probability) * (profile[i].perceivedArrivalTimes[sample] - (time * perceivedWaitingFactor[sample])));
                    }
                    probability = newProbability;
                    if (probability >= 1) break;
                }
                if (probability >= 1.0) {
                    for (uint32_t sample = 0; sample < NumberOfSamples; sample++) {
                        result[sample] = result[sample] + transferCosts[sample];
                    }
                } else if (probability > 0.0000001) {
                    for (uint32_t sample = 0; sample < NumberOfSamples; sample++) {
                        result[sample] = (result[sample] * (1 / probability)) + transferCosts[sample];
                    }
                } else {
                    result.fill(NotReachable);
                }
            }
            return result;
        }
        inline ParentSample evaluate(const int32_t minTime, const int32_t maxTime) {
            while ((index > 0) && (profile[index - 1].departureTime >= minTime)) {
                index--;
            }
            while ((index < profile.size()) && (profile[index].departureTime < minTime)) {
                index++;
            }
            ProfileEntry result;
            for (uint32_t i = index; i < profile.size(); i++) {
                const int32_t time = std::min(profile[i].departureTime, maxTime);
                updateLabel(result.perceivedArrivalTimes, result.connections, [&](const uint32_t sample){
                    return profile[i].perceivedArrivalTimes[sample] - (time * perceivedWaitingFactor[sample]);
                }, [&](const uint32_t sample){
                    return profile[i].connections[sample];
                });
                if (profile[i].departureTime >= maxTime) break;
            }
            return result.connections;
        }
        inline ParentSample evaluate(const int32_t minTime, const int32_t maxTime, const int32_t maxAdaptationTime, const double adaptationCost) {
            while ((index > 0) && (profile[index - 1].departureTime >= minTime - maxAdaptationTime)) {
                index--;
            }
            while ((index < profile.size()) && (profile[index].departureTime < minTime - maxAdaptationTime)) {
                index++;
            }
            ProfileEntry result;
            for (uint32_t i = index; i < profile.size(); i++) {
                const int32_t time = std::min(profile[i].departureTime, maxTime + maxAdaptationTime);
                const int32_t perceivedAdaptationTime = std::max(std::max(0, minTime - time), time - maxTime) * adaptationCost;
                updateLabel(result.perceivedArrivalTimes, result.connections, [&](const uint32_t sample){
                    return profile[i].perceivedArrivalTimes[sample] - (time * perceivedWaitingFactor[sample]) + perceivedAdaptationTime;
                }, [&](const uint32_t sample){
                    return profile[i].connections[sample];
                });
                if (profile[i].departureTime >= maxTime + maxAdaptationTime) break;
            }
            return result.connections;
        }
        uint32_t index;
        std::vector<ProfileEntry> profile;
        FactorSample perceivedWaitingFactor;
    };

    struct ConnectionLabel {
        ConnectionLabel() {
            parents.fill(noConnection);
        }
        ParentSample parents;
    };

    struct TripLabel {
        TripLabel() :
            parent(noConnection) {
            perceivedArrivalTimes.fill(NotReachable);
        }
        ConnectionId parent;
        TimeSample perceivedArrivalTimes;
    };

public:
    SimpleSamplePATComputation(const CSA::Data& data, const CSA::TransferGraph& reverseGraph, const Settings& settings, const Profiler& profiler = Profiler()) :
        data(data),
        reverseGraph(reverseGraph),
        settings(settings),
        zeroEdge(Edge(reverseGraph.numEdges())),
        connectionLabels(data.numberOfConnections()),
        vertexLabels(reverseGraph.numVertices()),
        tripLabels(data.numberOfTrips()),
        finalFootpaths(reverseGraph.numVertices(), noEdge),
        targetVertex(noVertex),
        profiler(profiler) {
        std::vector<std::vector<std::pair<int32_t, ConnectionId>>> tempProfiles(reverseGraph.numVertices());
        for (const ConnectionId i : data.connectionIds()) {
            const CSA::Connection& connection = data.connections[i];
            const StopId departureStop = connection.departureStopId;
            const int32_t departureTime = connection.departureTime - data.minTransferTime(departureStop);
            tempProfiles[departureStop].emplace_back(departureTime, i);
            for (const Edge edge : reverseGraph.edgesFrom(departureStop)) {
                tempProfiles[reverseGraph.get(ToVertex, edge)].emplace_back(departureTime - reverseGraph.get(TravelTime, edge), i);
            }
        }
        DynamicConnectionGraph tempGraph;
        tempGraph.addVertices(data.numberOfConnections());
        for (const Vertex to : reverseGraph.vertices()) {
            sort(tempProfiles[to]);
            // if (tempProfiles[to].empty()) {
            //     warning("empty ", to);
            // }
            vertexLabels[to].profile.resize(std::max<size_t>(1, tempProfiles[to].size()));
            for (uint32_t index = 0; index < tempProfiles[to].size(); index++) {
                const Vertex from = Vertex(tempProfiles[to][index].second);
                const CSA::Connection& connection = data.connections[from];
                //const int32_t travelTime = connection.departureTime - tempProfiles[to][index].first - data.minTransferTime(connection.departureStopId);
                const Edge originalEdge = (connection.departureStopId == to) ? (zeroEdge) : (reverseGraph.findEdge(connection.departureStopId, to));
                Ensure(originalEdge != noEdge, "There is no edge from " << connection.departureStopId << " to " << to << "!");
                tempGraph.addEdge(from, to).set(OriginalEdge, originalEdge).set(ProfileIndex, index);
                vertexLabels[to].profile[index].departureTime = tempProfiles[to][index].first;
            }
        }
        Graph::move(std::move(tempGraph), connectionGraph);
        // exit(0);
    }

    inline void sampleRandomUtilities() {
        if (settings.randomUtilityDistribution == 0) {
            LogitNormal networkDrivingCostsDistribution(settings.networkDrivingCostsMin, settings.networkDrivingCostsMax, settings.networkDrivingCostsA, settings.networkDrivingCostsB, settings.randomSeed);
            LogitNormal networkWaitingCostsDistribution(settings.networkWaitingCostsMin, settings.networkWaitingCostsMax, settings.networkWaitingCostsA, settings.networkWaitingCostsB, settings.randomSeed);
            LogitNormal networkTransferCostsDistribution(settings.networkTransferCostsMin, settings.networkTransferCostsMax, settings.networkTransferCostsA, settings.networkTransferCostsB, settings.randomSeed);
            LogitNormal networkWalkingCostsDistribution(settings.networkWalkingCostsMin, settings.networkWalkingCostsMax, settings.networkWalkingCostsA, settings.networkWalkingCostsB, settings.randomSeed);
            LogitNormal passengerDrivingCostsDistribution(settings.passengerDrivingCostsMin, settings.passengerDrivingCostsMax, settings.passengerDrivingCostsA, settings.passengerDrivingCostsB, settings.randomSeed);
            LogitNormal passengerWaitingCostsDistribution(settings.passengerWaitingCostsMin, settings.passengerWaitingCostsMax, settings.passengerWaitingCostsA, settings.passengerWaitingCostsB, settings.randomSeed);
            LogitNormal passengerTransferCostsDistribution(settings.passengerTransferCostsMin, settings.passengerTransferCostsMax, settings.passengerTransferCostsA, settings.passengerTransferCostsB, settings.randomSeed);
            LogitNormal passengerWalkingCostsDistribution(settings.passengerWalkingCostsMin, settings.passengerWalkingCostsMax, settings.passengerWalkingCostsA, settings.passengerWalkingCostsB, settings.randomSeed);
            sampleRandomUtilities(networkDrivingCostsDistribution, networkWaitingCostsDistribution, networkTransferCostsDistribution, networkWalkingCostsDistribution, passengerDrivingCostsDistribution, passengerWaitingCostsDistribution, passengerTransferCostsDistribution, passengerWalkingCostsDistribution);
        } else if (settings.randomUtilityDistribution == 1) {
            PERT networkDrivingCostsDistribution(settings.networkDrivingCostsMin, settings.networkDrivingCostsMax, settings.networkDrivingCostsA, settings.networkDrivingCostsB, settings.randomSeed);
            PERT networkWaitingCostsDistribution(settings.networkWaitingCostsMin, settings.networkWaitingCostsMax, settings.networkWaitingCostsA, settings.networkWaitingCostsB, settings.randomSeed);
            PERT networkTransferCostsDistribution(settings.networkTransferCostsMin, settings.networkTransferCostsMax, settings.networkTransferCostsA, settings.networkTransferCostsB, settings.randomSeed);
            PERT networkWalkingCostsDistribution(settings.networkWalkingCostsMin, settings.networkWalkingCostsMax, settings.networkWalkingCostsA, settings.networkWalkingCostsB, settings.randomSeed);
            PERT passengerDrivingCostsDistribution(settings.passengerDrivingCostsMin, settings.passengerDrivingCostsMax, settings.passengerDrivingCostsA, settings.passengerDrivingCostsB, settings.randomSeed);
            PERT passengerWaitingCostsDistribution(settings.passengerWaitingCostsMin, settings.passengerWaitingCostsMax, settings.passengerWaitingCostsA, settings.passengerWaitingCostsB, settings.randomSeed);
            PERT passengerTransferCostsDistribution(settings.passengerTransferCostsMin, settings.passengerTransferCostsMax, settings.passengerTransferCostsA, settings.passengerTransferCostsB, settings.randomSeed);
            PERT passengerWalkingCostsDistribution(settings.passengerWalkingCostsMin, settings.passengerWalkingCostsMax, settings.passengerWalkingCostsA, settings.passengerWalkingCostsB, settings.randomSeed);
            sampleRandomUtilities(networkDrivingCostsDistribution, networkWaitingCostsDistribution, networkTransferCostsDistribution, networkWalkingCostsDistribution, passengerDrivingCostsDistribution, passengerWaitingCostsDistribution, passengerTransferCostsDistribution, passengerWalkingCostsDistribution);
        } else if (settings.randomUtilityDistribution == 2) {
            Kumaraswamy networkDrivingCostsDistribution(settings.networkDrivingCostsMin, settings.networkDrivingCostsMax, settings.networkDrivingCostsA, settings.networkDrivingCostsB, settings.randomSeed);
            Kumaraswamy networkWaitingCostsDistribution(settings.networkWaitingCostsMin, settings.networkWaitingCostsMax, settings.networkWaitingCostsA, settings.networkWaitingCostsB, settings.randomSeed);
            Kumaraswamy networkTransferCostsDistribution(settings.networkTransferCostsMin, settings.networkTransferCostsMax, settings.networkTransferCostsA, settings.networkTransferCostsB, settings.randomSeed);
            Kumaraswamy networkWalkingCostsDistribution(settings.networkWalkingCostsMin, settings.networkWalkingCostsMax, settings.networkWalkingCostsA, settings.networkWalkingCostsB, settings.randomSeed);
            Kumaraswamy passengerDrivingCostsDistribution(settings.passengerDrivingCostsMin, settings.passengerDrivingCostsMax, settings.passengerDrivingCostsA, settings.passengerDrivingCostsB, settings.randomSeed);
            Kumaraswamy passengerWaitingCostsDistribution(settings.passengerWaitingCostsMin, settings.passengerWaitingCostsMax, settings.passengerWaitingCostsA, settings.passengerWaitingCostsB, settings.randomSeed);
            Kumaraswamy passengerTransferCostsDistribution(settings.passengerTransferCostsMin, settings.passengerTransferCostsMax, settings.passengerTransferCostsA, settings.passengerTransferCostsB, settings.randomSeed);
            Kumaraswamy passengerWalkingCostsDistribution(settings.passengerWalkingCostsMin, settings.passengerWalkingCostsMax, settings.passengerWalkingCostsA, settings.passengerWalkingCostsB, settings.randomSeed);
            sampleRandomUtilities(networkDrivingCostsDistribution, networkWaitingCostsDistribution, networkTransferCostsDistribution, networkWalkingCostsDistribution, passengerDrivingCostsDistribution, passengerWaitingCostsDistribution, passengerTransferCostsDistribution, passengerWalkingCostsDistribution);
        } else if (settings.randomUtilityDistribution == 3) {
            Beta networkDrivingCostsDistribution(settings.networkDrivingCostsA, settings.networkDrivingCostsB, settings.randomSeed);
            Beta networkWaitingCostsDistribution(settings.networkWaitingCostsA, settings.networkWaitingCostsB, settings.randomSeed);
            Beta networkTransferCostsDistribution(settings.networkTransferCostsA, settings.networkTransferCostsB, settings.randomSeed);
            Beta networkWalkingCostsDistribution(settings.networkWalkingCostsA, settings.networkWalkingCostsB, settings.randomSeed);
            Beta passengerDrivingCostsDistribution(settings.passengerDrivingCostsA, settings.passengerDrivingCostsB, settings.randomSeed);
            Beta passengerWaitingCostsDistribution(settings.passengerWaitingCostsA, settings.passengerWaitingCostsB, settings.randomSeed);
            Beta passengerTransferCostsDistribution(settings.passengerTransferCostsA, settings.passengerTransferCostsB, settings.randomSeed);
            Beta passengerWalkingCostsDistribution(settings.passengerWalkingCostsA, settings.passengerWalkingCostsB, settings.randomSeed);
            sampleRandomUtilities(networkDrivingCostsDistribution, networkWaitingCostsDistribution, networkTransferCostsDistribution, networkWalkingCostsDistribution, passengerDrivingCostsDistribution, passengerWaitingCostsDistribution, passengerTransferCostsDistribution, passengerWalkingCostsDistribution);
        }
    }

    inline void run(const Vertex target, const int minDepartureTime = -NotReachable) noexcept {
        profiler.startInitialization();
        Vector::fill(tripLabels);
        for (const Vertex vertex : reverseGraph.vertices()) {
            vertexLabels[vertex].setIndexBack();
        }
        clearFinalFootpaths();
        targetVertex = target;
        findFinalFootpaths();
        profiler.doneInitialization();
        for (ConnectionId i = ConnectionId(data.numberOfConnections() - 1); i < data.numberOfConnections(); i--) {
            profiler.scanConnection(i);
            const CSA::Connection& connection = data.connections[i];
            if (connection.departureTime < minDepartureTime) break;

            vertexLabels[connection.arrivalStopId].findIndex(connection.arrivalTime);
            TimeSample perceivedArrivalTimes;
            perceivedArrivalTimes.fill(NotReachable);
            connectionLabels[i].parents.fill(noConnection);

            // if (connection.arrivalStopId == target) { // The connection arrives at the target
            //     updateLabel(perceivedArrivalTimes, connectionLabels[i].parents, [&](const uint32_t){
            //         return connection.arrivalTime;
            //     }, [](const uint32_t){
            //         return noConnection;
            //     });
            // }
            if (finalFootpaths[connection.arrivalStopId] != noEdge) { // It is possible to walk from the connections arrival stop to the target
                updateLabel(perceivedArrivalTimes, connectionLabels[i].parents, [&](const uint32_t sample){
                    return connection.arrivalTime + perceivedTravelTimesOfEdge[finalFootpaths[connection.arrivalStopId]][sample];
                }, [](const uint32_t){
                    return noConnection;
                });
            }
            if (tripLabels[connection.tripId].parent != noConnection) { // It is possible to use another connection of the same trip
                updateLabel(perceivedArrivalTimes, connectionLabels[i].parents, [&](const uint32_t sample){
                    return tripLabels[connection.tripId].perceivedArrivalTimes[sample];
                }, [&](const uint32_t){
                    return tripLabels[connection.tripId].parent;
                });
            }
            if (vertexLabels[connection.arrivalStopId].entry().departureTime >= connection.arrivalTime) { // It is possible to change trips at the connections arrival stop
                const StopId arrivalStop = connection.arrivalStopId;
                const TimeSample newPerceivedArrivalTimes = vertexLabels[arrivalStop].evaluatePAT(connection.arrivalTime, settings.maxDelay, perceivedTransferCostsOfStop[arrivalStop]);
                updateLabel(perceivedArrivalTimes, connectionLabels[i].parents, [&](const uint32_t sample){
                    return newPerceivedArrivalTimes[sample];
                }, [&](const uint32_t sample){
                    return vertexLabels[arrivalStop].entry().connections[sample];
                });
            }

            tripLabels[connection.tripId].parent = i;
            for (uint32_t sample = 0; sample < NumberOfSamples; sample++) {
                if (perceivedArrivalTimes[sample] > NotReachable) continue;
                perceivedArrivalTimes[sample] += perceivedTravelTimesOfConnection[i][sample];
                tripLabels[connection.tripId].perceivedArrivalTimes[sample] = perceivedArrivalTimes[sample];
                perceivedArrivalTimes[sample] += data.minTransferTime(connection.departureStopId) * vertexLabels[connection.departureStopId].perceivedWaitingFactor[sample];
            }

            for (const Edge edge : connectionGraph.edgesFrom(Vertex(i))) {
                const Vertex vertex = connectionGraph.get(ToVertex, edge);
                const Edge originalEdge = Edge(connectionGraph.get(OriginalEdge, edge));
                const uint32_t index = connectionGraph.get(ProfileIndex, edge);
                for (uint32_t sample = 0; sample < NumberOfSamples; sample++) {
                    vertexLabels[vertex].profile[index].perceivedArrivalTimes[sample] = perceivedArrivalTimes[sample] + perceivedTravelTimesOfEdge[originalEdge][sample] + (vertexLabels[vertex].profile[index].departureTime * vertexLabels[vertex].perceivedWaitingFactor[sample]);
                    vertexLabels[vertex].profile[index].connections[sample] = i;
                }
            }
        }
        for (const Vertex vertex : reverseGraph.vertices()) {
            vertexLabels[vertex].findIndex(minDepartureTime);
        }
    }

    inline const VertexLabel& getProfile(const Vertex vertex) const noexcept {
        return vertexLabels[vertex];
    }

    inline VertexLabel& getProfile(const Vertex vertex) noexcept {
        return vertexLabels[vertex];
    }

    inline const ConnectionLabel& getConnectionLabel(const ConnectionId i) const noexcept {
        return connectionLabels[i];
    }

    inline ConnectionLabel& getConnectionLabel(const ConnectionId i) noexcept {
        return connectionLabels[i];
    }

private:
    template<typename TYPE, typename RANDOM>
    inline TYPE getSample(const RANDOM& random) const noexcept {
        TYPE result;
        for (uint32_t sample = 0; sample < NumberOfSamples; sample++) {
            result[sample] = random();
        }
        return result;
    }

    template<typename DISTRIBUTION>
    inline void sampleRandomUtilities(const DISTRIBUTION& networkDrivingCostsDistribution, const DISTRIBUTION& networkWaitingCostsDistribution, const DISTRIBUTION& networkTransferCostsDistribution, const DISTRIBUTION& networkWalkingCostsDistribution, const DISTRIBUTION& passengerDrivingCostsDistribution, const DISTRIBUTION& passengerWaitingCostsDistribution, const DISTRIBUTION& passengerTransferCostsDistribution, const DISTRIBUTION& passengerWalkingCostsDistribution) {
        FactorSample passengerWaitingFactor = getSample<FactorSample>(passengerWaitingCostsDistribution);
        for (const Vertex vertex : reverseGraph.vertices()) {
            for (uint32_t sample = 0; sample < NumberOfSamples; sample++) {
                vertexLabels[vertex].perceivedWaitingFactor[sample] = networkWaitingCostsDistribution() + passengerWaitingFactor[sample];
            }
        }
        TimeSample passengerDrivingCosts = getSample<TimeSample>(passengerDrivingCostsDistribution);
        perceivedTravelTimesOfConnection.resize(data.numberOfConnections());
        for (const ConnectionId i : data.connectionIds()) {
            for (uint32_t sample = 0; sample < NumberOfSamples; sample++) {
                perceivedTravelTimesOfConnection[i][sample] = data.connections[i].travelTime() * (networkDrivingCostsDistribution() + passengerDrivingCosts[sample]);
            }
        }
        TimeSample passengerTransferCosts = getSample<TimeSample>(passengerTransferCostsDistribution);
        perceivedTransferCostsOfStop.resize(data.numberOfStops());
        for (const StopId stop : data.stops()) {
            for (uint32_t sample = 0; sample < NumberOfSamples; sample++) {
                perceivedTransferCostsOfStop[stop][sample] = networkTransferCostsDistribution() + passengerTransferCosts[sample];
            }
        }
        TimeSample passengerWalkingCosts = getSample<TimeSample>(passengerWalkingCostsDistribution);
        perceivedTravelTimesOfEdge.resize(reverseGraph.numEdges() + 1);
        for (const Edge edge : reverseGraph.edges()) {
            for (uint32_t sample = 0; sample < NumberOfSamples; sample++) {
                perceivedTravelTimesOfEdge[edge][sample] = reverseGraph.get(TravelTime, edge) * (networkWalkingCostsDistribution() + passengerWalkingCosts[sample]);
            }
        }
        perceivedTravelTimesOfEdge[zeroEdge].fill(0);
    }

    inline void findFinalFootpaths() noexcept {
        if (!reverseGraph.isVertex(targetVertex)) return;
        finalFootpaths[targetVertex] = zeroEdge;
        for (const Edge edge : reverseGraph.edgesFrom(targetVertex)) {
            profiler.relaxEdge(edge);
            const Vertex stop = reverseGraph.get(ToVertex, edge);
            if (!data.isStop(stop)) continue;
            finalFootpaths[stop] = edge;
        }
    }

    inline void clearFinalFootpaths() noexcept {
        if (!reverseGraph.isVertex(targetVertex)) return;
        finalFootpaths[targetVertex] = noEdge;
        for (const Edge edge : reverseGraph.edgesFrom(targetVertex)) {
            profiler.relaxEdge(edge);
            const Vertex stop = reverseGraph.get(ToVertex, edge);
            if (!data.isStop(stop)) continue;
            finalFootpaths[stop] = noEdge;
        }
    }

private:
    const CSA::Data& data;
    const CSA::TransferGraph& reverseGraph;
    const Settings& settings;

    const Edge zeroEdge;

    ConnectionGraph connectionGraph;

    std::vector<TimeSample> perceivedTravelTimesOfConnection;
    std::vector<TimeSample> perceivedTransferCostsOfStop;
    std::vector<TimeSample> perceivedTravelTimesOfEdge;

    std::vector<ConnectionLabel> connectionLabels;
    std::vector<VertexLabel> vertexLabels;
    std::vector<TripLabel> tripLabels;
    std::vector<Edge> finalFootpaths;
    Vertex targetVertex;

    Profiler profiler;

};

}
