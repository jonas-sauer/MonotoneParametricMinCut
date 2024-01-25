#pragma once

#include <vector>

#include "../DecisionModels/DistanceLogit.h"

#include "../../DataStructures/Assignment/Settings.h"

#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/CSA/TimeExpandedNetwork.h"

#include "../../Helpers/Types.h"

namespace Assignment {

class ComputeSequentialPATs {

public:
    using NodeData = DecisionModels::NodeData;
    using DecisionModel = DecisionModels::DistanceLogit;

    struct ConnectionLabel {
        NodeData entryNode;
        NodeData exitNode;
        NodeData transferNode;
        PerceivedTime skipPAT{Unreachable};
        PerceivedTime tripPAT{Unreachable};
    };

public:
    ComputeSequentialPATs(const CSA::Data& data, const CSA::TransferGraph& reverseGraph, const CSA::TimeExpandedNetwork& timeExpandedNetwork, const Settings& settings, const DecisionModel& decisionModel) :
        data(data),
        reverseGraph(reverseGraph),
        timeExpandedNetwork(timeExpandedNetwork),
        settings(settings),
        decisionModel(decisionModel),
        connectionLabels(data.numberOfConnections()),
        tripOption(data.numberOfTrips()),
        skipOption(data.numberOfStops()),
        skipWaitingTime(data.numberOfConnections(), INFTY),
        transferDistanceToTarget(data.numberOfStops(), INFTY),
        targetVertex(noVertex) {
        for (const ConnectionId i : data.connectionIds()) {
            const ConnectionId next = timeExpandedNetwork.getNextDepartingConnection(i);
            if (next == noConnection) continue;
            skipWaitingTime[i] = settings.waitingCosts * (data.connections[next].departureTime - data.connections[i].departureTime);
        }
    }

    inline void run(const Vertex target) noexcept {
        initialize(target);
        for (ConnectionId i = ConnectionId(data.numberOfConnections() - 1); i < data.numberOfConnections(); i--) {
            const CSA::Connection& connection = data.connections[i];

            std::vector<NodeData> transferOptions;
            const int targetDistance = transferDistanceToTarget[connection.arrivalStopId];
            if (targetDistance != INFTY) {
                const PerceivedTime perceivedTargetArrivalTime = connection.arrivalTime + (settings.walkingCosts + 1) * targetDistance;
                const int targetArrivalTime = connection.arrivalTime + targetDistance;
                transferOptions.emplace_back(NodeData{perceivedTargetArrivalTime, targetArrivalTime});
            }
            for (const LinkId id : timeExpandedNetwork.getOutgoingLinkIDs(i)) {
                const CSA::Link& link = timeExpandedNetwork.getLink(id);
                const NodeData& headNode = connectionLabels[link.headConnection].entryNode;
                if (headNode.arrivalTime >= INFTY) continue;
                transferOptions.emplace_back(headNode);
                transferOptions.back().expectedPAT += settings.walkingCosts * link.travelTime + settings.waitingCosts * link.waitingTime;
            }
            NodeData& transferNode = connectionLabels[i].transferNode;
            decisionModel.merge(transferNode, transferOptions, connection.arrivalTime);

            NodeData& tripNode = tripOption[connection.tripId];
            connectionLabels[i].tripPAT = tripNode.expectedPAT;
            NodeData& exitNode = connectionLabels[i].exitNode;
            decisionModel.merge(exitNode, tripNode, transferNode, connection.arrivalTime);
            tripNode = exitNode;

            NodeData tripEntryNode = tripNode;
            tripEntryNode.expectedPAT += settings.transferCosts;
            NodeData& skipNode = skipOption[connection.departureStopId];
            if (skipNode.arrivalTime != INFTY) {
                skipNode.expectedPAT += skipWaitingTime[i];
            }
            connectionLabels[i].skipPAT = skipNode.expectedPAT;
            NodeData& entryNode = connectionLabels[i].entryNode;
            decisionModel.merge(entryNode, skipNode, tripEntryNode, connection.departureTime);
            skipNode = entryNode;
         }
    }

    inline const ConnectionLabel& connectionLabel(const ConnectionId i) const noexcept {
        return connectionLabels[i];
    }

    inline PerceivedTime targetPAT(const CSA::Connection& connection) const noexcept {
        const int distance = transferDistanceToTarget[connection.arrivalStopId];
        return (distance < INFTY) ? (connection.arrivalTime + (settings.walkingCosts + 1) * distance) : Unreachable;
    }

private:
    inline void initialize(const Vertex target) noexcept {
        clear();
        if (reverseGraph.isVertex(targetVertex)) cleanUp();
        targetVertex = target;
        for (const Edge edge : reverseGraph.edgesFrom(targetVertex)) {
            const Vertex stop = reverseGraph.get(ToVertex, edge);
            if (!data.isStop(stop)) continue;
            transferDistanceToTarget[stop] = reverseGraph.get(TravelTime, edge);
        }
        if (data.isStop(targetVertex)) transferDistanceToTarget[targetVertex] = 0;
    }

    inline void clear() noexcept {
        Vector::fill(tripOption);
        Vector::fill(skipOption);
    }

    inline void cleanUp() noexcept {
        for (const Edge edge : reverseGraph.edgesFrom(targetVertex)) {
            const Vertex stop = reverseGraph.get(ToVertex, edge);
            if (!data.isStop(stop)) continue;
            transferDistanceToTarget[stop] = INFTY;
        }
        if (data.isStop(targetVertex)) transferDistanceToTarget[targetVertex] = INFTY;
    }

private:
    const CSA::Data& data;
    const CSA::TransferGraph& reverseGraph;
    const CSA::TimeExpandedNetwork& timeExpandedNetwork;
    const Settings& settings;
    const DecisionModel& decisionModel;

    std::vector<ConnectionLabel> connectionLabels;
    std::vector<NodeData> tripOption;
    std::vector<NodeData> skipOption;
    std::vector<int> skipWaitingTime;
    std::vector<int> transferDistanceToTarget;
    Vertex targetVertex;
};

}
