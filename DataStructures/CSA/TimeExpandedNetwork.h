#pragma once

#include <vector>

#include "Data.h"
#include "../../Helpers/Ranges/SubRange.h"
#include "../../Helpers/Types.h"

namespace CSA {

struct Link {
    Link(const ConnectionId /*tail*/, const ConnectionId head, const int travelTime, const int waitingTime = 0) :
        //tailConnection(tail),
        headConnection(head),
        travelTime(travelTime),
        waitingTime(waitingTime) {
    }

    //ConnectionId tailConnection;
    ConnectionId headConnection;
    int travelTime;
    int waitingTime;
};

class TimeExpandedNetwork {
public:
    TimeExpandedNetwork(const Data& data, const bool useTransferBufferTimes) :
        firstOutgoingLink(data.numberOfConnections() + 1, LinkId(0)),
        nextDepartingConnection(data.numberOfConnections(), noConnection) {
        std::vector<std::vector<ConnectionId>> departingConnectionsPerStop(data.numberOfStops());
        for (const ConnectionId i : data.connectionIds()) {
            const Connection& connection = data.connections[i];
            const StopId departureStop = connection.departureStopId;
            if (!departingConnectionsPerStop[departureStop].empty()) {
                const ConnectionId previousDepartingConnection = departingConnectionsPerStop[departureStop].back();
                nextDepartingConnection[previousDepartingConnection] = i;
            }
            departingConnectionsPerStop[departureStop].emplace_back(i);
        }

        for (const ConnectionId i : data.connectionIds()) {
            const Connection& connection = data.connections[i];
            firstOutgoingLink[i] = LinkId(links.size());
            for (const Edge edge : data.transferGraph.edgesFrom(connection.arrivalStopId)) {
                const Vertex toStop = data.transferGraph.get(ToVertex, edge);
                if (!data.isStop(toStop)) continue;
                const int bufferTime = useTransferBufferTimes ? data.minTransferTime(StopId(toStop)) : 0;
                const int travelTime = data.transferGraph.get(TravelTime, edge) + bufferTime;
                const int arrivalTime = connection.arrivalTime + travelTime;
                for (size_t j = 0; j < departingConnectionsPerStop[toStop].size(); j++) {
                    const ConnectionId toConnection = departingConnectionsPerStop[toStop][j];
                    const int waitingTime = data.connections[toConnection].departureTime - arrivalTime;
                    if (waitingTime >= 0) {
                        links.emplace_back(i, toConnection, travelTime, waitingTime);
                        break;
                    }
                }
            }
            for (size_t j = 0; j < departingConnectionsPerStop[connection.arrivalStopId].size(); j++) {
                const ConnectionId toConnection = departingConnectionsPerStop[connection.arrivalStopId][j];
                const int bufferTime = data.minTransferTime(connection.arrivalStopId);
                const int arrivalTime = connection.arrivalTime + bufferTime;
                const int waitingTime = data.connections[toConnection].departureTime - arrivalTime;
                if (waitingTime >= 0) {
                    links.emplace_back(i, toConnection, bufferTime, waitingTime);
                    break;
                }
            }
        }
        firstOutgoingLink.back() = LinkId(links.size());
    }

    inline long long byteSize() const noexcept {
        long long result = Vector::byteSize(links);
        result += Vector::byteSize(firstOutgoingLink);
        result += Vector::byteSize(nextDepartingConnection);
        return result;
    }

    inline size_t numberOfLinks() const noexcept {
        return links.size();
    }

    inline Range<LinkId> getOutgoingLinkIDs(const ConnectionId connection) const noexcept {
        return Range<LinkId>(firstOutgoingLink[connection], firstOutgoingLink[connection + 1]);
    }

    inline const Link& getLink(const LinkId id) const noexcept {
        return links[id];
    }

    inline ConnectionId getNextDepartingConnection(const ConnectionId connection) const noexcept {
        return nextDepartingConnection[connection];
    }

private:
    std::vector<Link> links;
    std::vector<LinkId> firstOutgoingLink;
    std::vector<ConnectionId> nextDepartingConnection;
};

class CapacityTimeExpandedNetwork {

public:
    CapacityTimeExpandedNetwork(const Data& data, const bool useTransferBufferTimes) :
        firstOutgoingLink(data.numberOfConnections() + 1, LinkId(0)),
        firstOutgoingFailureLink(data.numberOfConnections(), LinkId(0)),
        firstIncomingLink(data.numberOfConnections() + 1, 0),
        previousDepartingConnection(data.numberOfConnections(), noConnection),
        previousConnectionInTrip(data.numberOfConnections(), noConnection) {
        std::vector<std::vector<ConnectionId>> departingConnectionsPerStop(data.numberOfStops());
        std::vector<ConnectionId> lastConnectionInTrip(data.numberOfTrips(), noConnection);
        for (const ConnectionId i : data.connectionIds()) {
            const Connection& connection = data.connections[i];
            const StopId departureStop = connection.departureStopId;
            if (!departingConnectionsPerStop[departureStop].empty()) {
                previousDepartingConnection[i] = departingConnectionsPerStop[departureStop].back();
            }
            departingConnectionsPerStop[departureStop].emplace_back(i);
            if (lastConnectionInTrip[connection.tripId] != noConnection) {
                previousConnectionInTrip[i] = lastConnectionInTrip[connection.tripId];
            }
            lastConnectionInTrip[connection.tripId] = i;
        }
        std::vector<std::vector<LinkId>> incomingLinksPerConnection(data.numberOfConnections());
        for (const ConnectionId i : data.connectionIds()) {
            const Connection& connection = data.connections[i];
            firstOutgoingLink[i] = LinkId(links.size());
            for (const Edge edge : data.transferGraph.edgesFrom(connection.arrivalStopId)) {
                const Vertex toStop = data.transferGraph.get(ToVertex, edge);
                if (!data.isStop(toStop)) continue;
                const int bufferTime = useTransferBufferTimes ? data.minTransferTime(StopId(toStop)) : 0;
                const int arrivalTime = connection.arrivalTime + data.transferGraph.get(TravelTime, edge) + bufferTime;
                for (size_t j = 0; j < departingConnectionsPerStop[toStop].size(); j++) {
                    const ConnectionId toConnection = departingConnectionsPerStop[toStop][j];
                    if (data.connections[toConnection].departureTime >= arrivalTime) {
                        incomingLinksPerConnection[toConnection].emplace_back(links.size());
                        links.emplace_back(i, toConnection, data.transferGraph.get(TravelTime, edge));
                        break;
                    }
                }
            }
            for (size_t j = 0; j < departingConnectionsPerStop[connection.arrivalStopId].size(); j++) {
                const ConnectionId toConnection = departingConnectionsPerStop[connection.arrivalStopId][j];
                if (data.connections[toConnection].departureTime >= connection.arrivalTime + data.minTransferTime(connection.arrivalStopId)) {
                    incomingLinksPerConnection[toConnection].emplace_back(links.size());
                    links.emplace_back(i, toConnection, 0);
                    break;
                }
            }
            firstOutgoingFailureLink[i] = LinkId(links.size());
            for (const Edge edge : data.transferGraph.edgesFrom(connection.departureStopId)) {
                const Vertex toStop = data.transferGraph.get(ToVertex, edge);
                if (!data.isStop(toStop)) continue;
                const int bufferTime = useTransferBufferTimes ? data.minTransferTime(StopId(toStop)) : 0;
                const int arrivalTime = connection.departureTime + data.transferGraph.get(TravelTime, edge) + bufferTime;
                for (size_t j = 0; j < departingConnectionsPerStop[toStop].size(); j++) {
                    const ConnectionId toConnection = departingConnectionsPerStop[toStop][j];
                    if (data.connections[toConnection].departureTime >= arrivalTime) {
                        incomingLinksPerConnection[toConnection].emplace_back(links.size());
                        links.emplace_back(i, toConnection, data.transferGraph.get(TravelTime, edge));
                        break;
                    }
                }
            }
            for (size_t j = 0; j < departingConnectionsPerStop[connection.departureStopId].size(); j++) {
                const ConnectionId toConnection = departingConnectionsPerStop[connection.departureStopId][j];
                if (toConnection != i && data.connections[toConnection].departureTime >= connection.departureTime + data.minTransferTime(connection.arrivalStopId)) {
                    incomingLinksPerConnection[toConnection].emplace_back(links.size());
                    links.emplace_back(i, toConnection, 0);
                    break;
                }
            }
        }
        firstOutgoingLink.back() = LinkId(links.size());
        incomingLinks.reserve(links.size());
        for (const ConnectionId i : data.connectionIds()) {
            firstIncomingLink[i] = incomingLinks.size();
            incomingLinks += incomingLinksPerConnection[i];
        }
        firstIncomingLink.back() = incomingLinks.size();
    }

    inline size_t numberOfLinks() const noexcept {
        return links.size();
    }

    inline Range<LinkId> getOutgoingLinkIDs(const ConnectionId connection) const noexcept {
        return Range<LinkId>(firstOutgoingLink[connection], firstOutgoingFailureLink[connection]);
    }

    inline Range<LinkId> getOutgoingFailureLinkIDs(const ConnectionId connection) const noexcept {
        return Range<LinkId>(firstOutgoingFailureLink[connection], firstOutgoingLink[connection + 1]);
    }

    inline const Link& getLink(const LinkId id) const noexcept {
        return links[id];
    }

    inline SubRange<std::vector<LinkId>> getIncomingLinkIDs(const ConnectionId connection) const noexcept {
        return SubRange<std::vector<LinkId>>(incomingLinks, firstIncomingLink[connection], firstIncomingLink[connection + 1]);
    }

    inline ConnectionId getPreviousDepartingConnection(const ConnectionId connection) const noexcept {
        return previousDepartingConnection[connection];
    }

    inline ConnectionId getPreviousConnectionInTrip(const ConnectionId connection) const noexcept {
        return previousConnectionInTrip[connection];
    }

private:
    std::vector<Link> links;
    std::vector<LinkId> firstOutgoingLink;
    std::vector<LinkId> firstOutgoingFailureLink;
    std::vector<size_t> firstIncomingLink;
    std::vector<LinkId> incomingLinks;
    std::vector<ConnectionId> previousDepartingConnection;
    std::vector<ConnectionId> previousConnectionInTrip;
};

}
