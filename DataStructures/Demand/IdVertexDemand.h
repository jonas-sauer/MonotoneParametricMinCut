#pragma once

#include <iostream>
#include <sstream>
#include <vector>
#include <string>

#include "Passenger.h"
#include "AccumulatedVertexDemand.h"

#include "../CSA/Data.h"

#include "../../Helpers/Assert.h"
#include "../../Helpers/IO/Serialization.h"
#include "../../Helpers/IO/ParserCSV.h"
#include "../../Helpers/String/String.h"

class IdVertexDemand {

public:
    static const std::string CSV_HEADER;

    struct Entry {
        Entry(const DestinationSpecificPassengerList ids = DestinationSpecificPassengerList(), const int departureTime = 0, const Vertex originVertex = noVertex, const Vertex destinationVertex = noVertex) :
            ids(ids),
            departureTime(departureTime),
            originVertex(originVertex),
            destinationVertex(destinationVertex) {
        }
        Entry(const AccumulatedVertexDemand::Entry& d) :
            ids(),
            departureTime(d.earliestDepartureTime),
            originVertex(d.originVertex),
            destinationVertex(d.destinationVertex) {
        }
        Entry(IO::Deserialization& deserialize) {
            this->deserialize(deserialize);
        }
        friend std::ostream& operator<<(std::ostream& out, const Entry& e) {
            return out << "IdVertexDemand::Entry{" << e.ids.size() << ", " << e.departureTime  << ", " << e.originVertex  << ", " << e.destinationVertex << "}";
        }
        inline void serialize(IO::Serialization& serialize) const noexcept {
            serialize(ids, departureTime, originVertex, destinationVertex);
        }
        inline void deserialize(IO::Deserialization& deserialize) noexcept {
            deserialize(ids, departureTime, originVertex, destinationVertex);
        }
        inline std::ostream& toCSV(std::ostream& out) const {
            for (const int id : ids) {
                out << id << "," << departureTime << "," << originVertex << "," << destinationVertex;
            }
            return out;
        }
        inline std::string toCSV() const {
            std::stringstream ss;
            toCSV(ss);
            return ss.str();
        }
        inline long long byteSize() const noexcept {
            return Vector::byteSize(ids) + sizeof(IdVertexDemand);
        }
        DestinationSpecificPassengerList ids{};
        int departureTime{-1};
        Vertex originVertex{noVertex};
        Vertex destinationVertex{noVertex};
    };

private:
    IdVertexDemand() : numIds(0), numberOfPassengers(0) {};

public:
    inline static IdVertexDemand FromBinary(const std::string& filename) noexcept {
        IdVertexDemand result;
        result.deserialize(filename);
        return result;
    }

    inline static IdVertexDemand FromAccumulatedVertexDemand(const AccumulatedVertexDemand& demand, const int multiplier = 10, const int timeStep = 300, const bool includeIntervalBorder = true) noexcept {
        IdVertexDemand result;
        std::vector<DestinationSpecificPassengerId> idByDestination;
        for (const AccumulatedVertexDemand::Entry& initialEntry : demand.entries) {
            std::vector<AccumulatedVertexDemand::Entry> entries;
            u_int32_t timeOffset = 0;
            for (size_t i = 0; i < (initialEntry.numberOfPassengers * multiplier); i++) {
                if (entries.size() <= timeOffset) entries.resize(timeOffset + 1, AccumulatedVertexDemand::Entry(initialEntry, (timeOffset * timeStep)));
                entries[timeOffset].numberOfPassengers++;
                timeOffset++;
                if (includeIntervalBorder) {
                    if ((initialEntry.earliestDepartureTime + int(timeOffset * timeStep)) > initialEntry.latestDepartureTime) timeOffset = 0;
                } else {
                    if ((initialEntry.earliestDepartureTime + int(timeOffset * timeStep)) >= initialEntry.latestDepartureTime) timeOffset = 0;
                }
            }
            if (initialEntry.destinationVertex >= idByDestination.size()) idByDestination.resize(initialEntry.destinationVertex + 1, 0);
            DestinationSpecificPassengerId& id = idByDestination[initialEntry.destinationVertex];
            for (const AccumulatedVertexDemand::Entry& processedEntry : entries) {
                if (processedEntry.numberOfPassengers == 0) continue;
                result.entries.emplace_back(processedEntry);
                const size_t idOffset = id + processedEntry.numberOfPassengers;
                while (id < idOffset) {
                    result.entries.back().ids.emplace_back(id);
                    result.numberOfPassengers++;
                    id++;
                    result.numIds = std::max(result.numIds, id);
                }
            }
        }
        result.passengerMultiplier = multiplier;
        return result;
    }

public:
    friend std::ostream& operator<<(std::ostream& out, const IdVertexDemand& d) {
        return out << "IdVertexDemand{" << d.entries.size()  << ", " << d.numIds  << ", " << d.numberOfPassengers  << ", " << d.passengerMultiplier << "}";
    }

    inline void serialize(const std::string& fileName) const noexcept {
        IO::serialize(fileName, entries, numIds, numberOfPassengers, passengerMultiplier);
    }

    inline void deserialize(const std::string& fileName) noexcept {
        IO::deserialize(fileName, entries, numIds, numberOfPassengers, passengerMultiplier);
    }

    inline std::ostream& toCSV(std::ostream& out) const {
        out << CSV_HEADER << "\n";
        for (const Entry entry : entries) {
            entry.toCSV(out);
            out << "\n";
        }
        return out;
    }

    inline void toCSV(const std::string& fileName) const {
        std::ofstream os(fileName);
        AssertMsg(os, "cannot open file: " << fileName);
        AssertMsg(os.is_open(), "cannot open file: " << fileName);
        toCSV(os);
    }

    inline long long byteSize() const noexcept {
        long long result = sizeof(IdVertexDemand);
        for (const Entry& e : entries) {
            result += e.byteSize();
        }
        return result;
    }

public:
    std::vector<Entry> entries;
    DestinationSpecificPassengerId numIds;
    size_t numberOfPassengers;
    int passengerMultiplier;

};

const std::string IdVertexDemand::CSV_HEADER = "id,departureTime,origin,destination";
