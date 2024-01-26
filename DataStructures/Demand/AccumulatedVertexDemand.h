#pragma once

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <random>
#include <tuple>

#include "../CSA/Data.h"

#include "../../Algorithms/CH/CH.h"

#include "../../Helpers/Assert.h"
#include "../../Helpers/IO/Serialization.h"
#include "../../Helpers/IO/ParserCSV.h"
#include "../../Helpers/String/String.h"

class AccumulatedVertexDemand {

public:
    static const std::string CSV_HEADER;

    struct Entry {
        Entry(const int earliestDepartureTime = -1, const int latestDepartureTime = -1, const Vertex originVertex = noVertex, const Vertex destinationVertex = noVertex, const size_t numberOfPassengers = 0) :
            earliestDepartureTime(earliestDepartureTime),
            latestDepartureTime(latestDepartureTime),
            originVertex(originVertex),
            destinationVertex(destinationVertex),
            numberOfPassengers(numberOfPassengers) {
        }
        Entry(const Entry& e, const int offset) :
            demandIndex(e.demandIndex),
            earliestDepartureTime(e.earliestDepartureTime + offset),
            latestDepartureTime(e.earliestDepartureTime + offset),
            originVertex(e.originVertex),
            destinationVertex(e.destinationVertex),
            numberOfPassengers(0) {
        }
        Entry(const Entry& e, const int earliestDepartureTime, const int latestDepartureTime) :
            demandIndex(e.demandIndex),
            earliestDepartureTime(earliestDepartureTime),
            latestDepartureTime(latestDepartureTime),
            originVertex(e.originVertex),
            destinationVertex(e.destinationVertex),
            numberOfPassengers(1) {
        }
        Entry(IO::Deserialization& deserialize) {
            this->deserialize(deserialize);
        }
        friend std::ostream& operator<<(std::ostream& out, const Entry& e) {
            return out << "AccumulatedVertexDemand::Entry{" << e.demandIndex << ", " << e.earliestDepartureTime  << ", " << e.latestDepartureTime  << ", " << e.originVertex  << ", " << e.destinationVertex << ", " << e.numberOfPassengers  << "}";
        }
        inline void serialize(IO::Serialization& serialize) const noexcept {
            serialize(demandIndex, earliestDepartureTime, latestDepartureTime, originVertex, destinationVertex, numberOfPassengers);
        }
        inline void deserialize(IO::Deserialization& deserialize) noexcept {
            deserialize(demandIndex, earliestDepartureTime, latestDepartureTime, originVertex, destinationVertex, numberOfPassengers);
        }
        inline std::ostream& toCSV(std::ostream& out) const {
            return out << demandIndex << "," << earliestDepartureTime << "," << latestDepartureTime << "," << originVertex.value() << "," << destinationVertex.value() << "," << numberOfPassengers;
        }
        inline std::string toCSV() const {
            std::stringstream ss;
            toCSV(ss);
            return ss.str();
        }
        inline std::tuple<Vertex, Vertex, int, int> toTuple() const noexcept {
            return std::make_tuple(destinationVertex, originVertex, earliestDepartureTime, latestDepartureTime);
        }
        inline bool operator<(const Entry& other) const noexcept {
            return toTuple() < other.toTuple();
        }
        size_t demandIndex{(size_t)-1};
        int earliestDepartureTime{-1};
        int latestDepartureTime{-1};
        Vertex originVertex{noVertex};
        Vertex destinationVertex{noVertex};
        size_t numberOfPassengers{0};
    };

public:
    AccumulatedVertexDemand() : numberOfPassengers(0) {}

    inline static AccumulatedVertexDemand FromBinary(const std::string& filename) noexcept {
        AccumulatedVertexDemand result;
        result.deserialize(filename);
        return result;
    }

    inline static AccumulatedVertexDemand FromStopsCSV(const std::string& filename, const CSA::Data& data) noexcept {
        std::cout << "Reading Demand from CSV file (" << filename << ")..." << std::endl;
        AccumulatedVertexDemand result;
        if (FileSystem::isFile(filename)) {
            Timer timer;
            size_t count = 0;
            int firstDeparture = std::numeric_limits<int>::max();
            int lastDeparture = std::numeric_limits<int>::min();
            IO::CSVReader<4, IO::TrimChars<>, IO::DoubleQuoteEscape<',','"'>> in(filename);
            in.readHeader("dep_time", "dep_stop", "arr_stop", "passenger_count");
            Entry demand(-1, -1, noVertex, noVertex, -1);
            while (in.readRow(demand.earliestDepartureTime, demand.originVertex, demand.destinationVertex, demand.numberOfPassengers)) {
                count++;
                if (demand.numberOfPassengers <= 0) continue;
                if (demand.originVertex == demand.destinationVertex) continue;
                if (!data.isStop(demand.originVertex)) continue;
                if (!data.isStop(demand.destinationVertex)) continue;
                demand.demandIndex = count - 1;
                demand.latestDepartureTime = demand.earliestDepartureTime;
                result.entries.emplace_back(demand);
                result.numberOfPassengers += demand.numberOfPassengers;
                if (firstDeparture > demand.earliestDepartureTime) firstDeparture = demand.earliestDepartureTime;
                if (lastDeparture < demand.latestDepartureTime) lastDeparture = demand.latestDepartureTime;
            }
            std::cout << " firstDeparture: " << String::secToTime(firstDeparture) << std::endl;
            std::cout << " lastDeparture:  " << String::secToTime(lastDeparture) << std::endl;
            std::cout << " numberOfPassengers:  " << String::prettyInt(result.numberOfPassengers) << std::endl;
            std::cout << " done (Using " << String::prettyInt(result.entries.size()) << " of " << String::prettyInt(count) << " entries in " << String::msToString(timer.elapsedMilliseconds()) << ")." << std::endl;
        } else {
            std::cout << " file not found." << std::endl;
        }
        return result;
    }

    inline static AccumulatedVertexDemand FromZoneCSV(const std::string& filename, const CSA::Data& data, const CSA::TransferGraph& reverseGraph, const size_t multiplier = 1) noexcept {
        std::cout << "Reading Demand from CSV file (" << filename << ")..." << std::endl;
        AccumulatedVertexDemand result;
        if (FileSystem::isFile(filename)) {
            Timer timer;
            size_t count = 0;
            int firstDeparture = std::numeric_limits<int>::max();
            int lastDeparture = std::numeric_limits<int>::min();
            IO::CSVReader<5, IO::TrimChars<>, IO::DoubleQuoteEscape<',','"'>> in(filename);
            in.readHeader(alias{"min_dep_time", "MINDEPARTURE[SEC]", "MIN_DEPARTURE[SEC]"}, alias{"max_dep_time", "MAXDEPARTURE[SEC]", "MAX_DEPARTURE[SEC]"}, alias{"dep_zone", "FROMZONENO[-]"}, alias{"arr_zone", "TOZONENO[-]"}, alias{"passenger_count", "DEMAND[-]"});
            Entry demand(-1, -1, noVertex, noVertex, -1);
            double passengerFlow;
            while (in.readRow(demand.earliestDepartureTime, demand.latestDepartureTime, demand.originVertex, demand.destinationVertex, passengerFlow)) {
                demand.numberOfPassengers = passengerFlow * multiplier;
                count++;
                if (demand.numberOfPassengers <= 0) continue;
                if (demand.latestDepartureTime < demand.earliestDepartureTime) continue;
                if (demand.originVertex == demand.destinationVertex) continue;
                demand.originVertex += data.numberOfStops();
                if (!data.transferGraph.isVertex(demand.originVertex)) continue;
                if (data.transferGraph.outDegree(demand.originVertex) == 0) continue;
                demand.destinationVertex += data.numberOfStops();
                if (!data.transferGraph.isVertex(demand.destinationVertex)) continue;
                if (reverseGraph.outDegree(demand.destinationVertex) == 0) continue;
                demand.demandIndex = result.entries.size();
                result.entries.emplace_back(demand);
                result.numberOfPassengers += demand.numberOfPassengers;
                if (firstDeparture > demand.earliestDepartureTime) firstDeparture = demand.earliestDepartureTime;
                if (lastDeparture < demand.latestDepartureTime) lastDeparture = demand.latestDepartureTime;
            }
            std::cout << " firstDeparture: " << String::secToTime(firstDeparture) << std::endl;
            std::cout << " lastDeparture:  " << String::secToTime(lastDeparture) << std::endl;
            std::cout << " numberOfPassengers:  " << String::prettyInt(result.numberOfPassengers) << std::endl;
            std::cout << " done (Using " << String::prettyInt(result.entries.size()) << " of " << String::prettyInt(count) << " entries in " << String::msToString(timer.elapsedMilliseconds()) << ")." << std::endl;
        } else {
            std::cout << " file not found." << std::endl;
        }
        return result;
    }

    inline static void MakeImpassableZones(const std::string& inputFileName, const std::string& outputFileName) noexcept {
        std::cout << "Reading Demand from CSV file (" << inputFileName << ")..." << std::endl;
        if (FileSystem::isFile(inputFileName)) {
            AccumulatedVertexDemand result;
            Timer timer;
            size_t count = 0;
            int firstDeparture = std::numeric_limits<int>::max();
            int lastDeparture = std::numeric_limits<int>::min();
            IO::CSVReader<5, IO::TrimChars<>, IO::DoubleQuoteEscape<',','"'>> in(inputFileName);
            in.readHeader(alias{"min_dep_time", "MINDEPARTURE[SEC]", "MIN_DEPARTURE[SEC]"}, alias{"max_dep_time", "MAXDEPARTURE[SEC]", "MAX_DEPARTURE[SEC]"}, alias{"dep_zone", "FROMZONENO[-]"}, alias{"arr_zone", "TOZONENO[-]"}, alias{"passenger_count", "DEMAND[-]"});
            Entry demand(-1, -1, noVertex, noVertex, -1);
            while (in.readRow(demand.earliestDepartureTime, demand.latestDepartureTime, demand.originVertex, demand.destinationVertex, demand.numberOfPassengers)) {
                count++;
                demand.originVertex = Vertex(demand.originVertex * 2);
                demand.destinationVertex = Vertex((demand.destinationVertex * 2) + 1);
                demand.demandIndex = count - 1;
                result.entries.emplace_back(demand);
                result.numberOfPassengers += demand.numberOfPassengers;
                if (firstDeparture > demand.earliestDepartureTime) firstDeparture = demand.earliestDepartureTime;
                if (lastDeparture < demand.latestDepartureTime) lastDeparture = demand.latestDepartureTime;
            }
            std::cout << " firstDeparture: " << String::secToTime(firstDeparture) << std::endl;
            std::cout << " lastDeparture:  " << String::secToTime(lastDeparture) << std::endl;
            std::cout << " numberOfPassengers:  " << String::prettyInt(result.numberOfPassengers) << std::endl;
            result.toCSV(outputFileName);
            std::cout << " done (Using " << String::prettyInt(result.entries.size()) << " of " << String::prettyInt(count) << " entries in " << String::msToString(timer.elapsedMilliseconds()) << ")." << std::endl;
        } else {
            std::cout << " file not found." << std::endl;
        }
    }

    inline static AccumulatedVertexDemand ForDestination(const CSA::Data& data, const Vertex destinationVertex, const int minDepartureTime = 0, const int maxDepartureTime = 86400) noexcept {
        std::cout << "Creating Demand for destination vertex " << destinationVertex << "..." << std::endl;
        Timer timer;
        AccumulatedVertexDemand result;
        Entry demand(-1, -1, noVertex, destinationVertex, 1);
        int firstDeparture = std::numeric_limits<int>::max();
        int lastDeparture = std::numeric_limits<int>::min();
        for (const CSA::Connection& connection : data.connections) {
            if (connection.departureTime < minDepartureTime) continue;
            if (connection.departureTime > maxDepartureTime) continue;
            demand.originVertex = connection.departureStopId;
            demand.earliestDepartureTime = connection.departureTime;
            demand.latestDepartureTime = connection.departureTime;
            result.entries.emplace_back(demand);
            result.numberOfPassengers += demand.numberOfPassengers;
            if (firstDeparture > demand.earliestDepartureTime) firstDeparture = demand.earliestDepartureTime;
            if (lastDeparture < demand.latestDepartureTime) lastDeparture = demand.latestDepartureTime;
        }
        std::cout << " firstDeparture: " << String::secToTime(firstDeparture) << std::endl;
        std::cout << " lastDeparture:  " << String::secToTime(lastDeparture) << std::endl;
        std::cout << " numberOfPassengers:  " << String::prettyInt(result.numberOfPassengers) << std::endl;
        std::cout << " done (Using " << String::prettyInt(result.entries.size()) << " entries in " << String::msToString(timer.elapsedMilliseconds()) << ")." << std::endl;
        return result;
    }

    inline static AccumulatedVertexDemand Random(const CSA::Data& data, const size_t size, const bool stopBased, const bool vertexBased, const int minDepartureTime, const int maxDepartureTime, const int minDepartureWindow, const int maxDepartureWindow, const int minGroupSize, const int maxGroupSize) noexcept {
        std::cout << "Creating random Demand..." << std::endl;
        Timer timer;
        std::vector<bool> canBeOrigin(data.transferGraph.numVertices(), false);
        std::vector<bool> canBeDestination(data.transferGraph.numVertices(), false);
        for (const Vertex from : data.transferGraph.vertices()) {
            for (const Edge edge : data.transferGraph.edgesFrom(from)) {
                const Vertex to = data.transferGraph.get(ToVertex, edge);
                canBeOrigin[from] = true;
                canBeDestination[to] = true;
            }
        }

        std::vector<Vertex> origins;
        std::vector<Vertex> destinations;
        for (const Vertex vertex : data.transferGraph.vertices()) {
            if (data.isStop(vertex)) {
                if (!stopBased) continue;
            } else {
                if (!vertexBased) continue;
            }
            if (canBeOrigin[vertex]) {
                origins.emplace_back(vertex);
            }
            if (canBeDestination[vertex]) {
                destinations.emplace_back(vertex);
            }
        }

        std::random_device randomDevice;
        std::mt19937 randomGenerator(randomDevice());
        std::uniform_int_distribution<size_t> origin(0, origins.size() - 1);
        std::uniform_int_distribution<size_t> destination(0, destinations.size() - 1);
        std::uniform_int_distribution<int> departureTime(minDepartureTime, maxDepartureTime);
        std::uniform_int_distribution<int> departureWindow(minDepartureWindow, maxDepartureWindow);
        std::uniform_int_distribution<int> groupSize(minGroupSize, maxGroupSize);

        AccumulatedVertexDemand result;
        Entry demand(-1, -1, noVertex, noVertex, -1);
        int firstDeparture = std::numeric_limits<int>::max();
        int lastDeparture = std::numeric_limits<int>::min();
        while (result.numberOfPassengers < size) {
            demand.earliestDepartureTime = departureTime(randomGenerator);
            demand.latestDepartureTime = demand.earliestDepartureTime + departureWindow(randomGenerator);
            demand.originVertex = origins[origin(randomGenerator)];
            if (!stopBased) demand.originVertex -= data.numberOfStops();
            demand.destinationVertex = destinations[destination(randomGenerator)];
            if (!stopBased) demand.destinationVertex -= data.numberOfStops();
            if (demand.originVertex == demand.destinationVertex) continue;
            demand.numberOfPassengers = groupSize(randomGenerator);
            demand.demandIndex = result.entries.size();
            result.entries.emplace_back(demand);
            result.numberOfPassengers += demand.numberOfPassengers;
            if (firstDeparture > demand.earliestDepartureTime) firstDeparture = demand.earliestDepartureTime;
            if (lastDeparture < demand.latestDepartureTime) lastDeparture = demand.latestDepartureTime;
        }
        std::cout << " firstDeparture: " << String::secToTime(firstDeparture) << std::endl;
        std::cout << " lastDeparture:  " << String::secToTime(lastDeparture) << std::endl;
        std::cout << " numberOfPassengers:  " << String::prettyInt(result.numberOfPassengers) << std::endl;
        std::cout << " done (Using " << String::prettyInt(result.entries.size()) << " entries in " << String::msToString(timer.elapsedMilliseconds()) << ")." << std::endl;
        return result;
    }

public:
    friend std::ostream& operator<<(std::ostream& out, const AccumulatedVertexDemand& d) {
        return out << "AccumulatedVertexDemand{" << d.entries.size()  << ", " << d.numberOfPassengers << "}";
    }

    //Distribute departure times of each demand entry uniformly across departure interval
    //keepIntervals: true = Retain latestDepartureTime for each new entry. false = Create intervals of length zero.
    //includeIntervalBorder: true = Include upper border in departure interval. false = exclude upper border
    inline void discretize(const int timeStep = 300, const bool keepIntervals = false, const bool includeIntervalBorder = true) noexcept {
        if (timeStep < 0) return;
        std::vector<Entry> newEntries;
        for (const Entry& oldEntry : entries) {
            size_t firstIndex = newEntries.size();
            int timeOffset = 0;
            for (size_t i = 0; i < oldEntry.numberOfPassengers; i++) {
                if (firstIndex + timeOffset >= newEntries.size()) {
                    Assert(newEntries.size() == firstIndex + timeOffset, "Index does not match number of items (Index: " << (firstIndex + timeOffset) << ", number of items: " << newEntries.size() << ")!");
                    const int earliestDepartureTime = oldEntry.earliestDepartureTime + (timeOffset * timeStep);
                    const int latestDepartureTime = (keepIntervals) ? (std::min(earliestDepartureTime + timeStep, oldEntry.latestDepartureTime)) : (earliestDepartureTime);
                    newEntries.emplace_back(oldEntry, earliestDepartureTime, latestDepartureTime);
                } else {
                    newEntries[firstIndex + timeOffset].numberOfPassengers++;
                }
                timeOffset++;
                if (includeIntervalBorder) {
                    if (oldEntry.earliestDepartureTime + (timeOffset * timeStep) > oldEntry.latestDepartureTime) timeOffset = 0;
                } else {
                    if (oldEntry.earliestDepartureTime + (timeOffset * timeStep) >= oldEntry.latestDepartureTime) timeOffset = 0;
                }
            }
        }
        entries.swap(newEntries);
    }

    inline void sanitize() noexcept {
        std::cout << "Sanitizing Demand..." << std::endl;
        Timer timer;
        numberOfPassengers = 0;
        int firstDeparture = std::numeric_limits<int>::max();
        int lastDeparture = std::numeric_limits<int>::min();
        for (size_t i = 0; i < entries.size(); i++) {
            entries[i].demandIndex = i;
            numberOfPassengers += entries[i].numberOfPassengers;
            if (firstDeparture > entries[i].earliestDepartureTime) firstDeparture = entries[i].earliestDepartureTime;
            if (lastDeparture < entries[i].latestDepartureTime) lastDeparture = entries[i].latestDepartureTime;
        }
        std::cout << " firstDeparture: " << String::secToTime(firstDeparture) << std::endl;
        std::cout << " lastDeparture:  " << String::secToTime(lastDeparture) << std::endl;
        std::cout << " numberOfPassengers:  " << String::prettyInt(numberOfPassengers) << std::endl;
        std::cout << " done (Using " << String::prettyInt(entries.size()) << " entries in " << String::msToString(timer.elapsedMilliseconds()) << ")." << std::endl;
    }

    inline void toZoneIDs(const CSA::Data& data) noexcept {
        for (Entry& entry : entries) {
            entry.originVertex -= data.numberOfStops();
            entry.destinationVertex -= data.numberOfStops();
        }
    }

    inline void sortByOrigin() noexcept {
        sort(entries, [](const Entry& a, const Entry& b){return a.originVertex < b.originVertex;});
    }

    inline void sortByDestination() noexcept {
        sort(entries, [](const Entry& a, const Entry& b){return a.destinationVertex < b.destinationVertex;});
    }

    inline void lexicographicalSort() noexcept {
        sort(entries);
    }

    inline void serialize(const std::string& fileName) const noexcept {
        IO::serialize(fileName, entries, numberOfPassengers);
    }

    inline void deserialize(const std::string& fileName) noexcept {
        IO::deserialize(fileName, entries, numberOfPassengers);
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
        Assert(os, "cannot open file: " << fileName);
        Assert(os.is_open(), "cannot open file: " << fileName);
        toCSV(os);
    }

public:
    std::vector<Entry> entries;
    size_t numberOfPassengers;

};

const std::string AccumulatedVertexDemand::CSV_HEADER = "index,min_dep_time,max_dep_time,dep_zone,arr_zone,passenger_count";
