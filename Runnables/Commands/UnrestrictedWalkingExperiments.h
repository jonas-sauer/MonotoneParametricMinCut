#pragma once

#include <iostream>
#include <algorithm>
#include <random>
#include <vector>
#include <string>

#include "../../Shell/Shell.h"

#include "../../Helpers/Assert.h"
#include "../../Helpers/String/String.h"
#include "../../Helpers/IO/File.h"
#include "../../Helpers/IO/CSVData.h"
#include "../../Helpers/ConfigFile.h"
#include "../../Helpers/Console/CommandLineParser.h"
#include "../../Helpers/Console/Progress.h"
#include "../../Helpers/Vector/Vector.h"

#include "../../Helpers/Timer.h"
#include "../../Helpers/MultiThreading.h"
#include "../../Helpers/IO/CSVData.h"
#include "../../Helpers/IO/ParserCSV.h"

#include "../../DataStructures/Intermediate/Data.h"
#include "../../DataStructures/RAPTOR/Data.h"
#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/Geometry/Functions.h"
#include "../../DataStructures/Geometry/CoordinateTree.h"
#include "../../DataStructures/RAPTOR/Entities/ArrivalLabel.h"

#include "../../Algorithms/Dijkstra/Dijkstra.h"
#include "../../Algorithms/CSA/CSA.h"
#include "../../Algorithms/RAPTOR/DijkstraRAPTOR.h"
#include "../../Algorithms/RAPTOR/RangeRAPTOR/RangeRAPTOR.h"
#include "../../Algorithms/RAPTOR/RAPTOR.h"
#include "../../Algorithms/RAPTOR/Profiler.h"
#include "../../Algorithms/CH/Preprocessing/CHData.h"
#include "../../Algorithms/CH/Query/CHQuery.h"
#include "../../Algorithms/CSA/OldProfileCSA.h"
#include "../../Algorithms/CSA/ULTRACSA.h"
#include "../../Algorithms/RAPTOR/AlternatingRAPTOR.h"
#include "../../Algorithms/RAPTOR/ULTRARAPTOR.h"

#include "../../Visualization/FunctionVisualization.h"
#include "../../Visualization/PDF.h"
#include "../../Visualization/PNG.h"
#include "../../Visualization/Color.h"

using namespace Shell;

class GenerateRankQueries : public ParameterizedCommand {

public:
    GenerateRankQueries(BasicShell& shell) :
        ParameterizedCommand(shell, "generateRankQueries", "Computes random rank queries for an intermediate network and saves them in a file.") {
        addParameter("Intermediate network");
        addParameter("Number of queries");
        addParameter("Output file");
        addParameter("Seed", "42");
    }

    virtual void execute() noexcept {
        const std::string networkFile = getParameter("Intermediate network");
        const Intermediate::Data inter = Intermediate::Data::FromBinary(networkFile);
        inter.printInfo();
        const size_t n = getParameter<size_t>("Number of queries");
        const std::string outputFile = getParameter("Output file");
        const int seed = getParameter<int>("Seed");
        srand(seed);
        IO::CSVData<size_t, IO::TrimChars<>, IO::DoubleQuoteEscape<',','"'>> result({"queryId", "sourceStop", "targetStop", "rank"});
        Dijkstra<Intermediate::TransferGraph> dijkstra(inter.transferGraph, inter.transferGraph[TravelTime]);
        Progress progress(n);
        for (size_t i = 0; i < n;) {
            const StopId source = StopId(rand() % inter.numberOfStops());
            size_t currentIndex = 0;
            size_t nextIndex = 100;
            std::vector<size_t> targets;
            dijkstra.run(source, noVertex, [&](const Vertex u) {
                currentIndex++;
                if (inter.isStop(u) && currentIndex >= nextIndex) {
                    targets.emplace_back(u);
                    nextIndex = nextIndex << 1;
                }
            });
            if (currentIndex > inter.transferGraph.numVertices() / 2) {
                for (size_t j = 0; j < targets.size(); j++) {
                    result.columnData[0].emplace_back(result.numRows());
                    result.columnData[1].emplace_back(source);
                    result.columnData[2].emplace_back(targets[j]);
                    result.columnData[3].emplace_back(j + 7);
                }
                i++;
                progress++;
            }
        }
        result.write(outputFile);
        std::cout << std::endl;
    }
};

class GenerateConnectionBasedRankQueries : public ParameterizedCommand {

public:
    GenerateConnectionBasedRankQueries(BasicShell& shell) :
        ParameterizedCommand(shell, "generateConnectionBasedRankQueries", "Computes random rank queries, weighted by number of connections, for an intermediate network and saves them in a file.") {
        addParameter("Intermediate network");
        addParameter("Number of queries");
        addParameter("Output file");
        addParameter("Seed", "42");
    }

    virtual void execute() noexcept {
        const std::string networkFile = getParameter("Intermediate network");
        const Intermediate::Data inter = Intermediate::Data::FromBinary(networkFile);
        inter.printInfo();
        const size_t n = getParameter<size_t>("Number of queries");
        const std::string outputFile = getParameter("Output file");
        const int seed = getParameter<int>("Seed");
        srand(seed);
        IO::CSVData<size_t, IO::TrimChars<>, IO::DoubleQuoteEscape<',','"'>> result({"queryId", "sourceStop", "targetStop", "rank"});
        Dijkstra<Intermediate::TransferGraph> dijkstra(inter.transferGraph, inter.transferGraph[TravelTime]);
        std::vector<size_t> numberOfStopEventsByStop(inter.numberOfStops(), 0);
        for (const Intermediate::Trip& trip : inter.trips) {
            for (const Intermediate::StopEvent& stopEvent : trip.stopEvents) {
                numberOfStopEventsByStop[stopEvent.stopId]++;
            }
        }
        std::vector<size_t> cumulativeDistribution;
        cumulativeDistribution.emplace_back(numberOfStopEventsByStop[0]);
        for (size_t i = 1; i < numberOfStopEventsByStop.size(); i++) {
            cumulativeDistribution.emplace_back(cumulativeDistribution.back() + numberOfStopEventsByStop[i]);
        }
        Progress progress(n);
        for (size_t i = 0; i < n;) {
            const StopId source = StopId(chooseRandomly(cumulativeDistribution));
            size_t currentIndex = 0;
            size_t minIndex = 107; //2^6.75
            size_t maxIndex = 153; //2^7.25
            int bestStop = -1;
            std::vector<size_t> targets;
            dijkstra.run(source, noVertex, [&](const Vertex u) {
                currentIndex++;
                if (inter.isStop(u) && currentIndex >= minIndex) {
                    if (bestStop == -1 || numberOfStopEventsByStop[u] > numberOfStopEventsByStop[bestStop]) {
                        bestStop = u;
                    }
                    if (currentIndex > maxIndex) {
                        targets.emplace_back(bestStop);
                        bestStop = -1;
                        minIndex = minIndex << 1;
                        maxIndex = maxIndex << 1;
                    }
                }
            });
            if (currentIndex > inter.transferGraph.numVertices() / 2) {
                for (size_t j = 0; j < targets.size(); j++) {
                    result.columnData[0].emplace_back(result.numRows());
                    result.columnData[1].emplace_back(source);
                    result.columnData[2].emplace_back(targets[j]);
                    result.columnData[3].emplace_back(j + 7);
                }
                i++;
                progress++;
            }
        }
        result.write(outputFile);
        std::cout << std::endl;
    }

private:
    inline size_t chooseRandomly(const std::vector<size_t>& cumulativeDistribution) const noexcept {
        size_t randomValue = rand() % cumulativeDistribution.back();
        for (size_t i = 0; i < cumulativeDistribution.size(); i++) {
            if (randomValue < cumulativeDistribution[i]) {
                return i;
            }
        }
        return cumulativeDistribution.size();
    }
};

class GenerateConnectionBasedRandomQueries : public ParameterizedCommand {

public:
    GenerateConnectionBasedRandomQueries(BasicShell& shell) :
        ParameterizedCommand(shell, "generateConnectionBasedRandomQueries", "Computes random queries, weighted by number of connections, for an intermediate network and saves them in a file.") {
        addParameter("Intermediate network");
        addParameter("Number of queries");
        addParameter("Output file");
        addParameter("Seed", "42");
    }

    virtual void execute() noexcept {
        const std::string networkFile = getParameter("Intermediate network");
        const Intermediate::Data inter = Intermediate::Data::FromBinary(networkFile);
        inter.printInfo();
        const size_t n = getParameter<size_t>("Number of queries");
        const std::string outputFile = getParameter("Output file");
        const int seed = getParameter<int>("Seed");
        srand(seed);
        IO::CSVData<size_t, IO::TrimChars<>, IO::DoubleQuoteEscape<',','"'>> result({"queryId", "sourceStop", "targetStop"});
        std::vector<size_t> numberOfStopEventsByStop(inter.numberOfStops(), 0);
        for (const Intermediate::Trip& trip : inter.trips) {
            for (const Intermediate::StopEvent& stopEvent : trip.stopEvents) {
                numberOfStopEventsByStop[stopEvent.stopId]++;
            }
        }
        std::vector<size_t> cumulativeDistribution;
        cumulativeDistribution.emplace_back(numberOfStopEventsByStop[0]);
        for (size_t i = 1; i < numberOfStopEventsByStop.size(); i++) {
            cumulativeDistribution.emplace_back(cumulativeDistribution.back() + numberOfStopEventsByStop[i]);
        }
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const StopId source = StopId(chooseRandomly(cumulativeDistribution));
            const StopId target = StopId(chooseRandomly(cumulativeDistribution));
            result.columnData[0].emplace_back(i);
            result.columnData[1].emplace_back(source);
            result.columnData[2].emplace_back(target);
        }
        result.write(outputFile);
        std::cout << std::endl;
    }

private:
    inline size_t chooseRandomly(const std::vector<size_t>& cumulativeDistribution) const noexcept {
        size_t randomValue = rand() % cumulativeDistribution.back();
        for (size_t i = 0; i < cumulativeDistribution.size(); i++) {
            if (randomValue < cumulativeDistribution[i]) {
                return i;
            }
        }
        return cumulativeDistribution.size();
    }
};

class CompareTravelTimes : public ParameterizedCommand {

private:
    struct Query {
        int queryId{-1};
        StopId source{noStop};
        StopId target{noStop};
        size_t rank{0};
        double restrictedTime{-1};
        double unrestrictedTime{-1};
        size_t restrictedSize{0};
        size_t unrestrictedSize{0};
        friend std::ostream& operator<<(std::ostream& out, const Query& q) noexcept {
            return out << q.queryId << "," << q.restrictedTime << "," << q.unrestrictedTime << "," << q.restrictedSize << "," << q.unrestrictedSize << "\n";
        }
    };

public:
    CompareTravelTimes(BasicShell& shell) :
        ParameterizedCommand(shell, "compareTravelTimes", "Compares travel times on transitive and unrestricted network for the given queries.") {
        addParameter("Queries input file");
        addParameter("Transitive CSA network");
        addParameter("Unrestricted RAPTOR network");
        addParameter("Min rank");
        addParameter("Max rank");
        addParameter("Output directory");
    }

    virtual void execute() noexcept {
        const std::string queriesFile = getParameter("Queries input file");
        const size_t minRank = getParameter<size_t>("Min rank");
        const size_t maxRank = getParameter<size_t>("Max rank");
        std::vector<std::vector<Query>> queries;
        IO::CSVReader<4, IO::TrimChars<>, IO::DoubleQuoteEscape<',','"'>> in(queriesFile);
        in.readHeader("queryId", "sourceStop", "targetStop", "rank");
        Query q;
        size_t numberOfQueries = 0;
        while (in.readRow(q.queryId, q.source, q.target, q.rank)) {
            if (q.rank < minRank || q.rank > maxRank) break;
            if (queries.size() <= q.rank) {
                queries.resize(q.rank + 1);
            }
            queries[q.rank].emplace_back(q);
            numberOfQueries++;
        }

        const std::string csaFile = getParameter("Transitive CSA network");
        CSA::Data csa_data = CSA::Data::FromBinary(csaFile);
        csa_data.sortConnectionsAscendingByDepartureTime();
        CSA::OldProfileCSA<8> csa(csa_data, csa_data.transferGraph);
        csa_data.printInfo();

        const std::string raptorFile = getParameter("Unrestricted RAPTOR network");
        RAPTOR::Data raptor_data = RAPTOR::Data::FromBinary(raptorFile);
        RAPTOR::Data reverse_data = raptor_data.reverseNetwork();
        RAPTOR::AlternatingDijkstraRAPTOR<RAPTOR::DijkstraInitialTransfers, false, true> raptor(raptor_data, reverse_data);
        raptor_data.printInfo();

        const std::string outputDirectory = getParameter("Output directory");
        IO::OFStream resultOverview(outputDirectory + "_overview.csv");
        resultOverview << "queryId,restrictedTime,unrestrictedTime,restrictedSize,unrestrictedSize\n";
        Progress progress(numberOfQueries);
        for (size_t rank = minRank; rank < queries.size(); rank++) {
            if (queries[rank].empty()) continue;
            IO::OFStream resultRank(outputDirectory + "_rank_" + std::to_string(rank) + ".csv");
            resultRank << "queryId,departureTime,restrictedTimeTravelTime,unrestrictedTravelTime\n";
            for (Query& q : queries[rank]) {
                const Geometry::Function::Piecewise restrictedProfile(run(q, csa, q.restrictedTime, q.restrictedSize));
                const Geometry::Function::Piecewise unrestrictedProfile(run(q, raptor, q.unrestrictedTime, q.unrestrictedSize));
                resultOverview << q;
                if (q.restrictedSize * q.unrestrictedSize != 0) {
                    for (int departureTime = 0; departureTime < 24 * 60 * 60; departureTime += 10) {
                        resultRank << q.queryId << "," << departureTime;
                        resultRank << "," << restrictedProfile(departureTime);
                        resultRank << "," << unrestrictedProfile(departureTime);
                        resultRank << "\n";
                    }
                }
                resultRank.flush();
                progress++;
            }
            resultOverview.flush();
        }
    }

private:
    template<typename ALGORITHM>
    inline std::vector<Geometry::Point> run(Query& q, ALGORITHM& algorithm, double& time, size_t& size) const noexcept {
        Timer timer;
        algorithm.run(q.source, q.target, 0, 24 * 60 * 60);
        time = timer.elapsedMilliseconds();
        const RAPTOR::ProfileHandle profile = algorithm.getProfileHandle();
        if (profile.empty()) {
            size = 0;
            return std::vector<Geometry::Point>();
        } else {
            const std::vector<Geometry::Point> function = algorithm.getTravelTimeFunction();
            size = function.size();
            return function;
        }
    }
};

class AnalyzeTravelTimeDifference : public ParameterizedCommand {

private:
    struct Data {
        std::vector<double> restrictedTravelTime;
        std::vector<double> unrestrictedTravelTime;
        std::vector<double> travelTimeDiff;
        int impossible{0};

        inline size_t size() const noexcept {
            return restrictedTravelTime.size();
        }
    };

public:
    AnalyzeTravelTimeDifference(BasicShell& shell) :
        ParameterizedCommand(shell, "analyzeTravelTimeDifference", "Reads travel time profiles and computes some statistics.") {
        addParameter("Input file base path");
        addParameter("Min rank");
        addParameter("Max rank");
        addParameter("Output file base path");
    }

    virtual void execute() noexcept {
        const std::string inputFileBase = getParameter("Input file base path");
        const size_t minRank = getParameter<size_t>("Min rank");
        const size_t maxRank = getParameter<size_t>("Max rank");
        const std::string outputFileBase = getParameter("Output file base path");

        for (size_t rank = minRank; rank <= maxRank; rank++) {
            std::vector<Data> data;
            IO::CSVReader<4, IO::TrimChars<>, IO::DoubleQuoteEscape<',','"'>> in(inputFileBase + std::to_string(rank) + ".csv");
            in.readHeader("queryId", "departureTime", "restrictedTravelTime", "unrestrictedTravelTime");
            int queryId;
            int departureTime;
            double restrictedTravelTime;
            double unrestrictedTravelTime;
            while (in.readRow(queryId, departureTime, restrictedTravelTime, unrestrictedTravelTime)) {
                const size_t chunk = departureTime / TimeResolution;
                if (data.size() <= chunk) data.resize(chunk + 1);
                if (restrictedTravelTime - unrestrictedTravelTime <= 100000) {
                    data[chunk].travelTimeDiff.emplace_back(restrictedTravelTime - unrestrictedTravelTime);
                    data[chunk].restrictedTravelTime.emplace_back(restrictedTravelTime);
                    data[chunk].unrestrictedTravelTime.emplace_back(unrestrictedTravelTime);
                } else {
                    data[chunk].impossible++;
                }
            }

            writeResults(data, outputFileBase + "_restricted_" + std::to_string(rank));
        }
    }

private:
    inline static constexpr int TimeResolution = 180;

    inline static void writeResults(std::vector<Data>& chunkedData, const std::string& fileNameBase) noexcept {
        const std::string dataName = fileNameBase + ".csv";
        const std::string tikzName = fileNameBase + ".txt";
        const std::string plotName = fileNameBase + ".png";
        Geometry::Point canvasSize = Geometry::Point(Construct::XY, 2880, 2880);
        Geometry::Rectangle bb = Geometry::Rectangle::BoundingBox(Geometry::Point(Construct::XY, 0, 0), Geometry::Point(Construct::XY, 24 * 60, 24 * 60 * 60));
        FunctionVisualization<PNG> plot = FunctionVisualization<PNG>::FromBoundingBox(plotName, bb, canvasSize);
        IO::OFStream tikz(tikzName);
        IO::OFStream data(dataName);
        data << "minute,travelTimeMean,optimalTravelTimeMean,diffMin,diffMax,diffMean,diffMedian,difflowerQuantile,diffupperQuantile,incorrectCount0,incorrectCount10,incorrectCount30,impossibleCount\n";

        std::vector<std::string> pgfTravelTime;
        std::vector<std::string> pgfOptimalTravelTime;
        std::vector<std::string> pgfDiffMedian;
        std::vector<std::string> pgfLower;
        std::vector<std::string> pgfUpper;
        std::vector<std::string> pgfIncorrect0;
        std::vector<std::string> pgfIncorrect10;

        for (size_t chunk = 0; chunk < chunkedData.size(); chunk++) {
            std::sort(chunkedData[chunk].travelTimeDiff.begin(), chunkedData[chunk].travelTimeDiff.end());
            const double restrictedMean = Vector::mean(chunkedData[chunk].restrictedTravelTime);
            const double unrestrictedMean = Vector::mean(chunkedData[chunk].unrestrictedTravelTime);
            const double diffMin = Vector::min(chunkedData[chunk].travelTimeDiff);
            const double diffMax = Vector::max(chunkedData[chunk].travelTimeDiff);
            const double diffMean = Vector::mean(chunkedData[chunk].travelTimeDiff);
            const double diffMedian = Vector::median(chunkedData[chunk].travelTimeDiff);
            const double diffLowerQuantile = Vector::percentile(chunkedData[chunk].travelTimeDiff, 0.25);
            const double diffUpperQuantile = Vector::percentile(chunkedData[chunk].travelTimeDiff, 0.75);
            const double incorrectCount0 = Vector::count(chunkedData[chunk].travelTimeDiff, [](const double d){return d > 1.0;}) / static_cast<double>(chunkedData[chunk].size());
            const double incorrectCount10 = Vector::count(chunkedData[chunk].travelTimeDiff, [](const double d){return d > 600.0;}) / static_cast<double>(chunkedData[chunk].size());
            const double incorrectCount30 = Vector::count(chunkedData[chunk].travelTimeDiff, [](const double d){return d > 3600.0;}) / static_cast<double>(chunkedData[chunk].size());
            const double impossibleCount = chunkedData[chunk].impossible / static_cast<double>(chunkedData[chunk].size());

            data << chunk << "," << restrictedMean << "," << unrestrictedMean << "," << diffMin << "," << diffMax << "," << diffMean << "," << diffMedian << "," << diffLowerQuantile << "," << diffUpperQuantile << "," << incorrectCount0 << "," << incorrectCount10 << "," << impossibleCount << "\n";

            const std::string x = (chunk + 1 >= chunkedData.size()) ? ("1440") : (String::lexicalCast<std::string>(chunk * 3));
            pgfTravelTime.emplace_back("(" + x + "," + String::lexicalCast<std::string>(restrictedMean / 3600) + ") ");
            pgfOptimalTravelTime.emplace_back("(" + x + "," + String::lexicalCast<std::string>(unrestrictedMean / 3600) + ") ");
            pgfDiffMedian.emplace_back("(" + x + "," + String::lexicalCast<std::string>(diffMedian / 3600) + ") ");
            pgfLower.emplace_back("(" + x + "," + String::lexicalCast<std::string>(diffLowerQuantile / 3600) + ") ");
            pgfUpper.emplace_back("(" + x + "," + String::lexicalCast<std::string>(diffUpperQuantile / 3600) + ") ");
            pgfIncorrect0.emplace_back("(" + x + "," + String::lexicalCast<std::string>(incorrectCount0) + ") ");
            pgfIncorrect10.emplace_back("(" + x + "," + String::lexicalCast<std::string>(incorrectCount30) + ") ");

            plot.fillRectangle(Geometry::Rectangle::BoundingBox(Geometry::Point(Construct::XY, chunk * 3, diffLowerQuantile - 120), Geometry::Point(Construct::XY, chunk * 3 + 3, diffUpperQuantile + 120)), Color::KITblue >> 0.30);
            plot.fillRectangle(Geometry::Rectangle::BoundingBox(Geometry::Point(Construct::XY, chunk * 3, diffMean - 120), Geometry::Point(Construct::XY, chunk * 3 + 3, diffMean + 120)), Color::KITblue);
            plot.fillRectangle(Geometry::Rectangle::BoundingBox(Geometry::Point(Construct::XY, chunk * 3, diffMedian - 120), Geometry::Point(Construct::XY, chunk * 3 + 3, diffMedian + 120)), Color::KITseablue);

            plot.fillRectangle(Geometry::Rectangle::BoundingBox(Geometry::Point(Construct::XY, chunk * 3, incorrectCount0*24*60*60 - 120), Geometry::Point(Construct::XY, chunk * 3 + 3, incorrectCount0*24*60*60 + 120)), Color::KITred);
            plot.fillRectangle(Geometry::Rectangle::BoundingBox(Geometry::Point(Construct::XY, chunk * 3, incorrectCount10*24*60*60 - 120), Geometry::Point(Construct::XY, chunk * 3 + 3, incorrectCount10*24*60*60 + 120)), Color::KITred >> 0.60);
            plot.fillRectangle(Geometry::Rectangle::BoundingBox(Geometry::Point(Construct::XY, chunk * 3, incorrectCount30*24*60*60 - 120), Geometry::Point(Construct::XY, chunk * 3 + 3, incorrectCount30*24*60*60 + 120)), Color::KITred >> 0.30);
            plot.fillRectangle(Geometry::Rectangle::BoundingBox(Geometry::Point(Construct::XY, chunk * 3, impossibleCount*24*60*60 - 120), Geometry::Point(Construct::XY, chunk * 3 + 3, impossibleCount*24*60*60 + 120)), Color::Red);

            plot.fillRectangle(Geometry::Rectangle::BoundingBox(Geometry::Point(Construct::XY, chunk * 3, restrictedMean - 120), Geometry::Point(Construct::XY, chunk * 3 + 3, restrictedMean + 120)), Color::KITorange);
            plot.fillRectangle(Geometry::Rectangle::BoundingBox(Geometry::Point(Construct::XY, chunk * 3, unrestrictedMean - 120), Geometry::Point(Construct::XY, chunk * 3 + 3, unrestrictedMean + 120)), Color::KITgreen);
        }

        tikz << "\\newcommand{\\plotWidth}{0.951\\textwidth}\n"
             << "\n"
             << "\\begin{tikzpicture}\n"
             << "\\pgfplotsset{\n"
             << "   grid style={KITblack20,line width = 0.2pt,dash pattern = on 2pt off 1pt},\n"
             << "}\n"
             << "\n"
             << "\\begin{axis}[\n"
             << "   height=6.7cm,\n"
             << "   width=\\plotWidth,\n"
             << "   xlabel={Maximum walking time [min]},\n"
             << "   xmin=0,\n"
             << "   xmax=1440,\n"
             << "   xtick={0, 360, 720, 1080, 1440},\n"
             << "   xticklabel=\\pgfmathparse{(\\tick / 60)}${\\pgfmathprintnumber{\\pgfmathresult}}$:$00$,\n"
             << "   xtick pos=left,\n"
             << "   xtick align=outside,\n"
             << "   ymin=0,\n"
             << "   ymax=8,\n"
             << "   ytick={0, 2, 4, 6, 8},\n"
             << "   yticklabel=\\pgfmathparse{\\tick}${\\pgfmathprintnumber{\\pgfmathresult}}$h,\n"
             << "   ytick align=outside\n"
             << "]\n"
             << "\n"
             << "\\addplot [fill=KITseablue30,draw=none] coordinates {\n"
             << "   ";
        for (size_t i = 0; i < pgfUpper.size(); i++) {
            tikz << pgfUpper[i];
        }
        for (size_t i = 0; i < pgfLower.size(); i++) {
            tikz << pgfLower[pgfLower.size() - i - 1];
        }
        tikz << "};\n"
             << "\n"
             << "\\end{axis}\n"
             << "\n"
             << "\\begin{axis}[\n"
             << "   height=6.7cm,\n"
             << "   width=\\plotWidth,\n"
             << "   ticklabel pos=right,\n"
             << "   xmin=0,\n"
             << "   xmax=1440,\n"
             << "   xmajorticks=false,\n"
             << "   xtick={0, 360, 720, 1080, 1440},\n"
             << "   xticklabel=\\pgfmathparse{(\\tick)}${}$,\n"
             << "   minor x tick num={5},\n"
             << "   ymin=0,\n"
             << "   ymax=1,\n"
             << "   ytick={0, 0.25, 0.5, 0.75, 1},\n"
             << "   yticklabel=\\pgfmathparse{\\tick*100}${\\pgfmathprintnumber{\\pgfmathresult}}\\%$,\n"
             << "   yticklabel shift=2pt,\n"
             << "   minor y tick num={1},\n"
             << "   grid=both\n"
             << "]\n"
             << "\n"
             << "\\addplot [color=KITred50,line width=1.5pt,dash pattern = on 1pt off 2pt,line cap=round] coordinates {\n"
             << "   ";
        for (size_t i = 0; i < pgfIncorrect10.size(); i++) {
            tikz << pgfIncorrect10[i];
        }
        tikz << "};\n"
             << "\n"
             << "\\addplot [color=KITred,line width=1.5pt,dash pattern = on 1pt off 2pt,line cap=round] coordinates {\n"
             << "   ";
        for (size_t i = 0; i < pgfIncorrect0.size(); i++) {
            tikz << pgfIncorrect0[i];
        }
        tikz << "};\n"
             << "\n"
             << "\\end{axis}\n"
             << "\n"
             << "\\begin{axis}[\n"
             << "   height=6.7cm,\n"
             << "   width=\\plotWidth,\n"
             << "   xmin=0,\n"
             << "   xmax=1440,\n"
             << "   xticklabel=\\pgfmathparse{(\\tick)}${}$,\n"
             << "   xtick style={draw=none},\n"
             << "   ymin=0,\n"
             << "   ymax=8,\n"
             << "   ytick={0, 2, 4, 6, 8},\n"
             << "   yticklabel=\\pgfmathparse{(\\tick)}${}$,\n"
             << "   ytick style={draw=none}\n"
             << "]\n"
             << "\n"
             << "\\addplot [color=KITseablue,line width=1.5pt] coordinates {\n"
             << "   ";
        for (size_t i = 0; i < pgfDiffMedian.size(); i++) {
            tikz << pgfDiffMedian[i];
        }
        tikz << "};\n"
             << "\n"
             << "\\addplot [color=KITorange,line width=1.5pt] coordinates {\n"
             << "   ";
        for (size_t i = 0; i < pgfTravelTime.size(); i++) {
            tikz << pgfTravelTime[i];
        }
        tikz << "};\n"
             << "\n"
             << "\\addplot [color=KITgreen,line width=1.5pt] coordinates {\n"
             << "   ";
        for (size_t i = 0; i < pgfOptimalTravelTime.size(); i++) {
            tikz << pgfOptimalTravelTime[i];
        }
        tikz << "};\n"
             << "\n"
             << "\\end{axis}\n"
             << "\n"
             << "\\end{tikzpicture}\n";

        plot.drawXAxis(60, [](const int x){return std::to_string(x / 60) + "h";});
        plot.drawYAxis(60 * 60, [](const int y){return std::to_string(y / 3600) + "h";});
    }
};

class FindUnrestrictedWalkingExamples : public ParameterizedCommand {

public:
    FindUnrestrictedWalkingExamples(BasicShell& shell) :
        ParameterizedCommand(shell, "findUnrestrictedWalkingExamples", "Finds examples where unrestricted walking improves the arrival time.") {
        addParameter("ULTRA-RAPTOR data");
        addParameter("Bucket CH data");
        addParameter("Transitive RAPTOR data");
        addParameter("Number of examples");
        addParameter("Arrival time threshold");
        addParameter("Walking threshold");
        addParameter("Output file");
        addParameter("Restrict initial walking?", "true", {"true", "false"});
    }

    virtual void execute() noexcept {
        RAPTOR::Data ultraRaptorData(getParameter("ULTRA-RAPTOR data"));
        ultraRaptorData.useImplicitDepartureBufferTimes();
        ultraRaptorData.printInfo();
        CH::CH ch(getParameter("Bucket CH data"));
        RAPTOR::Data transitiveRaptorData(getParameter("Transitive RAPTOR data"));
        transitiveRaptorData.useImplicitDepartureBufferTimes();
        transitiveRaptorData.printInfo();
        ULTRARAPTOR ultraRaptor(ultraRaptorData, ch);

        if (getParameter<bool>("Restrict initial walking?")) {
            TransitiveRAPTOR transitiveRaptor(transitiveRaptorData);
            run(ultraRaptor, transitiveRaptor, transitiveRaptorData);
        } else {
            ULTRARAPTOR transitiveRaptor(transitiveRaptorData, ch);
            run(ultraRaptor, transitiveRaptor, transitiveRaptorData);
        }
    }

private:
    using ULTRARAPTOR = RAPTOR::ULTRARAPTOR<RAPTOR::NoProfiler, false>;
    using TransitiveRAPTOR = RAPTOR::RAPTOR<true, RAPTOR::NoProfiler, true, false, false>;

    template<typename TRANSITIVE_RAPTOR>
    inline void run(ULTRARAPTOR& ultraRaptor, TRANSITIVE_RAPTOR& transitiveRaptor, const RAPTOR::Data& raptorData) const noexcept {
        const size_t numberOfExamples = getParameter<bool>("Number of examples");
        const int arrivalTimeThreshold = getParameter<int>("Arrival time threshold");
        const int walkingThreshold = getParameter<int>("Walking threshold");
        const std::string outputFile = getParameter("Output file");

        IO::OFStream out(outputFile + ".csv");
        IO::OFStream unlimitedJourneys(outputFile + "_journeys_unlimited.txt");
        IO::OFStream limitedJourneys(outputFile + "_journeys_limited.txt");
        out << "From;From ID;To;To ID;Departure Time;Arrival Time Unlimited;Arrival Time Limited;Time Saved;Walking Time Unlimited;Walking Time Limited\n";

        srand(42);
        Progress progress(numberOfExamples);
        size_t foundExamples = 0;

        while (foundExamples < numberOfExamples) {
            const StopId source(rand() % raptorData.numberOfStops());
            const StopId target(rand() % raptorData.numberOfStops());
            const int departureTime(rand() % (24 * 60 * 60));
            ultraRaptor.run(source, departureTime, target);
            transitiveRaptor.run(source, departureTime, target);
            const int ultraArrival = ultraRaptor.getEarliestArrivalTime(target);
            const int transitiveArrival = transitiveRaptor.getEarliestArrivalTime(target);
            const RAPTOR::Journey& ultraJourney = ultraRaptor.getEarliestJourney(target);
            const RAPTOR::Journey& transitiveJourney = transitiveRaptor.getEarliestJourney(target);
            const int ultraWalkingTime = totalTransferTime(ultraJourney);
            const int transitiveWalkingTime = totalTransferTime(transitiveJourney);
            if (transitiveArrival - ultraArrival >= arrivalTimeThreshold && ultraWalkingTime <= walkingThreshold) {
                out << raptorData.stopData[source].name;
                out << ";" << int(source);
                out << ";" << raptorData.stopData[target].name;
                out << ";" << int(target);
                out << ";" << String::secToTime(departureTime, never, true);
                out << ";" << String::secToTime(ultraArrival, never, true);
                out << ";" << String::secToTime(transitiveArrival, never, true);
                out << ";" << String::secToTime(transitiveArrival - ultraArrival, never, true);
                out << ";" << String::secToTime(ultraWalkingTime, never, true);
                out << ";" << String::secToTime(transitiveWalkingTime, never, true);
                out << "\n";
                out.flush();
                unlimitedJourneys << raptorData.journeyToText(ultraJourney) << "\n";
                unlimitedJourneys.flush();
                limitedJourneys << raptorData.journeyToText(transitiveJourney) << "\n";
                limitedJourneys.flush();
                foundExamples++;
                progress++;
            }
        }
    }
};
