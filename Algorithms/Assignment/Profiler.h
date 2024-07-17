#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/Assignment/ConnectionLoadData.h"

#include "../../Helpers/Types.h"
#include "../../Helpers/Timer.h"
#include "../../Helpers/String/String.h"
#include "../../Helpers/HighlightText.h"

namespace Assignment {

class NoPATProfiler {

public:
    inline void start() const noexcept {}
    inline void done() const noexcept {}

    inline void startInitialization() const noexcept {}
    inline void doneInitialization() const noexcept {}

    inline void scanConnection(const ConnectionId) const noexcept {}
    inline void relaxEdge(const Edge) const noexcept {};

    inline void addToProfile() const noexcept {}
    inline void insertToProfile() const noexcept {}
    inline void evaluateProfile() const noexcept {}

    inline void copyEntry() const noexcept {}

    inline std::string getStatistics() const noexcept {return "";}
};

class PATProfiler : public NoPATProfiler {

public:
    inline void start() noexcept {
        queries = 0;
        scannedConnections = 0;
        relaxedEdges = 0;
        profileAdds = 0;
        profileInsertions = 0;
        profileEvaluations = 0;
        copiedEntries = 0;
        timer.restart();
    }
    inline void done() noexcept {
        std::cout << getStatistics() << std::endl;
    }

    inline void startInitialization() noexcept {
        queries++;
    }

    inline void scanConnection(const ConnectionId) noexcept {
        scannedConnections++;
    }
    inline void relaxEdge(const Edge) noexcept {
        relaxedEdges++;
    };

    inline void addToProfile() noexcept {
        profileAdds++;
    }
    inline void insertToProfile() noexcept {
        profileInsertions++;
    }
    inline void evaluateProfile() noexcept {
        profileEvaluations++;
    }

    inline void copyEntry() noexcept {
        copiedEntries++;
    }

    inline std::string getStatistics() noexcept {
        double time = timer.elapsedMilliseconds();
        std::stringstream result;
        result << "\nComputePATs Statistics:\n"
               << "   Number of queries: " << String::prettyInt(queries) << "\n"
               << "   Total time: " << String::msToString(time) << " (" << String::msToString(time / queries) << " per query)\n"
               << "\n"
               << "   Variable            " << std::setw(15) << "Value"                               << std::setw(15) << "per query"                                                             << std::setw(12) << "per ms" << "\n"
               << "   scannedConnections: " << std::setw(15) << String::prettyInt(scannedConnections) << std::setw(15) << String::prettyDouble(scannedConnections / static_cast<double>(queries)) << std::setw(12) << String::prettyDouble(scannedConnections / time) << "\n"
               << "   relaxedEdges:       " << std::setw(15) << String::prettyInt(relaxedEdges)       << std::setw(15) << String::prettyDouble(relaxedEdges / static_cast<double>(queries))       << std::setw(12) << String::prettyDouble(relaxedEdges / time) << "\n"
               << "   profileAdds:        " << std::setw(15) << String::prettyInt(profileAdds)        << std::setw(15) << String::prettyDouble(profileAdds / static_cast<double>(queries))        << std::setw(12) << String::prettyDouble(profileAdds / time) << "\n"
               << "   profileInsertions:  " << std::setw(15) << String::prettyInt(profileInsertions)  << std::setw(15) << String::prettyDouble(profileInsertions / static_cast<double>(queries))  << std::setw(12) << String::prettyDouble(profileInsertions / time) << "\n"
               << "   profileEvaluations: " << std::setw(15) << String::prettyInt(profileEvaluations) << std::setw(15) << String::prettyDouble(profileEvaluations / static_cast<double>(queries)) << std::setw(12) << String::prettyDouble(profileEvaluations / time) << "\n"
               << "   copiedEntries:      " << std::setw(15) << String::prettyInt(copiedEntries)      << std::setw(15) << String::prettyDouble(copiedEntries / static_cast<double>(queries))      << std::setw(12) << String::prettyDouble(copiedEntries / time) << "\n";
        return result.str();
    }

private:
    Timer timer;

    long long queries{0};
    long long scannedConnections{0};
    long long relaxedEdges{0};
    long long profileAdds{0};
    long long profileInsertions{0};
    long long profileEvaluations{0};
    long long copiedEntries{0};

};

class NoProfiler {

public:
    inline void initialize() const noexcept {}
    inline void initialize(const CSA::Data&) const noexcept {}

    inline void start() const noexcept {}
    inline void done() const noexcept {}

    inline void startPATComputation() const noexcept {}
    inline void donePATComputation() const noexcept {}

    inline void startAssignment() const noexcept {}
    inline void doneAssignment() const noexcept {}

    inline void startInitialWalking() const noexcept {}
    inline void doneInitialWalking() const noexcept {}

    inline void startCycleElimination() const noexcept {}
    inline void doneCycleElimination() const noexcept {}

    inline void setPathsPerPassenger(const double) const noexcept {}

    inline double getSetupTime() const noexcept {return 0;}
    inline double getPATComputationTime() const noexcept {return 0;}
    inline double getAssignmentTime() const noexcept {return 0;}
    inline double getInitialWalkingTime() const noexcept {return 0;}
    inline double getCycleEliminationTime() const noexcept {return 0;}
    inline double getTotalTime() const noexcept {return 0;}

    inline void printStatistics() const noexcept {}

    inline void startAssignmentForDestination(const int) const noexcept {}
    inline void assignConnection(const ConnectionId) const noexcept {}
    inline void doneAssignmentForDestination(const int) const noexcept {}

    inline void moveGroups(const std::string&, const std::string&) const noexcept {}
    inline void moveGroupsDestination(const int) const noexcept {}
    inline void moveGroupsPATs(const double, const double) const noexcept {}
    inline void moveGroupsProbabilities(const std::array<int, 3>) const noexcept {}
    inline void moveGroupsSizes(const std::array<int, 2>) const noexcept {}

    inline void distributePassengersPATs(const std::vector<int>&, const std::vector<int>&) const noexcept {}
    inline void distributePassengersPATs(const std::vector<double>&, const std::vector<ConnectionId>&) const noexcept {}
    inline void distributePassengersProbabilities(const std::vector<int>&) const noexcept {}
    inline void distributePassengersSizes(const std::vector<size_t>&) const noexcept {}

    inline void printConnectionStatistics(const ConnectionId, const double, const double, const double, const BoardingData&, const BoardingData&, const BoardingData&) const noexcept {}

    inline NoProfiler& operator+=(const NoProfiler&) noexcept {
        return *this;
    }
};

class TimeProfiler : public NoProfiler {

public:
    inline void initialize() noexcept {
        numberOfCalculations = 0;
        numberOfPATComputations = 0;
        timeForEverything = 0;
        timeForPATComputation = 0;
        timeForAssignment = 0;
        timeForInitialWalking = 0;
        timeForCycleElimination = 0;
        pathsPerPassenger = 0;
    }
    inline void initialize(const CSA::Data&) noexcept {
        initialize();
    }

    inline void start() noexcept {
        numberOfCalculations++;
        timerForEverything.restart();
    }
    inline void done() noexcept {
        timeForEverything += timerForEverything.elapsedMicroseconds();
    }

    inline void startPATComputation() noexcept {
        numberOfPATComputations++;
        timerForPATComputation.restart();
    }
    inline void donePATComputation() noexcept {
        timeForPATComputation += timerForPATComputation.elapsedMicroseconds();
    }

    inline void startAssignment() noexcept {
        timerForAssignment.restart();
    }
    inline void doneAssignment() noexcept {
        timeForAssignment += timerForAssignment.elapsedMicroseconds();
    }

    inline void startInitialWalking() noexcept {
        timerForInitialWalking.restart();
    }
    inline void doneInitialWalking() noexcept {
        timeForInitialWalking += timerForInitialWalking.elapsedMicroseconds();
    }

    inline void startCycleElimination() noexcept {
        timerForCycleElimination.restart();
    }
    inline void doneCycleElimination() noexcept {
        timeForCycleElimination += timerForCycleElimination.elapsedMicroseconds();
    }

    inline void setPathsPerPassenger(const double n) noexcept {
        pathsPerPassenger = n;
    }
    inline double getSetupTime() const noexcept {
        return (timeForEverything - (timeForPATComputation + timeForAssignment + timeForCycleElimination)) / numberOfCalculations;
    }
    inline double getPATComputationTime() const noexcept {
        return timeForPATComputation / numberOfCalculations;
    }
    inline double getAssignmentTime() const noexcept {
        return timeForAssignment / numberOfCalculations;
    }
    inline double getInitialWalkingTime() const noexcept {
        return timeForInitialWalking / numberOfCalculations;
    }
    inline double getCycleEliminationTime() const noexcept {
        return timeForCycleElimination / numberOfCalculations;
    }
    inline double getTotalTime() const noexcept {
        return timeForEverything / numberOfCalculations;
    }

    inline void printStatistics() const noexcept {
        std::cout << "Setup:           " << String::musToString(getSetupTime()) << std::endl;
        std::cout << "PAT:             " << String::musToString(getPATComputationTime()) << std::endl;
        std::cout << "Initial walking: " << String::musToString(getInitialWalkingTime()) << std::endl;
        std::cout << "Assignment:      " << String::musToString(getAssignmentTime()) << std::endl;
        std::cout << "Cycle removal:   " << String::musToString(getCycleEliminationTime()) << std::endl;
        std::cout << "Total time:      " << String::musToString(getTotalTime()) << std::endl;
        std::cout << "#Targets:        " << String::prettyInt(numberOfPATComputations) << std::endl;
        std::cout << "#Paths:          " << String::prettyDouble(pathsPerPassenger) << std::endl;
    }

    inline TimeProfiler& operator+=(const TimeProfiler& other) noexcept {
        numberOfPATComputations += other.numberOfPATComputations;
        timeForPATComputation += other.timeForPATComputation;
        timeForAssignment += other.timeForAssignment;
        timeForInitialWalking += other.timeForInitialWalking;
        timeForCycleElimination += other.timeForCycleElimination;
        return *this;
    }

private:
    int numberOfCalculations;
    int numberOfPATComputations;

    Timer timerForEverything;
    Timer timerForPATComputation;
    Timer timerForAssignment;
    Timer timerForInitialWalking;
    Timer timerForCycleElimination;

    double timeForEverything;
    double timeForPATComputation;
    double timeForAssignment;
    double timeForInitialWalking;
    double timeForCycleElimination;

    double pathsPerPassenger;
};

class DecisionProfiler : public NoProfiler {

public:
    inline void initialize(const CSA::Data& data) noexcept {
        this->data = &data;
    }

    inline void startAssignmentForDestination(const int destination) const noexcept {
        blue("Current Destination: ", destination) << std::endl;
    }
    inline void assignConnection(const ConnectionId connection) const noexcept {
        const CSA::Connection& c = data->connections[connection];
        blue("\rConnection ", connection, ", time: ", String::secToTime(c.departureTime), "-", String::secToTime(c.arrivalTime), "   (", c.departureTime, ", ", c.arrivalTime, ", ", c.departureStopId, " -> ", c.arrivalStopId, ")                 ");
    }
    inline void doneAssignmentForDestination(const int) const noexcept {
        std::cout << std::endl << std::endl;
    }

    inline void moveGroups(const std::string& a, const std::string& b) const noexcept {
        yellow("\n   Decision:    ", std::setw(16), a, std::setw(16), b) << std::endl;
    }

    inline void moveGroupsDestination(const int destination) const noexcept {
        blue("   Destination: ", destination) << std::endl;
    }

    inline void moveGroupsPATs(const double a, const double b) const noexcept {
        std::cout << "   PATs:        " << std::setw(16) << String::prettyInt(a) << std::setw(16) << String::prettyInt(b) << grey(std::setw(16), String::prettyInt(abs(a - b))) << std::endl;
    }
    inline void moveGroupsProbabilities(const std::array<int, 3> a) const noexcept {
        std::cout << "   Probability: " << std::setw(16) << String::percent(a[0] / static_cast<double>(a[2])) << std::setw(16) << String::percent(a[1] / static_cast<double>(a[2])) << std::endl;
    }
    inline void moveGroupsSizes(const std::array<int, 2> a) const noexcept {
        std::cout << "   Passengers:  " << std::setw(16) << String::prettyInt(a[0]) << std::setw(16) << String::prettyInt(a[1]) << std::endl;
    }

    inline void distributePassengersPATs(const std::vector<int>& v, const std::vector<int>& t) const noexcept {
        if (v.empty()) return;
        yellow("\n   Distribution:") << std::endl;
        std::cout << "   Time:        ";
        for (const int i : t) {
            std::cout << std::setw(16) << String::secToTime(i);
        }
        std::cout << std::endl;
        std::cout << "                ";
        for (const int i : t) {
            std::cout << grey(std::setw(16), String::prettyInt(i));
        }
        std::cout << std::endl;
        std::cout << "   PATs:        ";
        for (const int i : v) {
            std::cout << std::setw(16) << String::prettyInt(i);
        }
        std::cout << std::endl;
    }
    inline void distributePassengersPATs(const std::vector<double>& v, const std::vector<ConnectionId>& c) const noexcept {
        if (v.empty()) return;
        yellow("\n   Distribution:") << std::endl;
        std::cout << "   Time:        ";
        for (const ConnectionId i : c) {
            std::cout << std::setw(16) << ((i == noConnection) ? "target" : String::secToTime(data->connections[i].departureTime));
        }
        std::cout << std::endl;
        std::cout << "                ";
        for (const ConnectionId i : c) {
            std::cout << grey(std::setw(16), (i == noConnection) ? "target" : String::prettyInt(data->connections[i].departureTime));
        }
        std::cout << std::endl;
        std::cout << "   PATs:        ";
        for (const double i : v) {
            std::cout << std::setw(16) << String::prettyDouble(i);
        }
        std::cout << std::endl;
    }
    inline void distributePassengersProbabilities(const std::vector<int>& v) const noexcept {
        if (v.size() <= 1) return;
        std::cout << "   Probability: ";
        for (size_t i = 0; i < v.size() - 1; i++) {
            std::cout << std::setw(16) << String::percent(v[i] / static_cast<double>(v.back()));
        }
        std::cout << std::endl;
    }
    inline void distributePassengersSizes(const std::vector<size_t>& v) const noexcept {
        if (v.empty()) return;
        std::cout << "   Passengers:  ";
        for (const int i : v) {
            std::cout << std::setw(16) << String::prettyInt(i);
        }
        std::cout << std::endl;
    }
    inline void printConnectionStatistics(const ConnectionId i, const double oldLoad, const double currentLoad, const double newLoad, const BoardingData& oldBoardingData, const BoardingData& currentBoardingData, const BoardingData& newBoardingData) const noexcept {
        std::cout << "Connection " << i << ":" << std::endl;
        std::cout << "\tLoad: " << oldLoad << " " << currentLoad << " " << newLoad << std::endl;
        std::cout << "\tBoarding demand: " << oldBoardingData.boardingDemand << " " << currentBoardingData.boardingDemand << " " << newBoardingData.boardingDemand << std::endl;
        std::cout << "\tBoarding capacity: " << oldBoardingData.boardingCapacity << " " << currentBoardingData.boardingCapacity << " " << newBoardingData.boardingCapacity << std::endl;
        std::cout << "\tBoarding probability: " << oldBoardingData.boardingProbability() << " " << currentBoardingData.boardingProbability() << " " << newBoardingData.boardingProbability() << std::endl;
    }
    inline DecisionProfiler& operator+=(const DecisionProfiler&) noexcept {
        return *this;
    }

private:
    const CSA::Data* data{nullptr};

};

struct NoConnectionProfiler {
    inline void initialize(const size_t) noexcept {}
    inline void startCollectIncomingGroups() noexcept {}
    inline void stopCollectIncomingGroups() noexcept {}
    inline void startCalculateDecisionPATs() noexcept {}
    inline void stopCalculateDecisionPATs() noexcept {}
    inline void startBoarding() noexcept {}
    inline void stopBoarding() noexcept {}
    inline void startAssignment() noexcept {}
    inline void stopAssignment() noexcept {}
    inline void startAlighting() noexcept {}
    inline void stopAlighting() noexcept {}
    inline void startTarget() noexcept {}
    inline void stopTarget() noexcept {}
    inline void startWalking() noexcept {}
    inline void stopWalking() noexcept {}
    inline void doneConnection() noexcept {}
    inline void print() noexcept {}
};

struct ConnectionProfiler {
   inline void initialize(const size_t frequency) noexcept {
       collectIncomingGroupsTime = 0.0;
       calculateDecisionPATsTime = 0.0;
       boardingTime = 0.0;
       assignmentTime = 0.0;
       alightingTime = 0.0;
       targetTime = 0.0;
       walkingTime = 0.0;
       processedConnections = 0;
       outputFrequency = frequency;
       totalTimer.restart();
   }

   inline void startCollectIncomingGroups() noexcept {
       timer.restart();
   }

   inline void stopCollectIncomingGroups() noexcept {
       collectIncomingGroupsTime += timer.elapsedMicroseconds();
   }

   inline void startCalculateDecisionPATs() noexcept {
       timer.restart();
   }

   inline void stopCalculateDecisionPATs() noexcept {
       calculateDecisionPATsTime += timer.elapsedMicroseconds();
   }

   inline void startBoarding() noexcept {
       timer.restart();
   }

   inline void stopBoarding() noexcept {
       boardingTime += timer.elapsedMicroseconds();
   }

   inline void startAssignment() noexcept {
       timer.restart();
   }

   inline void stopAssignment() noexcept {
       assignmentTime += timer.elapsedMicroseconds();
   }

   inline void startAlighting() noexcept {
       timer.restart();
   }

   inline void stopAlighting() noexcept {
       alightingTime += timer.elapsedMicroseconds();
   }

   inline void startTarget() noexcept {
       timer.restart();
   }

   inline void stopTarget() noexcept {
       targetTime += timer.elapsedMicroseconds();
   }

   inline void startWalking() noexcept {
       timer.restart();
   }

   inline void stopWalking() noexcept {
       walkingTime += timer.elapsedMicroseconds();
   }

   inline void doneConnection() noexcept {
       processedConnections++;
       if (processedConnections % outputFrequency == 0) {
           print();
       }
   }

   inline void print() noexcept {
       std::cout << "Processed " << processedConnections << " connections" << std::endl;
       std::cout << "\tCollect incoming groups: " << String::musToString(collectIncomingGroupsTime) << std::endl;
       std::cout << "\tCalculate decision PATs: " << String::musToString(calculateDecisionPATsTime) << std::endl;
       std::cout << "\tBoard vs. wait: " << String::musToString(boardingTime) << std::endl;
       std::cout << "\tConnection assignment: " << String::musToString(assignmentTime) << std::endl;
       std::cout << "\tContinue vs. alight: " << String::musToString(alightingTime) << std::endl;
       std::cout << "\tTarget vs. walk: " << String::musToString(targetTime) << std::endl;
       std::cout << "\tWalking after alighting: " << String::musToString(walkingTime) << std::endl;
       std::cout << "\tTotal: " << String::musToString(totalTimer.elapsedMicroseconds()) << std::endl;
   }

   double collectIncomingGroupsTime;
   double calculateDecisionPATsTime;
   double boardingTime;
   double assignmentTime;
   double alightingTime;
   double targetTime;
   double walkingTime;
   size_t processedConnections;
   size_t outputFrequency;
   Timer timer;
   Timer totalTimer;
};

}
