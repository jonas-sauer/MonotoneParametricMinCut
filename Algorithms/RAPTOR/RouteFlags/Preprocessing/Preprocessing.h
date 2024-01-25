#pragma once

#include <algorithm>

#include "../../../../DataStructures/RAPTOR/RouteFlags/RouteFlagsData.h"
#include "../../../../Helpers/MultiThreading.h"
#include "../../../../Helpers/Timer.h"
#include "RangeSearch.h"

namespace RAPTOR::RouteFlags {

template<typename FLAGS, typename FINE_PARTITION_TYPE, bool USE_ALTERNATIVE_BACKWARD_SEARCH = false>
class Preprocessing {

public:
    using Flags = FLAGS;
    using FinePartitionType = FINE_PARTITION_TYPE;
    static constexpr bool UseAlternativeBackwardSearch = USE_ALTERNATIVE_BACKWARD_SEARCH;
    using Type = Preprocessing<Flags, FinePartitionType, UseAlternativeBackwardSearch>;
    using RouteFlagsDataType = RouteFlagsData<Flags, FinePartitionType>;
    using ForwardRangeRAPTOR = RangeSearch<Flags, FinePartitionType, FORWARD, UseAlternativeBackwardSearch, false>;
    using BackwardRangeRAPTOR = RangeSearch<Flags, FinePartitionType, BACKWARD, UseAlternativeBackwardSearch, false>;

    Preprocessing(RouteFlagsDataType& flagsData) :
        flagsData(flagsData) {
    }

    inline void run(ThreadPinning& threadPinner, const int minTime = -never, const int maxTime = never, const bool verbose = true) noexcept {
        if (verbose) std::cout << "Running Route-Flags preprocessing (" << threadPinner.numberOfThreads << " threads, " << flagsData.getCoarsePartition().numberOfBoundaryVertices() << " boundary vertices)... " << std::flush;

        Timer totalTimer;
        totalTimer.restart();

        flagsData.markIntraCellStopEvents();

        const size_t coarseCells = flagsData.getCoarsePartition().numberOfCells();
        const size_t fineCells = flagsData.getFinePartition().numberOfCells();
        ShortcutFlags forwardShortcutFlags(coarseCells, fineCells, flagsData.getShortcutGraph(FORWARD));
        ShortcutFlags backwardShortcutFlags(coarseCells, fineCells, flagsData.getShortcutGraph(BACKWARD));

        omp_set_num_threads(threadPinner.numberOfThreads);
        #pragma omp parallel
        {
            threadPinner.pinThread();
            #pragma omp flush

            ForwardRangeRAPTOR forwardQuery(flagsData);
            BackwardRangeRAPTOR backwardQuery(flagsData);
            //Timer vertexTimer;

            const std::vector<Vertex>& boundaryVertices = flagsData.getCoarsePartition().boundaryVertices();
            #pragma omp for schedule(dynamic)
            for (size_t s = 0; s < boundaryVertices.size(); s++) {
                if (verbose) std::cout << "Vertex " << boundaryVertices[s] << " (" << s << "/" << boundaryVertices.size() << ")" << std::endl;
                //vertexTimer.restart();
                forwardQuery.run(boundaryVertices[s], minTime, maxTime);
                backwardQuery.run(boundaryVertices[s], -maxTime, -minTime);
                //if (verbose) std:: cout << "Took " << String::msToString(vertexTimer.elapsedMilliseconds()) << std::endl;
            }

            #pragma omp critical
            {
                forwardShortcutFlags.incorporate(backwardQuery.getShortcutFlags());
                backwardShortcutFlags.incorporate(forwardQuery.getShortcutFlags());
            }
        }

        if (verbose) std::cout << "Finished" << std::endl;

        flagsData.computeExtraFlags();
        flagsData.buildShortcuts(forwardShortcutFlags, backwardShortcutFlags);

        if (verbose) {
            std::cout << "Time: " << String::msToString(totalTimer.elapsedMilliseconds()) << std::endl;
            flagsData.printAverageStatistics();
        }
    }

private:
    RouteFlagsDataType& flagsData;
};

}
