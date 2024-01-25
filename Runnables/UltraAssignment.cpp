#include <iostream>
#include <algorithm>
#include <random>
#include <vector>
#include <string>

#include "Commands/Assignment.h"
#include "Commands/NetworkAnalysis.h"
#include "Commands/ULTRAPreprocessing.h"

#include "../Shell/Shell.h"

#include "../Helpers/Debug.h"
#include "../Helpers/MultiThreading.h"
#include "../Helpers/Console/CommandLineParser.h"

using namespace Shell;

int main(int argc, char** argv) {
    CommandLineParser clp(argc, argv);
    pinThreadToCoreId(clp.value<int>("core", 1));
    checkAsserts();

    ::Shell::Shell shell;
    new AssignmentGraph(shell);
    new AnalyzeNetwork(shell);
    new AssignmentTest(shell);
    new AssignmentStatistics(shell);
    new GroupAssignment(shell);
    new CSVStatistic(shell);
    new DemandHistogram(shell);
    new BuildMultiModalAssignmentData(shell);
    new BuildMultiModalDemandData(shell);
    new RandomDemand(shell);
    new UltraGroupAssignment(shell);
    new ComputeStopToStopShortcuts(shell);
    new TransferShortcutsForAssignment(shell);
    new BuildStations(shell);
    shell.run();
    return 0;
}
