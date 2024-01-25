#include <iostream>
#include <algorithm>
#include <random>
#include <vector>
#include <string>

#include "Commands/Assignment.h"
#include "Commands/NetworkIO.h"
#include "Commands/Visualization.h"

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
    new ParseCSAFromCSV(shell);
    new SanitizeDemand(shell);
    new GroupAssignment(shell);
    new SequentialAssignment(shell);
    new CapacityAssignment(shell);
    new LockstepCapacityAssignment(shell);
    new FixedCapacityAssignment(shell);
    new ConnectionBasedLockstepCapacityAssignment(shell);
    new ParseCapacities(shell);
    new BuildCapacityAssignmentTestNetwork(shell);
    new AssignmentStatistics(shell);
    new DrawAssignment(shell);
    shell.run();
    return 0;
}
