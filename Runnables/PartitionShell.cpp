#include <sched.h>

#include "../Helpers/Console/CommandLineParser.h"
#include "../Shell/Shell.h"
#include "Commands/Partition.h"

using namespace Shell;

int main(int argc, char** argv) {
    CommandLineParser clp(argc, argv);
    checkAsserts();
    Shell::Shell shell;
    new DrawGreedyCells(shell);
    new ComputeNetworkNestedDissection(shell);
    new RunGraphInertialFlow(shell);
    new RunNetworkInertialFlow(shell);
    new DrawVertexPartition(shell);
    new DrawNestedDissection(shell);
    new ComputeSampleGraph(shell);
    shell.run();
    return 0;
}
