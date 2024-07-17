#include <sched.h>

#include "../Helpers/Console/CommandLineParser.h"
#include "../Shell/Shell.h"
#include "Commands/MultimodalCoverage.h"

using namespace Shell;

int main(int argc, char** argv) {
    CommandLineParser clp(argc, argv);
    checkAsserts();
    Shell::Shell shell;
    new ComputeMultimodalParetoSets(shell);
    new ComputeMultimodalRestrictedParetoSets(shell);
    new ComputeMultimodalULTRAMcRAPTORCoverage(shell);
    new ComputeMultimodalUBMRAPTORCoverage(shell);
    new ComputeMultimodalBoundedMcRAPTORCoverage(shell);
    new ComputeMultimodalUBMHydRACoverage(shell);
    shell.run();
    return 0;
}
