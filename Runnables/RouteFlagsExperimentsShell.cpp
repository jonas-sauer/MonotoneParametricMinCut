#include <sched.h>

#include "../Helpers/Debug.h"
#include "../Helpers/Console/CommandLineParser.h"
#include "../Shell/Shell.h"
#include "Commands/RouteFlagsExperiments.h"

using namespace Shell;

int main(int argc, char** argv) {
    CommandLineParser clp(argc, argv);
    checkAsserts();
    Shell::Shell shell;
    new BuildPartitionCoreCH(shell);
    new BuildRouteFlagsData(shell);
    new RunRouteFlagsPreprocessing(shell);
    new VisualizeIntraCellConnections(shell);
    new VisualizeTransferShortcuts(shell);
    new VisualizeRouteFlagsPreprocessing(shell);
    new RunRouteFlagsQuery(shell);
    //new VisualizeRouteFlagsQuery(shell);
    new RunRouteFlagsQueries(shell);
    new CheckFlags(shell);
    new CheckShortcut(shell);
    new CellsOfRoute(shell);
    new CellsOfFootpath(shell);
    shell.run();
    return 0;
}
