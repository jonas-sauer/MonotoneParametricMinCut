#include <sched.h>

#include "../Helpers/Console/CommandLineParser.h"
#include "../Shell/Shell.h"
#include "Commands/NetworkOrder.h"

using namespace Shell;

int main(int argc, char** argv) {
    CommandLineParser clp(argc, argv);
    checkAsserts();
    Shell::Shell shell;
    new MeasureDijkstra(shell);
    new MeasureBFS(shell);
    new MeasurePHAST(shell);
    new VisualizeVertexOrder(shell);
    new ReorderGraph(shell);
    new ReorderNetwork(shell);
    new ReorderNetworkStops(shell);
    new CreateRandomVertexOrder(shell);
    new CreateDFSVertexOrder(shell);
    new CreateDFSStopOrder(shell);
    new CreateHilbertCurveVertexOrder(shell);
    shell.run();
    return 0;
}
