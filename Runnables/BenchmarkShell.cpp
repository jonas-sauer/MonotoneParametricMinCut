#include <sched.h>

#include "../Helpers/Console/CommandLineParser.h"
#include "../Shell/Shell.h"
#include "Commands/BenchmarkMcULTRA.h"
#include "Commands/BenchmarkMultimodal.h"
#include "Commands/BenchmarkULTRA.h"
#include "Commands/BenchmarkULTRAPHAST.h"

using namespace Shell;

int main(int argc, char** argv) {
    CommandLineParser clp(argc, argv);
    checkAsserts();
    Shell::Shell shell;

    //ULTRA
    new GenerateRandomVertexQueries(shell);
    new GenerateRandomStopQueries(shell);
    new RunTransitiveCSAQueries(shell);
    new RunDijkstraCSAQueries(shell);
    new RunHLCSAQueries(shell);
    new RunULTRACSAQueries(shell);
    new RunParetoCSAQueries(shell);
    new RunParetoULTRACSAQueries(shell);
    new RunTransitiveProfileCSAQueries(shell);
    new RunTransitiveRAPTORQueries(shell);
    new RunDijkstraRAPTORQueries(shell);
    new RunHLRAPTORQueries(shell);
    new RunULTRARAPTORQueries(shell);
    new RunTransitiveTBQueries(shell);
    new RunULTRATBQueries(shell);
    new RunHLULTRATBQueries(shell);

    //ULTRA-PHAST
    new RunOneToAllDijkstraCSAQueriesToVertices(shell);
    new RunOneToManyDijkstraCSAQueriesToStops(shell);
    new RunUPCSAQueries(shell);
    new RunOneToAllDijkstraRAPTORQueriesToVertices(shell);
    new RunOneToManyDijkstraRAPTORQueriesToStops(shell);
    new RunUPRAPTORQueries(shell);
    new RunUPTBQueries(shell);
    new CreateBallTargetSets(shell);
    new BuildCoreCHForTargetSets(shell);
    new BuildUPCHForTargetSets(shell);
    new RunOneToManyDijkstraCSAQueriesToBall(shell);
    new RunUPCSAQueriesToBall(shell);
    new RunOneToManyDijkstraRAPTORQueriesToBall(shell);
    new RunUPRAPTORQueriesToBall(shell);

    //McULTRA
    new RunTransitiveMcRAPTORQueries(shell);
    new RunMCRQueries(shell);
    new RunULTRAMcRAPTORQueries(shell);
    new RunULTRAMcTBQueries(shell);
    new RunTransitiveBoundedMcRAPTORQueries(shell);
    new RunUBMRAPTORQueries(shell);
    new RunNewUBMRAPTORQueries(shell);
    new RunUBMTBRAPTORQueries(shell);
    new RunUBMTBQueries(shell);
    new RunUBMHydRAQueries(shell);
    new ComputeTransferTimeSavings(shell);

    //Multiple transfer modes
    new RunMultimodalMCRQueries(shell);
    new RunMultimodalULTRAMcRAPTORQueries(shell);
    new RunMultimodalUBMRAPTORQueries(shell);
    new RunMultimodalUBMHydRAQueries(shell);

    shell.run();
    return 0;
}
