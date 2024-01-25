#include <sched.h>

#include "../Helpers/Console/CommandLineParser.h"
#include "../Shell/Shell.h"
#include "Commands/PublicTransitQueries.h"

using namespace Shell;

int main(int argc, char** argv) {
    CommandLineParser clp(argc, argv);
    checkAsserts();
    Shell::Shell shell;
    new RunUnrestrictedRAPTORQuery(shell);
    new RunOneToAllDijkstraRAPTORQuery(shell);
    new RunOneToManyDijkstraRAPTORQuery(shell);
    new RunTransitiveRAPTORQuery(shell);
    new RunTransitiveMcRAPTORQuery(shell);
    new RunMCRQuery(shell);
    new RunMultimodalMCRQuery(shell);
    new RunULTRARAPTORQuery(shell);
    new RunUPRAPTORQuery(shell);
    new RunULTRAMcRAPTORQuery(shell);
    new RunMultimodalULTRAMcRAPTORQuery(shell);
    new RunMultimodalUBMRAPTORQuery(shell);
    new RunMultimodalUBMHydRAQuery(shell);
    new RunHLRAPTORQuery(shell);
    new RunTransitiveTripBasedQuery(shell);
    new RunULTRATripBasedQuery(shell);
    new RunHLULTRATripBasedQuery(shell);
    new RunULTRAMcTripBasedQuery(shell);
    new RunUnrestrictedProfileRAPTORQuery(shell);
    new RunOneToAllRAPTORQuery(shell);
    new RunOneToAllRangeRAPTORQuery(shell);
    new RunUPRangeRAPTORQuery(shell);
    new RunUnrestrictedCSAQuery(shell);
    new RunTransitiveCSAQuery(shell);
    new RunTransitiveProfileCSAQuery(shell);
    new RunULTRACSAQuery(shell);
    new RunUPCSAQuery(shell);
    new RunOneToAllDijkstraCSAQuery(shell);
    new RunDijkstraQuery(shell);
    new RunCHQuery(shell);
    new RunPHASTQuery(shell);
    shell.run();
    return 0;
}
