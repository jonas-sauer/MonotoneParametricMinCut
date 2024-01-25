#include <sched.h>

#include "../Helpers/Console/CommandLineParser.h"
#include "../Shell/Shell.h"
#include "Commands/QueryValidation.h"

using namespace Shell;

int main(int argc, char** argv) {
    CommandLineParser clp(argc, argv);
    checkAsserts();
    Shell::Shell shell;
    new ValidateOneToAllDijkstraRAPTOR(shell);
    new ValidateOneToManyDijkstraRAPTOR(shell);
    new ValidateULTRARAPTOR(shell);
    new ValidateUPRAPTOR(shell);
    new ValidateHL(shell);
    new ValidateHLRAPTOR(shell);
    new ValidateULTRATripBased(shell);
    new ValidateHLULTRATripBased(shell);
    new ValidateTransitiveProfileRAPTOR(shell);
    new ValidateDijkstraProfileRAPTOR(shell);
    new ValidateULTRAProfileRAPTOR(shell);
    new ValidateBoundedMcRAPTOR(shell);
    new ValidateUBMRAPTOR(shell);
    new ValidateNewUBMRAPTOR(shell);
    new ValidateUBMTBRAPTOR(shell);
    new ValidateUBMHydRA(shell);
    new ValidateULTRAMcRAPTOR(shell);
    new ValidateBoundedULTRAMcTripBased(shell);
    new ValidateULTRAMcTripBased(shell);
    new ValidateTransitiveCSA(shell);
    new ValidateTransitiveProfileCSA(shell);
    new ValidateParetoCSA(shell);
    new ValidateDijkstraCSA(shell);
    new ValidateOneToAllDijkstraCSA(shell);
    new ValidateDijkstraProfileCSA(shell);
    new ValidateOneToAllDijkstraProfileCSA(shell);
    new ValidateULTRACSA(shell);
    new ValidateUPCSA(shell);
    new ValidateHLCSA(shell);
    new ValidateULTRAProfileCSA(shell);
    new ValidateUPProfileCSA(shell);
    new ValidateParetoULTRACSA(shell);
    new ValidateParetoUPCSA(shell);
    new ValidateTransitiveTripBased(shell);
    new ValidateUPTB(shell);
    new ComputeBoundedMcRAPTORCoverage(shell);
    shell.run();
    return 0;
}
