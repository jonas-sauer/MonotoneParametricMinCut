#include "Commands/NetworkAnalysis.h"
#include "Commands/ULTRAPreprocessing.h"

#include "../Helpers/Console/CommandLineParser.h"

#include "../Shell/Shell.h"
using namespace Shell;

int main(int argc, char** argv) {
    CommandLineParser clp(argc, argv);
    pinThreadToCoreId(clp.value<int>("core", 1));
    checkAsserts();
    ::Shell::Shell shell;
    new BuildFreeTransferGraph(shell);
    new ComputeStopToStopShortcuts(shell);
    new ComputeMcStopToStopShortcuts(shell);
    new ComputeMultimodalMcStopToStopShortcuts(shell);
    new RAPTORToTripBased(shell);
    new ComputeEventToEventShortcuts(shell);
    new ComputeTransitiveEventToEventShortcuts(shell);
    new ComputeDelayEventToEventShortcuts(shell);
    new ComputeMcEventToEventShortcuts(shell);
    new ComputeMultimodalMcEventToEventShortcuts(shell);
    new AugmentTripBasedShortcuts(shell);
    new ValidateStopToStopShortcuts(shell);
    new ValidateEventToEventShortcuts(shell);

    new PrintStopEventID(shell);
    new PrintOutgoingTransfers(shell);

    shell.run();
    return 0;
}
