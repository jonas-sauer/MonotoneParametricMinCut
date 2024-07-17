#include "Commands/BikeSharing.h"

#include "../Helpers/Console/CommandLineParser.h"

#include "../Shell/Shell.h"
using namespace Shell;

int main(int argc, char** argv) {
    CommandLineParser clp(argc, argv);
    pinThreadToCoreId(clp.value<int>("core", 1));
    checkAsserts();
    ::Shell::Shell shell;
    new SanitizeBikeSharingStations(shell);
    new BuildBikeSharingData(shell);
    new BuildExtendedBikeSharingData(shell);
    new BuildPartialBikeSharingData(shell);
    new GenerateBikeSharingQueries(shell);
    new GenerateBikeSharingRankQueries(shell);
    new RunBikeSharingQueries(shell);
    shell.run();
    return 0;
}