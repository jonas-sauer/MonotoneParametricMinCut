#include <sched.h>

#include "../Helpers/Console/CommandLineParser.h"
#include "../Shell/Shell.h"
#include "Commands/UnrestrictedWalkingExperiments.h"

using namespace Shell;

int main(int argc, char** argv) {
    CommandLineParser clp(argc, argv);
    checkAsserts();
    Shell::Shell shell;
    new GenerateRankQueries(shell);
    new GenerateConnectionBasedRankQueries(shell);
    new GenerateConnectionBasedRandomQueries(shell);
    new CompareTravelTimes(shell);
    new AnalyzeTravelTimeDifference(shell);
    new FindUnrestrictedWalkingExamples(shell);
    shell.run();
    return 0;
}
