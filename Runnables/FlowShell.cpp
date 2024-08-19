#include <sched.h>

#include "../Helpers/Console/CommandLineParser.h"
#include "../Shell/Shell.h"
#include "Commands/Experiments.h"

using namespace Shell;

int main(int argc, char** argv) {
    CommandLineParser clp(argc, argv);
    checkAsserts();
    Shell::Shell shell;
    new RunPushRelabel(shell);
    new RunIBFS(shell);
    new RunExcessesIBFS(shell);
    new RunParametricIBFS(shell);
    new RunChordScheme(shell);
    new PrecisionExperiment(shell);
    shell.run();
    return 0;
}
