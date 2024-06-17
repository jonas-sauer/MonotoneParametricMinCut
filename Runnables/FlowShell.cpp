#include <sched.h>

#include "../Helpers/Console/CommandLineParser.h"
#include "../Shell/Shell.h"
#include "Commands/Flow.h"

using namespace Shell;

int main(int argc, char** argv) {
    CommandLineParser clp(argc, argv);
    checkAsserts();
    Shell::Shell shell;
    new LoadMaxFlowInstanceFromDimacs(shell);
    new RunPushRelabel(shell);
    new RunIBFS(shell);
    shell.run();
    return 0;
}
