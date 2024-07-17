#include <sched.h>

#include "../Helpers/Console/CommandLineParser.h"
#include "../Shell/Shell.h"

#include "Commands/CH.h"

int main(int argc, char** argv) {
    CommandLineParser clp(argc, argv);
    pinThreadToCoreId(clp.value<int>("core", 1));
    checkAsserts();
    ::Shell::Shell shell;
    new BuildCH(shell);
    new ResumeCH(shell);
    new BuildCHFromOrder(shell);
    new BuildCoreCH(shell);
    new BuildPartialKeyCH(shell);
    new BuildMinLevelKeyCH(shell);
    new BuildFactorKeyCH(shell);
    new RunCHQueries(shell);
    new RunBucketCHQueries(shell);
    new RunBidirectionalRPHASTQueries(shell);
    new BuildBlockingKeyCH(shell);
    shell.run();
    return 0;
}
