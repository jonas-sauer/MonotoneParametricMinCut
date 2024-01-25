#include "Commands/Evaluation.h"

#include "../Helpers/MultiThreading.h"
#include "../Helpers/Console/CommandLineParser.h"

#include "../Shell/Shell.h"
using namespace Shell;

int main(int argc, char** argv) {
    CommandLineParser clp(argc, argv);
    pinThreadToCoreId(clp.value<int>("core", 1));
    checkAsserts();
    ::Shell::Shell shell;

    new TikzPlot(shell);
    new AccumulateDijkstraRank(shell);
    new DijkstraRankPlot(shell);
    new BoxCode(shell);
    new AnalyzeData(shell);
    new AnalyzeHeuristics(shell);
    new BuildDfgScript(shell);

    shell.run();
    return 0;
}
