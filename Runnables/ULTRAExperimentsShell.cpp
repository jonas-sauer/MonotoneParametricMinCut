#include <sched.h>

#include "../Helpers/Debug.h"
#include "../Helpers/Console/CommandLineParser.h"
#include "../Shell/Shell.h"
#include "Commands/ULTRAExperiments.h"
#include "Commands/ULTRAPreprocessing.h"

using namespace Shell;

int main(int argc, char** argv) {
    CommandLineParser clp(argc, argv);
    checkAsserts();
    Shell::Shell shell;
    new ComputeStopToStopShortcuts(shell);
    new CoreDegreeWitnessLimitExperiment(shell);
    new ParallelizationExperiment(shell);
    new TransferSpeedPreprocessingExperiment(shell);
    new TransferSpeedMcPreprocessingExperiment(shell);
    new RAPTORQueryExperiment(shell);
    new CSAQueryExperiment(shell);
    new TBQueryExperiment(shell);
    new TransferSpeedQueryExperiment(shell);
    new TransferSpeedMcQueryExperiment(shell);
    new TransferSpeedTravelTimeExperiment(shell);
    new GeoRankExperiment(shell);
    new EdgeLengthHistogram(shell);
    new EdgeGeoDistanceHistogram(shell);
    new StopDegreeHistogram(shell);
    new ComputeSizeOfTransitiveClosure(shell);
    new ComputeTransitiveClosure(shell);
    new CompareGraphs(shell);
    shell.run();
    return 0;
}
