#include <iostream>
#include <vector>
#include <string>

#include "../Helpers/Debug.h"
#include "../Helpers/MultiThreading.h"
#include "../Helpers/Console/CommandLineParser.h"

#include "Commands/CH.h"
#include "Commands/Graphs.h"
#include "Commands/NetworkAnalysis.h"
#include "Commands/NetworkIO.h"
#include "Commands/NetworkTools.h"
#include "Commands/Visualization.h"

#include "../Shell/Shell.h"
using namespace Shell;

int main(int argc, char** argv) {
    CommandLineParser clp(argc, argv);
    pinThreadToCoreId(clp.value<int>("core", 1));
    checkAsserts();
    ::Shell::Shell shell;
    new LoadDimacsGraph(shell);
    new DynamicToStaticTransferGraph(shell);
    new ReverseTransferGraph(shell);
    new BenToDimacs(shell);
    new SaveDimacsGraph(shell);

    new ParseGTFS(shell);
    new GTFSToIntermediate(shell);
    new IntermediateToCSA(shell);
    new IntermediateToRAPTOR(shell);
    new ParseCSAFromCSV(shell);
    new CSAToCSV(shell);
    new CSAToExchangeFormat(shell);
    new CSAToIntermediate(shell);
    new RAPTORToCSV(shell);
    new RAPTORToIntermediate(shell);
    new BuildMultimodalRAPTORData(shell);
    new AddModeToMultimodalRAPTORData(shell);
    new BuildMultimodalTripBasedData(shell);
    new AddModeToMultimodalTripBasedData(shell);

    new DuplicateTrips(shell);
    new ReverseRAPTORNetwork(shell);
    new AddGraph(shell);
    new ReplaceGraph(shell);
    new ReduceGraph(shell);
    new ReduceToMaximumConnectedComponent(shell);
    new ReduceToMaximumConnectedComponentWithTransitive(shell);
    new ApplyBoundingBox(shell);
    new MakeOneHopTransfers(shell);
    new MakeOneHopTransfersByGeoDistance(shell);
    new MakeTransitiveStopGraph(shell);
    new ApplyMaxTransferSpeed(shell);
    new ApplyConstantTransferSpeed(shell);
    new ApplyMinTransferTravelTime(shell);

    new AnalyzeNetwork(shell);
    new AnalyzeConnectivity(shell);
    new AnalyzeConnectedComponents(shell);
    new AnalyzeMinTransferTimes(shell);
    new PrintRoute(shell);
    new PrintTrip(shell);

    new DrawRAPTOR(shell);
    new DrawCSA(shell);

    new BuildCH(shell);
    new BuildCoreCH(shell);
    shell.run();
    return 0;
}
