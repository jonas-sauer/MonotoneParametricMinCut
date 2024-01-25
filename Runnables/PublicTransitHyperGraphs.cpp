#include <iostream>
#include <algorithm>
#include <random>
#include <vector>
#include <string>
#include <cmath>

#include "../Helpers/Assert.h"
#include "../Helpers/Debug.h"
#include "../Helpers/Timer.h"
#include "../Helpers/Console/CommandLineParser.h"

#include "../DataStructures/RAPTOR/Data.h"
#include "../DataStructures/RAPTOR/HyperGraph/HyperGraph.h"
#include "../DataStructures/RAPTOR/HyperGraph/Partition.h"

#include "../Visualization/PNG.h"
#include "../Visualization/PDF.h"
#include "../Visualization/SVG.h"
#include "../Visualization/Visualization.h"
#include "../Visualization/TimeTableVisualization.h"

void usage(char *program) {
    std::cout << "Usage: " << program << " <options>" << std::endl
         << std::endl
         << "-export <options>" << std::endl
         << "    -pt    <file>      -- root directory containing the public transit network files (with extension)" << std::endl
         << "    -hgr   <file>      -- Filename of the hyper graph (without extension)" << std::endl
         << "    -rw    <int>       -- Weight for the route-hyper-vertices [default: 0]" << std::endl
         << "                          0 - Defines a constant weight of 1 for all route-hyper-vertices" << std::endl
         << "                          1 - Defines a weight equal to the number of trips within the route" << std::endl
         << "                          2 - Defines a weight equal to the number of stops contained in the route" << std::endl
         << "                          3 - Defines a weight equal to the number of stop-events in the route" << std::endl
         << "                          4 - Defines a weight equal to the sum of trips and stops in the route" << std::endl
         << "    -ew    <int>       -- Defines the constant weight used for edge-hyper-vertices [default: 1]" << std::endl
         << "    -vw    <int>       -- Weight for the vertex-hyper-edges [default: 1]" << std::endl
         << "                          A positive value defines a constant weight for all hyper-edges" << std::endl
         << "                          A negative value defines a weight equal to the number of trips reachable from the vertex within 5 minutes" << std::endl
         << std::endl
         << "-visualize <options>" << std::endl
         << "    -pt    <file>      -- root directory to be served (with extension)" << std::endl
         << "    -part  <file>      -- Filename of the partition (with extension)" << std::endl
         << "    -vis   <file>      -- Filename for the visualization" << std::endl
         << "    -f     <string>    -- Format of the visualization (pdf, png, or svg) [default: pdf]" << std::endl
         << std::endl;
    exit(0);
}

template<typename DOCUMENT>
inline void draw(DOCUMENT doc, const RAPTOR::Partition& part) noexcept {
    doc.drawRoutesByType();
    doc.newPage();
    doc.drawPartition(part);
}

int main(int argc, char** argv) {
    checkAsserts();

    CommandLineParser clp(argc, argv);
    if (!(clp.isSet("export") || clp.isSet("visualize"))) usage(argv[0]);

    if (clp.isSet("export")) {
        const std::string ptFilename = clp.get<std::string>("pt", "");
        const std::string hgrFilename = clp.get<std::string>("hgr", "");
        const int rw = clp.get<int>("rw", 0);
        const int ew = clp.get<int>("ew", 1);
        const int vw = clp.get<int>("vw", 1);
        if (ptFilename == "" || hgrFilename == "") usage(argv[0]);

        RAPTOR::Data pt = RAPTOR::Data::FromBinary(ptFilename);
        pt.printInfo();
        RAPTOR::HyperGraph hgr(pt, rw, ew, vw);
        hgr.printInfo();
        hgr.toHMetis(hgrFilename);
    }

    if (clp.isSet("visualize")) {
        const std::string ptFilename = clp.get<std::string>("pt", "");
        const std::string partFilename = clp.get<std::string>("part", "");
        const std::string visFilename = clp.get<std::string>("vis", "");
        const std::string format = clp.get<std::string>("f", "pdf");
        if (ptFilename == "" || partFilename == "" || visFilename == "") usage(argv[0]);

        RAPTOR::Data pt = RAPTOR::Data::FromBinary(ptFilename);
        pt.printInfo();
        RAPTOR::HyperGraph hgr(pt);
        RAPTOR::Partition part(hgr, partFilename, true);
        if (format == "pdf") draw(TimeTableVisualization<PDF>::FromRAPTOR(visFilename, pt, 0.3), part);
        if (format == "png") draw(TimeTableVisualization<PNG>::FromRAPTOR(visFilename, pt, 0.3), part);
        if (format == "svg") draw(TimeTableVisualization<SVG>::FromRAPTOR(visFilename, pt, 0.3), part);
    }

}
