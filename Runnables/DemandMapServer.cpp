#include <iostream>
#include <vector>
#include <string>

#include "../MapServer/HTTPServer.h"
#include "../MapServer/Servlets/AlgorithmsAndParameterInformationServlet.h"
#include "../MapServer/Servlets/BoundingBoxServlet.h"
#include "../MapServer/Servlets/JSONServlet.h"
#include "../MapServer/Servlets/FileServlet.h"
#include "../MapServer/Servlets/QueryServlet.h"
#include "../MapServer/Servlets/QueryServletWithVertexInput.h"

#include "../DataStructures/Result.h"
#include "../DataStructures/Parameter.h"
#include "../DataStructures/Geometry/CoordinateTree.h"
#include "../DataStructures/Geometry/Metric.h"
#include "../DataStructures/Geometry/Point.h"

#include "../DataStructures/Graph/Graph.h"

#include "../DataStructures/CSA/Data.h"
#include "../DataStructures/Demand/PassengerData.h"

#include "../Helpers/Console/CommandLineParser.h"

#include "../Algorithms/MapServer/HelloWorld.h"
#include "../Algorithms/MapServer/FindVertex.h"
#include "../Algorithms/MapServer/ShowGraph.h"
#include "../Algorithms/MapServer/RaptorQuery.h"
#include "../Algorithms/MapServer/ProfileRaptorQuery.h"
#include "../Algorithms/MapServer/IsolatedStops.h"
#include "../Algorithms/MapServer/ShowAssignedPath.h"

using CoordinateTreeType = CoordinateTree<Geometry::GeoMetric>;

constexpr int defaultPort = 8080;
constexpr char defaultCSA[] = "/home/tobias/Data/stuttgart/Transitive/csa.binary";
constexpr char defaultAssignment[] = "/home/tobias/Data/stuttgart/Transitive/assignment.binary";

void usage(char *program) {
    std::cout << "Usage: " << program << " <options>" << std::endl
         << std::endl
         << "-d   <directory> -- root directory to be served" << std::endl
         << "-p   <integer>   -- port number of server (Default: " << defaultPort << ")" << std::endl
         << "-c   <file>      -- CSA binary (Default: " << defaultCSA << ")" << std::endl
         << "-a   <file>      -- Assignment binary (Default: /" << defaultAssignment << ")" << std::endl
         << "-h               -- print help" << std::endl
         << std::endl;
    exit(0);
}

int main(int argc, char **argv) {
    checkAsserts();

    CommandLineParser clp(argc, argv);
    if (clp.isSet("h")) usage(argv[0]);

    const int port(clp.value<int>("p", defaultPort));
    const std::string directory(clp.value<std::string>("d", ""));
    if (directory == "") usage(argv[0]);

    const std::string csaFile(clp.value<std::string>("c", defaultCSA));
    const std::string assignmentFile(clp.value<std::string>("a", defaultAssignment));

    MapServer::HTTPServer server(false);
    server.setRootDirectory(directory);

    CSA::Data csaData = CSA::Data::FromBinary(csaFile);
    csaData.sortConnectionsAscendingByDepartureTime();
    csaData.printInfo();
    std::cout << std::endl;

    CSA::TransferGraph graph = csaData.transferGraph;
    Graph::printInfo(graph);
    graph.printAnalysis();
    std::cout << std::endl;

    PassengerData assignmentData = PassengerData::FromBinary(assignmentFile);
    assignmentData.printInfo();
    std::cout << std::endl;

    std::vector<Geometry::Point> coordinates = graph[Coordinates];
    CoordinateTreeType ct(coordinates);

    BoundingBoxServlet bbs(coordinates);
    server.registerServlet(&bbs);

    AlgorithmsAndParameterInformationServlet ais;
    server.registerServlet(&ais);

    JSONServlet js;
    ais.registerAlgorithm(js.Handler(), "JSON Texteingabe", js.parameters);
    server.registerServlet(&js);

    ShowAssignedPath sap(csaData, assignmentData);
    QueryServlet<CoordinateTreeType> sapServlet("showAssignedPath", sap, coordinates, ct);
    ais.registerAlgorithm(sapServlet, "Show Assigned Path");
    server.registerServlet(&sapServlet);

    FindVertex fv;
    QueryServlet<CoordinateTreeType> fvServlet("findVertex", fv, coordinates, ct);
    ais.registerAlgorithm(fvServlet, "Find Vertex");
    server.registerServlet(&fvServlet);

    ShowGraph<CSA::TransferGraph> sg(graph);
    QueryServlet<CoordinateTreeType> sgServlet("showGraph", sg, coordinates, ct);
    ais.registerAlgorithm(sgServlet, "Show Graph");
    server.registerServlet(&sgServlet);

    HelloWorld hw;
    QueryServlet<CoordinateTreeType> hwServlet("hwServlet", hw, coordinates, ct);
    ais.registerAlgorithm(hwServlet, "Hello World!");
    server.registerServlet(&hwServlet);

    server.start(port);
    do {
        getchar();
    } while(!server.stop());

    return 0;
}
