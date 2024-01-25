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

#include "../Helpers/Console/CommandLineParser.h"

#include "../Algorithms/MapServer/ShortestPath.h"
#include "../Algorithms/MapServer/HelloWorld.h"
#include "../Algorithms/MapServer/FindVertex.h"
#include "../Algorithms/MapServer/ShowGraph.h"
#include "../Algorithms/MapServer/ShowEdges.h"
#include "../Algorithms/MapServer/ShowRoute.h"
#include "../Algorithms/MapServer/RaptorQuery.h"
#include "../Algorithms/MapServer/ProfileRaptorQuery.h"
#include "../Algorithms/MapServer/IsolatedStops.h"

#include "../Algorithms/Dijkstra/Dijkstra.h"

using CoordinateTreeType = CoordinateTree<Geometry::GeoMetric>;

void usage(char *program) {
    std::cout << "Usage: " << program << " <options>" << std::endl
         << std::endl
         << "-d                  <directory> -- root directory to be served" << std::endl
         << "-p                  <integer>   -- port number of server (Default: 8080)" << std::endl
         << "-ben                <file>      -- graph to be used in Ben's format" << std::endl
         << "-dimacs             <file>      -- graph to be used in dimacs format" << std::endl
         << "-binary             <file>      -- graph to be used in binary format" << std::endl
         << "-secondsPerTimeUnit <double>    -- conversion factor from time unit to seconds" << std::endl
         << "-h                              -- print help" << std::endl
         << std::endl;
    exit(0);
}

int main(int argc, char **argv) {
    checkAsserts();

    CommandLineParser clp(argc, argv);
    if (clp.isSet("h")) usage(argv[0]);

    const int port(clp.value<int>("p", 8080));
    const std::string directory(clp.value<std::string>("d", ""));
    if (directory == "") usage(argv[0]);

    const std::string benFile(clp.value<std::string>("ben", ""));
    const std::string dimacsFile(clp.value<std::string>("dimacs", ""));
    const std::string binaryFile(clp.value<std::string>("binary", ""));
    const double secondsPerTimeUnit = clp.value<double>("secondsPerTimeUnit", 1);

    MapServer::HTTPServer server(false);
    server.setRootDirectory(directory);

    TransferGraph graph;
    if (benFile != "") Graph::fromStrasserBinary(benFile, graph);
    if (dimacsFile != "") Graph::fromDimacs(dimacsFile, graph);
    if (binaryFile != "") graph.readBinary(binaryFile);
    Graph::printInfo(graph);
    graph.printAnalysis();
    std::cout << std::endl;

    Dijkstra<RAPTOR::TransferGraph, true> dijkstra(graph, graph[TravelTime]);

    std::vector<Geometry::Point> coordinates = graph[Coordinates];
    CoordinateTreeType ct(coordinates);

    BoundingBoxServlet bbs(coordinates);
    server.registerServlet(&bbs);

    AlgorithmsAndParameterInformationServlet ais;
    server.registerServlet(&ais);

    JSONServlet js;
    ais.registerAlgorithm(js.Handler(), "JSON Texteingabe", js.parameters);
    server.registerServlet(&js);

    ShortestPath<Dijkstra<RAPTOR::TransferGraph, true>, true> shortestPathDijkstra(dijkstra, secondsPerTimeUnit);
    QueryServletWithVertexInput<CoordinateTreeType> dijkstraServlet("dijkstraServlet", shortestPathDijkstra, coordinates, ct);
    ais.registerAlgorithm(dijkstraServlet, "Dijkstra");
    server.registerServlet(&dijkstraServlet);

    FindVertex fv;
    QueryServlet<CoordinateTreeType> fvServlet("findVertex", fv, coordinates, ct);
    ais.registerAlgorithm(fvServlet, "Find Vertex");
    server.registerServlet(&fvServlet);

    ShowGraph<RAPTOR::TransferGraph> sg(graph);
    QueryServlet<CoordinateTreeType> sgServlet("showGraph", sg, coordinates, ct);
    ais.registerAlgorithm(sgServlet, "Show Graph");
    server.registerServlet(&sgServlet);

    ShowEdges<RAPTOR::TransferGraph> se(graph);
    QueryServlet<CoordinateTreeType> seServlet("showEdges", se, coordinates, ct);
    ais.registerAlgorithm(seServlet, "Show Edges");
    server.registerServlet(&seServlet);

    server.start(port);
    do {
        getchar();
    } while(!server.stop());

    return 0;
}
