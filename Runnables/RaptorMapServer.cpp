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

#include "../Algorithms/CSA/OldProfileCSA.h"
#include "../Algorithms/Dijkstra/Dijkstra.h"
#include "../Algorithms/RAPTOR/AlternatingRAPTOR.h"
#include "../Algorithms/RAPTOR/Profiler.h"
#include "../Algorithms/RAPTOR/RAPTOR.h"
#include "../Algorithms/RAPTOR/DijkstraRAPTOR.h"

using CoordinateTreeType = CoordinateTree<Geometry::GeoMetric>;
using DijkstraType = Dijkstra<RAPTOR::TransferGraph, true>;
using ProfileCSAType = CSA::OldProfileCSA<12>;

const std::string raptorDefault = "/algoDaten/jsauer/Networks/Switzerland/Walking/Full/raptor.binary";
const std::string transitiveRaptorDefault = "/algoDaten/jsauer/Networks/Switzerland/Walking/Transitive/raptor.binary";

inline void usage(const char *program) noexcept {
    std::cout << "Usage: " << program << " <options>" << std::endl
         << std::endl
         << "-d   <directory> -- root directory to be served (Default: ../Sites/MapServer)" << std::endl
         << "-p   <integer>   -- port number of server (Default: 8080)" << std::endl
         << "-r   <file>      -- RAPTOR binary (Default: " << raptorDefault << ")" << std::endl
         << "-t   <file>      -- transitive RAPTOR binary (Default: " << transitiveRaptorDefault << ")" << std::endl
         << "-c   <file>      -- CSA binary (Optional)" << std::endl
         << "-nobuf           -- don't use departure buffer times" << std::endl
         << "-h               -- print help" << std::endl
         << std::endl;
    exit(0);
}

template<bool USE_MINIMUM_TRANSFER_TIMES>
inline void setup(CommandLineParser& clp) noexcept {
    using DijkstraRaptorType = RAPTOR::DijkstraRAPTOR<RAPTOR::DijkstraInitialTransfers, RAPTOR::SimpleProfiler, true, USE_MINIMUM_TRANSFER_TIMES, false>;
    using DijkstraRaptorTypeNoTP = RAPTOR::DijkstraRAPTOR<RAPTOR::DijkstraInitialTransfers, RAPTOR::SimpleProfiler, false, USE_MINIMUM_TRANSFER_TIMES, false>;
    using TransitiveRaptorType = RAPTOR::RAPTOR<true, RAPTOR::SimpleProfiler, true, USE_MINIMUM_TRANSFER_TIMES, false>;
    using ProfileRaptorType = RAPTOR::AlternatingDijkstraRAPTOR<RAPTOR::DijkstraInitialTransfers, true, USE_MINIMUM_TRANSFER_TIMES>;

    const int port(clp.value<int>("p", 8080));
    const std::string directory(clp.value<std::string>("d", "../Sites/MapServer"));
    const std::string raptorFile(clp.value<std::string>("r", raptorDefault));
    const std::string transitiveRaptorFile(clp.value<std::string>("t", transitiveRaptorDefault));

    MapServer::HTTPServer server(false);
    server.setRootDirectory(directory);

    RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(raptorFile);
    if constexpr (!USE_MINIMUM_TRANSFER_TIMES) raptorData.useImplicitDepartureBufferTimes();
    RAPTOR::Data reverseData = raptorData.reverseNetwork();
    raptorData.printInfo();
    std::cout << std::endl;

    RAPTOR::TransferGraph graph = raptorData.transferGraph;
    Graph::printInfo(graph);
    graph.printAnalysis();
    std::cout << std::endl;

    RAPTOR::Data transitiveRaptorData = RAPTOR::Data::FromBinary(transitiveRaptorFile);
    if constexpr (!USE_MINIMUM_TRANSFER_TIMES) transitiveRaptorData.useImplicitDepartureBufferTimes();
    transitiveRaptorData.printInfo();
    std::cout << std::endl;
    Graph::printInfo(transitiveRaptorData.transferGraph);
    transitiveRaptorData.transferGraph.printAnalysis();
    std::cout << std::endl;

    DijkstraType dijkstra(graph, graph[TravelTime]);

    std::vector<Geometry::Point> coordinates = graph[Coordinates];
    CoordinateTreeType ct(coordinates);
    std::vector<Geometry::Point> stopCoordinates = coordinates;
    stopCoordinates.resize(raptorData.numberOfStops());
    CoordinateTreeType sct(stopCoordinates);

    DijkstraRaptorTypeNoTP raptorNoTP(raptorData, raptorData.transferGraph, reverseData.transferGraph);
    DijkstraRaptorType raptorTP(raptorData, raptorData.transferGraph, reverseData.transferGraph);
    TransitiveRaptorType transitiveRaptor(transitiveRaptorData);
    ProfileRaptorType profileRaptor(raptorData, reverseData);

    BoundingBoxServlet bbs(coordinates);
    server.registerServlet(&bbs);

    AlgorithmsAndParameterInformationServlet ais;
    server.registerServlet(&ais);

    JSONServlet js;
    ais.registerAlgorithm(js.Handler(), "JSON Texteingabe", js.parameters);
    server.registerServlet(&js);

    ShortestPath<DijkstraType, true> shortestPathDijkstra(dijkstra);
    QueryServletWithVertexInput<CoordinateTreeType> dijkstraServlet("dijkstraServlet", shortestPathDijkstra, coordinates, ct);
    ais.registerAlgorithm(dijkstraServlet, "Dijkstra");
    server.registerServlet(&dijkstraServlet);

    RaptorQuery<DijkstraRaptorType> raptorTPQuery(raptorTP, raptorData);
    QueryServletWithVertexInput<CoordinateTreeType> raptorTPServlet("raptorTPServlet", raptorTPQuery, coordinates, ct);
    ais.registerAlgorithm(raptorTPServlet, "RAPTOR");
    server.registerServlet(&raptorTPServlet);

    RaptorQuery<DijkstraRaptorTypeNoTP> raptorNoTPQuery(raptorNoTP, raptorData);
    QueryServletWithVertexInput<CoordinateTreeType> raptorNoTPServlet("raptorNoTPServlet", raptorNoTPQuery, coordinates, ct);
    ais.registerAlgorithm(raptorNoTPServlet, "RAPTOR (No Target Pruning)");
    server.registerServlet(&raptorNoTPServlet);

    RaptorQuery<TransitiveRaptorType, true> transitiveRaptorQuery(transitiveRaptor, raptorData);
    QueryServletWithVertexInput<CoordinateTreeType> transitiveRaptorServlet("transitiveRaptorServlet", transitiveRaptorQuery, stopCoordinates, sct);
    ais.registerAlgorithm(transitiveRaptorServlet, "Transitive RAPTOR");
    server.registerServlet(&transitiveRaptorServlet);

    ProfileRaptorQuery<ProfileRaptorType> profileRaptorQuery(profileRaptor, raptorData);
    QueryServletWithVertexInput<CoordinateTreeType> profileRaptorServlet("profileRaptorServlet", profileRaptorQuery, coordinates, ct);
    ais.registerAlgorithm(profileRaptorServlet, "Profile RAPTOR");
    server.registerServlet(&profileRaptorServlet);

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

    ShowRoute sr(raptorData);
    QueryServlet<CoordinateTreeType> srServlet("showRoute", sr, stopCoordinates, sct);
    ais.registerAlgorithm(srServlet, "Show Route");
    server.registerServlet(&srServlet);

    ShowRoutes srs(raptorData);
    QueryServlet<CoordinateTreeType> srsServlet("showRoutes", srs, stopCoordinates, sct);
    ais.registerAlgorithm(srsServlet, "Show Routes");
    server.registerServlet(&srsServlet);

    IsolatedStops is(raptorData);
    QueryServlet<CoordinateTreeType> isServlet("isolatedStops", is, stopCoordinates, sct);
    ais.registerAlgorithm(isServlet, "Isolated Stops");
    server.registerServlet(&isServlet);

    HelloWorld hw;
    QueryServlet<CoordinateTreeType> hwServlet("hwServlet", hw, coordinates, ct);
    ais.registerAlgorithm(hwServlet, "Hello World!");
    server.registerServlet(&hwServlet);

    if (clp.isSet("c")) {
        const std::string csaFile(clp.value<std::string>("c", ""));

        CSA::Data csaData = CSA::Data::FromBinary(csaFile);
        csaData.sortConnectionsAscendingByDepartureTime();
        csaData.printInfo();
        std::cout << std::endl;

        CSA::TransferGraph reverseGraph = csaData.transferGraph;
        reverseGraph.revert();
        Graph::printInfo(reverseGraph);
        reverseGraph.printAnalysis();
        std::cout << std::endl;

        ProfileCSAType profileCSA(csaData, reverseGraph);

        ProfileRaptorQuery<ProfileCSAType, true> profileCSAQuery(profileCSA, raptorData);
        QueryServletWithVertexInput<CoordinateTreeType> profileCSAServlet("profileCsaServlet", profileCSAQuery, coordinates, ct);
        ais.registerAlgorithm(profileCSAServlet, "Profile CSA");
        server.registerServlet(&profileCSAServlet);

        server.start(port);
        do {
            getchar();
        } while(!server.stop());
    } else {
        server.start(port);
        do {
            getchar();
        } while(!server.stop());
    }
}

int main(int argc, char **argv) {
    checkAsserts();
    CommandLineParser clp(argc, argv);
    if (clp.isSet("h")) usage(argv[0]);
    if (clp.isSet("nobuf")) {
        setup<true>(clp);
    } else {
        setup<false>(clp);
    }
    return 0;
}
