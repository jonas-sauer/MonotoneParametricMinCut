#include <iostream>
#include <vector>
#include <string>

#include "HTTPServer.h"

#include "Servlets/AlgorithmsAndParameterInformationServlet.h"
#include "Servlets/BoundingBoxServlet.h"
#include "Servlets/JSONServlet.h"
#include "Servlets/FileServlet.h"
#include "Servlets/QueryServlet.h"

#include "../DataStructures/Result.h"
#include "../DataStructures/Parameter.h"
#include "../DataStructures/Geometry/CoordinateTree.h"
#include "../DataStructures/Geometry/Point.h"

#include "../Helpers/Console/CommandLineParser.h"

#include "../Algorithms/MapServer/HelloWorld.h"

using CoordinateTreeType = CoordinateTree<Geometry::GeoMetric>;

void usage(char *program) {
    std::cout << "Usage: " << program << " <options>" << std::endl
              << std::endl
              << "-d   <directory> -- root directory to be served" << std::endl
              << "-p   <integer>   -- port number of server (Default: 8080)" << std::endl
              << "-h               -- print help" << std::endl
              << std::endl
              << "(Try " << program << " -d ../Sites/MapServer)" << std::endl
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

    // Use "server(true)" in release mode and "server(false)" for debugging!
    MapServer::HTTPServer server(true);
    server.setRootDirectory(directory);

    // Load a list of coordinates
    // Clicking on the map will select the nearest neighbor from this vector
    std::vector<Geometry::Point> coordinates;
    coordinates.push_back(Geometry::Point(Construct::LatLong, 49.013747, 8.420266));
    coordinates.push_back(Geometry::Point(Construct::LatLong, 49.012540, 8.415863));
    coordinates.push_back(Geometry::Point(Construct::LatLong, 49.011860, 8.416898));
    coordinates.push_back(Geometry::Point(Construct::LatLong, 49.013932, 8.404513));
    coordinates.push_back(Geometry::Point(Construct::LatLong, 49.006878, 8.403641));
    coordinates.push_back(Geometry::Point(Construct::LatLong, 49.009232, 8.403902));
    coordinates.push_back(Geometry::Point(Construct::LatLong, 49.015442, 8.425663));
    coordinates.push_back(Geometry::Point(Construct::LatLong, 49.009774, 8.412496));
    CoordinateTreeType ct(coordinates);

    BoundingBoxServlet bbs(coordinates);
    server.registerServlet(&bbs);

    AlgorithmsAndParameterInformationServlet ais;
    server.registerServlet(&ais);

    JSONServlet js;
    ais.registerAlgorithm(js.Handler(), "JSON Texteingabe", js.parameters);
    server.registerServlet(&js);

    // Add an algorithm to the server in 4 steps
    // 1. Create an instance of a class that inherits from class Algorithm in Algorithms/MapServer/Algorithm.h (here HelloWorld).
    // 2. Create a QueryServlet serving your algorithm (First parameter (here "halloWorldServlet") is used via HTTP GET to call the algorithm).
    // 3. Register the QueryServlet with the AlgorithmsAndParameterInformationServlet so that the algorithm can be selected within the browser (Second parameter (here "Hello World!") is displayed to the user as name of the algorithm).
    // 4. Register the QueryServlet with the Server so that the server can invoke execution of the Algorithm if requested.
    HelloWorld halloWorld;
    QueryServlet<CoordinateTreeType> halloWorldServlet("halloWorldServlet", halloWorld, coordinates, ct);
    ais.registerAlgorithm(halloWorldServlet, "Hello World!");
    server.registerServlet(&halloWorldServlet);

    server.start(port);
    do {
        getchar();
    } while(!server.stop());

    return 0;
}
