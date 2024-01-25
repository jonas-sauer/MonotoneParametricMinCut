#pragma once

#include <string>
#include <sstream>
#include <map>

#include "Servlet.h"
#include "../Responses/ResponseJSON.h"

namespace MapServer {

class AlgorithmsInformationServlet : public Servlet {

public:
    AlgorithmsInformationServlet() : Servlet("algorithms_information") {}
    virtual ~AlgorithmsInformationServlet() {}

    virtual void Process(mg_connection* conn, const mg_request_info*) {
        std::cout << "Dumping information about algorithms... " << std::flush;

        std::stringstream json;
        json << "{" << std::endl
             << "\"NumberOfAlgorithms\":" << algorithms.size() << "," << std::endl
             << "\"Algorithms\": [" << std::endl;
        for (auto iterator = algorithms.begin(); iterator != algorithms.end(); ++iterator) {
            if (iterator != algorithms.begin())
                json << ", ";
            json << "{\"Key\":\"" << iterator->first << "\", \"Name\":\"" << iterator->second << "\"}";
        }
        json << "]}" << std::endl;

        std::cout << "done (" << json.str().length() << " bytes)" << std::endl;
        ResponseJSON(json.str()).Write(conn);
    }

    inline void registerAlgorithm(const std::string key, const std::string name) {
        algorithms[key] = name;
    }

protected:
    std::map<std::string, std::string> algorithms;

};

}
