#pragma once

#include <string>
#include <sstream>
#include <map>

#include "Servlet.h"
#include "QueryServlet.h"
#include "AlgorithmsInformationServlet.h"

#include "../Responses/ResponseJSON.h"

#include "../../DataStructures/Parameter.h"

class AlgorithmsAndParameterInformationServlet : public MapServer::AlgorithmsInformationServlet {
public:

    AlgorithmsAndParameterInformationServlet() : AlgorithmsInformationServlet() {}

    virtual ~AlgorithmsAndParameterInformationServlet() {}

    virtual void Process(mg_connection* conn, const mg_request_info*) {
        std::cout << "Dumping information about algorithms... " << std::flush;

        std::stringstream json;
        json << "{" << std::endl
             << "\"NumberOfAlgorithms\":" << algorithms.size() << "," << std::endl
             << "\"Algorithms\": [" << std::endl;
        for (std::map<std::string, std::string>::iterator iterator = algorithms.begin(); iterator != algorithms.end(); iterator++) {
            if (iterator != algorithms.begin()) json << ", ";
            std::map<std::string, const std::map<std::string, Parameter>*>::iterator parameters = algorithmsParameters.find(iterator->first);
            json << "{\"Key\":\"" << iterator->first << "\", \"Name\":\"" << iterator->second << "\"";
            if (parameters != algorithmsParameters.end()) {
                json << ", \"Parameter\":[";
                for (std::map<std::string, Parameter>::const_iterator parameter = parameters->second->begin(); parameter != parameters->second->end(); parameter++) {
                    if (parameter != parameters->second->begin()) json << ", ";
                    parameter->second.writeParameter(json);
                }
                json << "]";
            }
            json << "}";
        }
        json << "]}" << std::endl;

        std::cout << "done (" << json.str().length() << " bytes)" << std::endl;

        MapServer::ResponseJSON(json.str()).Write(conn);
    }

    template<typename ALGORITHM_TYPE>
    inline void registerAlgorithm(const ALGORITHM_TYPE& hls, const std::string name) {
        registerAlgorithm(hls.Handler(), name, hls.getParameters());
    }

    inline void registerAlgorithm(const std::string& key, const std::string& name) {
        algorithms[key] = name;
    }

    inline void registerAlgorithm(const std::string& key, const std::string& name, const std::map<std::string, Parameter>& parameters) {
        algorithms[key] = name;
        algorithmsParameters[key] = &parameters;
    }

protected:

    std::map<std::string, const std::map<std::string, Parameter>*> algorithmsParameters;

};
