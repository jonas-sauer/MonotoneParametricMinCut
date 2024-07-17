#pragma once

#include <map>
#include <string>
#include <iostream>

#include "Servlet.h"

#include "../Responses/ResponseJSON.h"
#include "../QueryStringParser.h"
#include "..//mongoose/mongoose.h"

#include "../../DataStructures/Geometry/Point.h"

class FileServlet : public MapServer::Servlet {

public:

    FileServlet(std::string handler = "fileservlet") : Servlet(handler) {
        parameters["filename"] = Parameter("string", "filename", "JSON Filename", "", true);
    }

    virtual ~FileServlet() {}

    virtual void Process(mg_connection *conn, const mg_request_info *requestInfo) {
        MapServer::QueryStringParser qsp(requestInfo);
        std::string filename = qsp.getValue<std::string>("filename", "");
        if (filename.find("file://") == 0) filename = filename.substr(7, filename.size() - 7);
        std::cout << "Requesting file: " << filename << std::endl;
        std::string json("");
        if (filename != "") {
            std::ifstream in(filename);
            if (in) {
                std::cout << "File found!" << std::endl;
                in.seekg(0, std::ios::end);
                json.resize(in.tellg());
                in.seekg(0, std::ios::beg);
                in.read(&json[0], json.size());
                in.close();
            }
        } else {
            json = "{}";
        }
        std::cout << "Result:" << std::endl;
        if (json.size() <= 200) {
            std::cout << json << std::endl;
        } else {
            std::cout << json.substr(0, 199) << "... (" << std::to_string(json.size() - 200) << " more characters)" << std::endl;
        }
        MapServer::ResponseJSON jsonresponse(json);
        jsonresponse.Write(conn);
        std::cout << "The calculation has been completed." << std::endl << std::endl;
    }

    std::map<std::string, Parameter> parameters;

};
