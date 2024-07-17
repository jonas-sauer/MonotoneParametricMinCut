#pragma once

#include <map>
#include <string>
#include <iostream>

#include "Servlet.h"

#include "../Responses/ResponseJSON.h"
#include "../QueryStringParser.h"
#include "../mongoose/mongoose.h"

#include "../../DataStructures/Geometry/Point.h"

class JSONServlet : public MapServer::Servlet {

public:

    JSONServlet(std::string handler = "jsonservlet") : Servlet(handler) {
        parameters["input"] = Parameter("string", "input", "JSON String", "", true);
    }

    virtual ~JSONServlet() {}

    virtual void Process(mg_connection *conn, const mg_request_info *requestInfo) {
        int maxLength = strlen(requestInfo->query_string);
        char* valueBuffer = new char[maxLength];
        int length = mg_get_var(requestInfo->query_string, maxLength, "input", valueBuffer, maxLength);
        if (length <= 0) {
            std::cout << "no result" << std::endl;
            MapServer::ResponseJSON jsonresponse("{}");
            jsonresponse.Write(conn);
        } else {
            std::cout << "returned " << length << " chars." << std::endl;
            std::string json(valueBuffer);
            MapServer::ResponseJSON jsonresponse(json);
            jsonresponse.Write(conn);
        }
        delete[] valueBuffer;
        std::cout << "The calculation has been completed." << std::endl << std::endl;
    }

    std::map<std::string, Parameter> parameters;

};
