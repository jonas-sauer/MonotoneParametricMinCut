#pragma once

#include <string>

#include "Response.h"

namespace MapServer {

class ResponseText : public Response {

public:
    ResponseText(const std::string contents) : Response("text/plain", contents) {}

};

}
