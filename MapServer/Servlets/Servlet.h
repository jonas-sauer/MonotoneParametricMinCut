#pragma once

#include <string>

#include "../mongoose/mongoose.h"

namespace MapServer {

class Servlet {

public:
    Servlet(const std::string h) : handler(h) {}
    virtual ~Servlet() {}

    virtual void Process(mg_connection *conn, const mg_request_info *requestInfo) = 0;

    virtual std::string Handler() const {return handler; }

protected:
    std::string handler;

};

}
