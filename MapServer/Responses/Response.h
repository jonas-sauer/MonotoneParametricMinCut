#pragma once

#include <string>

#include "../mongoose/mongoose.h"

namespace MapServer {

class Response {

public:
    Response(const std::string mt = "text/plain", const std::string contents = "") : mimeType(mt), buffer(contents) {}
    virtual ~Response() {}

    virtual void Write(mg_connection *conn) {
        mg_printf(conn, "HTTP/1.1 200 OK\r\n");
        mg_printf(conn, "Content-Type: %s\r\n", mimeType.c_str());
        mg_printf(conn, "Content-Length: %d\r\n\r\n", buffer.length());
        mg_write(conn, buffer.c_str(), buffer.length());
    }

    virtual void SetContents(const std::string &contents) {
        buffer = contents;
    }

protected:
    const std::string mimeType;
    std::string buffer;

};

}
