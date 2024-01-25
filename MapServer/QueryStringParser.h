#pragma once

#include <string>
#include <cstring>

#include "mongoose/mongoose.h"

#include "../Helpers/String/String.h"

namespace MapServer {

class QueryStringParser {

public:
    QueryStringParser(const mg_request_info *ri) : requestInfo(ri) {}

    template<typename T>
    inline T getValue(const std::string key, const T defaultValue = T()) {
        int length = mg_get_var(requestInfo->query_string, strlen(requestInfo->query_string), key.c_str(), valueBuffer, sizeof(valueBuffer));
        if (length == -1) return defaultValue;
        std::string value(valueBuffer);
        return String::lexicalCast<T>(value);
    }

private:
    const mg_request_info *requestInfo;
    char valueBuffer[1024];

};

}
