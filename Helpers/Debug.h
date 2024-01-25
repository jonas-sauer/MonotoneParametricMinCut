#pragma once

#include <iostream>
#include <vector>
#include <string>

#include <unistd.h>
#include <stdlib.h>
#include <execinfo.h>
#include <cxxabi.h>
#include <signal.h>
#include <stdio.h>

#include "Assert.h"
#include "Types.h"
#include "String/String.h"

#undef AssertMsg
#define AssertMsg(assumption, msg) assert((assumption) || (std::cout << "\n\033[31mASSERTION FAILED: " << msg << "\033[0m\nFile: " << __FILE__ << "\nLine: " << __LINE__ << "\n" << std::flush && printStackTrace()))

#define __Print_Here__ std::cout << "\033[1;36mFile: " << __FILE__ << ", Line: " << __LINE__ << "\033[0m" << std::endl

inline std::string shortType(const std::string& input) noexcept {
    if (input.empty()) return input;
    std::string cleanType = Meta::Implementation::type(input);
    cleanType = String::replaceAll(cleanType, Meta::type<Vertex>(), "Vertex");
    cleanType = String::replaceAll(cleanType, Meta::type<Edge>(), "Edge");
    cleanType = String::replaceAll(cleanType, Meta::type<StopId>(), "StopId");
    cleanType = String::replaceAll(cleanType, Meta::type<RouteId>(), "RouteId");
    cleanType = String::replaceAll(cleanType, Meta::type<TripId>(), "TripId");
    cleanType = String::replaceAll(cleanType, Meta::type<StopIndex>(), "StopIndex");
    cleanType = String::replaceAll(cleanType, Meta::type<ConnectionId>(), "ConnectionId");
    cleanType = String::replaceAll(cleanType, Meta::type<HyperVertexId>(), "HyperVertexId");
    cleanType = String::replaceAll(cleanType, "basic_string<char, std::char_traits<char>, std::allocator<char> >", "string");
    cleanType = String::replaceAll(cleanType, " >", ">");
    cleanType = String::replaceAll(cleanType, "Meta::", "");
    cleanType = String::replaceAll(cleanType, "std::", "");
    cleanType = String::replaceAll(cleanType, "__gnu_cxx::", "");
    cleanType = String::replaceAll(cleanType, "__cxx11::", "");
    cleanType = String::replaceAll(cleanType, "__debug::", "");
    cleanType = String::replaceAll(cleanType, "__ops::", "");
    cleanType = String::replaceAll(cleanType, "__normal_", "");
    cleanType = String::replaceAll(cleanType, "Implementation", "");
    while (cleanType.size() > 130) {
        size_t beginIndex = cleanType.find_last_of('<');
        if (beginIndex >= cleanType.size()) break;
        size_t endIndex = cleanType.find_first_of('>', beginIndex);
        if (endIndex >= cleanType.size()) break;
        if (endIndex - beginIndex + 1 < 6) {
            cleanType[beginIndex] = '{';
            cleanType[endIndex] = '}';
        } else {
            cleanType.replace(beginIndex, endIndex - beginIndex + 1, ";");
        }
    }
    cleanType = String::replaceAll(cleanType, ";", "<...>");
    cleanType = String::replaceAll(cleanType, "}", ">");
    cleanType = String::replaceAll(cleanType, "{", "<");
    return cleanType;
}

inline bool printStackTrace(std::ostream& out = std::cout, size_t maxFrames = 63) {
    out << "Stack trace: " << std::endl;
    void** addressList = new void*[maxFrames + 1];

    size_t numFrames = backtrace(addressList, maxFrames);
    if (numFrames <= 0) {
        out << "    ERROR! No stack trace available!" << std::endl;
        return false;
    }

    char** symbolList = backtrace_symbols(addressList, numFrames);

    size_t functionNameSize = 256;
    char* functionName = (char*)malloc(functionNameSize);

    for (size_t i = 0; i < numFrames; i++) {
        char* name = 0;
        char* caller = 0;
        for (char* j = symbolList[i]; *j; j++) {
            if (*j == '(') {
                *j = '\0';
                name = j + 1;
            } else if (*j == '+') {
                *j = '\0';
                caller = j + 1;
            } else if (*j == ')' && caller) {
                *j = '\0';
                break;
            }
        }
        if (name && caller && name < caller) {
            int status;
            char* demangledFunctionName = abi::__cxa_demangle(name, functionName, &functionNameSize, &status);
            if (status == 0) {
                functionName = demangledFunctionName;
                out << "   " << shortType(functionName) << std::endl;
            }
            else {
                out << "   " << name << std::endl;
            }
        }

    }

    free(functionName);
    free(symbolList);
    delete[] addressList;
    return false;
}

std::vector<int> debugStack;

void segfaultSigaction(int, siginfo_t* si, void*) {
    printf("Caught segfault at address %p\n", si->si_addr);
    std::cout << "Debug Stack:" << std::endl;
    for (const int line : debugStack) {
        std::cout << "   " << line << std::endl;
    }
    exit(0);
}

#define __Push__ debugStack.push_back(__LINE__)
#define __Pop__ debugStack.pop_back()
#define __Change__ debugStack.back() = __LINE__
#define __Debug_Stack__ struct sigaction __sa__; memset(&__sa__, 0, sizeof(struct sigaction)); sigemptyset(&__sa__.sa_mask); __sa__.sa_sigaction = segfaultSigaction; __sa__.sa_flags = SA_SIGINFO; sigaction(SIGSEGV, &__sa__, NULL)
