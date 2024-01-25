#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <sys/wait.h>

#include "String.h"

#include "../Assert.h"
#include "../Debug.h"

#include "../IO/File.h"
#include "../FileSystem/FileSystem.h"

namespace TextFile {

inline std::string read(std::ifstream& file) {
    std::stringstream ss;
    while (!file.eof()) {
        std::string line;
        getline(file, line);
        ss << line;
        if (!file.eof()) ss << '\n';
    }
    return ss.str();
}

inline std::string read(const std::vector<std::string>& fileNameAliases) {
    return read(IO::IFStream(fileNameAliases));
}

inline std::string read(const std::string& fileName) {
    return read(IO::IFStream(fileName));
}

inline void write(const std::string& fileName, const std::string& text) {
    std::ofstream file(fileName);
    AssertMsg(file.is_open(), "Could not open file: " << fileName);
    file << text;
    file.flush();
    file.close();
}

inline void trim(const std::string& fromFile, std::ostream& to, const size_t minWordLength = 1) {
    std::ifstream from(fromFile);
    Assert(from.is_open());
    while (!from.eof()) {
        std::string line;
        getline(from, line);
        line = String::trim(line);
        line = String::toLower(line);
        std::vector<std::string> words = String::split(line, ' ');
        for (std::string& word : words) {
            if (word.size() < minWordLength) continue;
            for (char c : word) {
                if (c >= 'a' && c <= 'z') {
                    to << c;
                }
            }
        }
    }
}

inline void trim(const std::string& fromFile, const std::string& toFile, const size_t minWordLength = 1) {
    std::ofstream to(toFile);
    Assert(to);
    Assert(to.is_open());
    trim(fromFile, to, minWordLength);
}

inline std::string trim(const std::string& fromFile, const size_t minWordLength = 1) {
    std::stringstream to;
    trim(fromFile, to, minWordLength);
    return to.str();
}

inline void pdfToText(std::string fromFile, std::string toFile) {
    pid_t pid = fork();
    switch (pid) {
    case -1:
        error("Failure during fork()!");
        exit(1);
    case 0:
        fromFile = FileSystem::extendPath(FileSystem::getWorkingDirectory(), fromFile);
        toFile = FileSystem::extendPath(FileSystem::getWorkingDirectory(), toFile);
        execl("/usr/bin/pdftotext", "pdftotext", fromFile.c_str(), toFile.c_str(), nullptr);
        error("/usr/bin/pdftotext failed!");
        exit(1);
    default:
        int status;
        while (!WIFEXITED(status)) {
            waitpid(pid, &status, 0);
        }
    }
}

}
