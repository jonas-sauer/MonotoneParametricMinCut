#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <map>

#include "String/String.h"

class ConfigFile {

public:
    friend std::ostream& operator<<(std::ostream& out, const ConfigFile& cf);

public:
    ConfigFile(const std::string& filename, const bool verbose = true) :
        modified(false) {
        setFilename(filename);
        read(verbose);
    }

    template<typename T>
    inline T get(const std::string& key, const T& defaultValue) {
        if (!contains(key)) {
            data[key] = String::lexicalCast<std::string>(defaultValue);
            modified = true;
        }
        return get<T>(key);
    }

    inline bool get(const std::string& key, const bool defaultValue) {
        return String::lexicalCast<bool>(get<std::string>(key, defaultValue ? "true" : "false"));
    }

    template<typename T>
    inline T get(const std::string& key) const {
        if (!contains(key)) return T();
        return String::lexicalCast<T>(data.find(key)->second);
    }

    inline std::string get(const std::string& key) const {
        if (!contains(key)) return "";
        return data.find(key)->second;
    }

    template<typename T>
    inline void set(const std::string& key, const T& value) {
        if (!contains(key) || data[key] != String::lexicalCast<std::string>(value)) {
            data[key] = String::lexicalCast<std::string>(value);
            modified = true;
        }
    }

    inline void set(const std::string& key, const bool value) {
        set<std::string>(key, value ? "true" : "false");
    }

    template<typename T>
    inline T operator[](const std::string& key) const {
        if (!contains(key)) return T();
        return String::lexicalCast<T>(data.find(key)->second);
    }

    inline std::string operator[](const std::string& key) const {
        if (!contains(key)) return "";
        return data.find(key)->second;
    }

    inline std::string& operator[](const std::string& key) {
        modified = true;
        return data[key];
    }

    inline bool contains(const std::string& key) const {
        return data.count(key) == 1;
    }

    inline bool isSet(const std::string& key) const {
        if (!contains(key)) return false;
        return String::lexicalCast<bool>(data.find(key)->second);
    }

    inline void read(const bool verbose = true) {
        clear();
        std::ifstream inputFile(getFilename());
        if (!inputFile.is_open()) return;
        while (!inputFile.eof()) {
            std::string line;
            getline(inputFile, line);
            line = String::trim(line);
            int splitPos = line.find(": ");
            if (splitPos < 0) continue;
            set(line.substr(0,splitPos), line.substr(splitPos + 2));
        }
        inputFile.close();
        modified = false;

        if (verbose) std::cout << std::endl << "Reading configuration file from " << getFilename() << ":" << std::endl << *this << std::endl;
    }

    inline void write(const bool verbose = true) const {
        std::ofstream outputFile(getFilename());
        outputFile << *this;
        outputFile.close();
        modified = false;

        if (verbose) std::cout << std::endl << "Writing configuration file to " << getFilename() << ":" << std::endl << *this << std::endl;
    }

    inline void writeIfModified(const bool verbose = true) const {
        if (modified) {
            write(verbose);
        }
    }

    inline void clear() {
        data.clear();
        modified = false;
    }

    inline const std::string& getFilename() const {
        return fname;
    }

    inline void setFilename(const std::string& filename) {
        if (!String::endsWith(filename, ".conf")) {
            setFilename(filename + ".conf");
        } else if (fname != filename) {
            fname = filename;
            modified = true;
        }
    }

private:
    inline int maxKeySize() const {
        int maxSize = 0;
        for (auto const& entry : data) {
            int size = entry.first.size();
            if (maxSize < size) maxSize = size;
        }
        return maxSize;
    }

    std::map<std::string, std::string> data;
    std::string fname;
    mutable bool modified;

};

inline std::ostream& operator<<(std::ostream& out, const ConfigFile& cf) {
    int width = cf.maxKeySize() + 2;
    for (auto const& entry : cf.data) {
        out << std::setw(width) << std::left << (entry.first + ": ") << entry.second << std::endl;
    }
    return out;
}
