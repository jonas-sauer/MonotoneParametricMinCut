#pragma once

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include <termios.h>
#include <unistd.h>
#include <stdio.h>
#include <limits.h>

#include <dirent.h>

#include <sys/stat.h>
#include <sys/wait.h>
#include <sys/types.h>

#include "../Assert.h"
#include "../String/String.h"

namespace FileSystem {
    inline bool isDirectory(const std::string& path) noexcept {
        DIR* dir;
        if ((dir = opendir(path.c_str())) != NULL) {
            closedir(dir);
            return true;
        } else {
            return false;
        }
    }

    inline std::string getParentDirectory(const std::string& fileName) noexcept {
        if (fileName.size() < 2) return "";
        size_t directoryEnd = fileName.find_last_of('/', fileName.size() - 2);
        if (directoryEnd >= fileName.size()) return "";
        return fileName.substr(0, directoryEnd);
    }

    inline std::string ensureExtension(const std::string& fileName, const std::string& extension) noexcept {
        if (String::endsWith(fileName, extension)) {
            return fileName;
        } else {
            return fileName + extension;
        }
    }
    inline void makeDirectory(const std::string& path) noexcept {
        if (path.empty()) return;
        if (isDirectory(path)) return;
        makeDirectory(getParentDirectory(path));
        mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }

    inline const std::string& ensureDirectoryExists(const std::string& fileName) noexcept {
        const std::string parentDirectory = getParentDirectory(fileName);
        if (parentDirectory == "") return fileName;
        if (isDirectory(parentDirectory)) return fileName;
        makeDirectory(parentDirectory);
        return fileName;
    }

}
