#include "Protocols/gtfs-realtime-OVapi.pb.h"

#include "../Helpers/String/String.h"
#include "../Helpers/String/Enumeration.h"
#include "../Helpers/String/TextFileUtils.h"
#include "../Helpers/Vector/Vector.h"
#include "../Helpers/Vector/Permutation.h"
#include "../Shell/Shell.h"

#include "../Helpers/ConfigFile.h"
#include "../Helpers/Console/CommandLineParser.h"

#include "../Helpers/Assert.h"
#include "../Helpers/Debug.h"
#include "../Helpers/Timer.h"
#include "../Helpers/Calendar.h"
#include "../Helpers/MultiThreading.h"
#include "../Helpers/IO/ParserCSV.h"
#include "../Helpers/FileSystem/Curl.h"

#include "../DataStructures/Container/Set.h"
#include "../DataStructures/GTFS/Data.h"
#include "../DataStructures/CSA/Data.h"
#include "../DataStructures/Intermediate/Data.h"
#include "../DataStructures/Intermediate/Data.h"

#include <iostream>
#include <fstream>
#include <random>
#include <thread>
#include <vector>
#include <string>
#include <chrono>
#include <set>

#include <google/protobuf/text_format.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

using namespace Shell;
Shell::Shell shell;

void usage(char *program) {
    std::cout << "Usage: " << program << " <options>" << std::endl
         << "   -gtfs <url>      -- URL of the static GTFS data" << std::endl
         << "   -pb <url>        -- URL of the real time GTFS feed" << std::endl
         << "   -data <folder>   -- Working directory of the program" << std::endl
         << "   -interactive     -- Start interactive shell" << std::endl
         << "   -v               -- Verbose mode" << std::endl
         << "   -h               -- Print help" << std::endl
         << std::endl;
    exit(0);
}

class Test : public Command {

public:
    Test(BasicShell& shell) : shell(shell) {
        shell.addCommand(this);
    }

    virtual std::string name() const noexcept {
        return "test";
    }

    virtual std::string helpText() const noexcept {
        return name() + " <input protocol buffer file> <output text file>: Converts a protocol buffer to a human readable text.";
    }

    virtual void execute(const std::string& parameter) noexcept {
        std::vector<std::string> tokens = String::split(String::trim(parameter), ' ');
        if (tokens.size() < 1) tokens.emplace_back(shell.ask("Input protocol buffer file A> "));
        tokens[0] = "/home/tobias/Data/examples/csa";
        CSA::Data csa = CSA::Data::FromBinary(tokens[0]);
        csa.printInfo();
        Intermediate::Data inter1 = Intermediate::Data::FromCSA(csa);
        inter1.writeCSV("/home/tobias/Data/examples/inter1");
        Intermediate::Data inter2 = Intermediate::Data::FromCSV("/home/tobias/Data/examples/inter1");
        inter2.writeCSV("/home/tobias/Data/examples/inter2");

        // if (tokens[0] == "") return;
        // std::vector<std::string> entitiesA = getEntities(tokens[0]);
        // std::cout << "Entities in A: " << String::prettyInt(entitiesA.size()) << std::endl;
        // if (tokens.size() < 2) tokens.emplace_back(shell.ask("Input protocol buffer file A> "));
        // if (tokens[1] == "") return;
        // std::vector<std::string> entitiesB = getEntities(tokens[1]);
        // std::cout << "Entities in B: " << String::prettyInt(entitiesB.size()) << std::endl;
        // std::cout << "Diff: " << String::prettyInt(entitiesB.size() - entitiesA.size()) << std::endl;
        // Set<std::string> entitySet;
        // for (const std::string& entity : entitiesA) {
        //     entitySet.insert(entity);
        // }
        // std::cout << "Entities in A: " << String::prettyInt(entitySet.size()) << std::endl;
        // std::vector<std::string> entitiesDiff;
        // for (const std::string& entity : entitiesB) {
        //     if (!entitySet.contains(entity)) {
        //         entitiesDiff.emplace_back(entity);
        //     }
        // }
        // std::cout << "Diff: " << String::prettyInt(entitiesDiff.size()) << std::endl;
        // if (tokens.size() < 3) tokens.emplace_back(shell.ask("Output text file> "));
        // if (tokens[2] == "") return;
        // std::ofstream out(tokens[2]);
        // for (const std::string& entity : entitiesDiff) {
        //     out << entity << std::endl;
        // }
    }

private:
    std::vector<std::string> getEntities(const std::string& fileName) const noexcept {
        std::vector<std::string> result;
        std::ifstream in(fileName, std::ios::binary);
        transit_realtime::FeedMessage feed_message;
        if(!feed_message.ParseFromIstream(&in)) {
            error("Error while parsing protocol buffer!");
        }
        std::string text;
        google::protobuf::TextFormat::PrintToString(feed_message, &text);
        std::stringstream ss;
        int braces = 0;
        bool copy = false;
        for (size_t i = 0; i < text.size(); i++) {
            if (braces == 0) {
                if (String::containsSubString(text, i, "entity {")) copy = true;
            }
            if (copy && !String::isWhiteSpace(text[i])) {
                ss << text[i];
            }
            if (text[i] == '{') {
                braces++;
            } else if (text[i] == '}') {
                braces--;
                if (braces <= 0 && copy) {
                    result.emplace_back(ss.str());
                    ss.str("");
                    ss.clear();
                    copy = false;
                }
            }
        }
        return result;
    }

    std::vector<std::string> getEntityIDs(const std::string& fileName) const noexcept {
        std::vector<std::string> result;
        std::ifstream in(fileName, std::ios::binary);
        transit_realtime::FeedMessage feed_message;
        if(!feed_message.ParseFromIstream(&in)) {
            error("Error while parsing protocol buffer!");
        }
        std::string text;
        google::protobuf::TextFormat::PrintToString(feed_message, &text);
        std::stringstream ss;
        int braces = 0;
        bool copy = false;
        bool id = false;
        bool quote = false;
        for (size_t i = 0; i < text.size(); i++) {
            if (braces == 0) {
                if (String::containsSubString(text, i, "entity {")) copy = true;
            }
            if (copy &&  (braces == 1)) {
                if (String::containsSubString(text, i, "id: ")) id = true;
            }
            if (id && quote && (text[i] != '"')) {
                ss << text[i];
            }
            if (text[i] == '{') {
                braces++;
                id = false;
            } else if (text[i] == '}') {
                braces--;
                id = false;
                if (braces <= 0 && copy) {
                    result.emplace_back(ss.str());
                    ss.str("");
                    ss.clear();
                    copy = false;
                }
            } else if (text[i] == '"') {
                quote = !quote;
            }
        }
        return result;
    }

private:
    BasicShell& shell;
};

class RealTimeGtfsToText : public Command {

public:
    RealTimeGtfsToText(BasicShell& shell) : shell(shell) {
        shell.addCommand(this);
    }

    virtual std::string name() const noexcept {
        return "realTimeGtfsToText";
    }

    virtual std::string helpText() const noexcept {
        return name() + " <input protocol buffer file> <output text file>: Converts a protocol buffer to a human readable text.";
    }

    virtual void execute(const std::string& parameter) noexcept {
        std::vector<std::string> tokens = String::split(String::trim(parameter), ' ');
        if (tokens.size() < 1) tokens.emplace_back(shell.ask("Input protocol buffer file> "));
        if (tokens[0] == "") return;
        std::ifstream in(tokens[0], std::ios::binary);
        transit_realtime::FeedMessage feed_message;
        if(!feed_message.ParseFromIstream(&in)) {
            error("Error while parsing protocol buffer!");
        }
        if (tokens.size() < 2) tokens.emplace_back(shell.ask("Output text file> "));
        if (tokens[1] == "") return;
        std::ofstream out(tokens[1]);
        google::protobuf::io::OstreamOutputStream gout(&out);
        google::protobuf::TextFormat::Print(feed_message, &gout);
    }

private:
    BasicShell& shell;
};

class Download : public Command {

public:
    Download(BasicShell& shell) : shell(shell) {
        shell.addCommand(this);
    }

    virtual std::string name() const noexcept {
        return "download";
    }

    virtual std::string helpText() const noexcept {
        return name() + " <url> <output file> <Unzip directory>: Downloads a url to a file and unzips it.";
    }

    virtual void execute(const std::string& parameter) noexcept {
        std::vector<std::string> tokens = String::split(String::trim(parameter), ' ');
        if (tokens.size() < 1) tokens.emplace_back(shell.ask("Url> "));
        if (tokens[0] == "") return;
        if (tokens.size() < 2) tokens.emplace_back(shell.ask("Output file> "));
        if (tokens[1] == "") return;
        Curl::downloadFile(tokens[0], tokens[1]);
        if (tokens.size() < 3) tokens.emplace_back(shell.ask("Unzip directory> "));
        if (tokens[2] == "") return;
        FileSystem::unzip(tokens[1], tokens[2]);
    }

private:
    BasicShell& shell;
};

inline void terminate(bool& running) noexcept {
    running = true;
    shell.blue("Real time GTFS parser is running. Press any key to stop.").endl();
    getchar();
    running = false;
    shell.blue("Terminating real time GTFS parser...").endl();
}

inline void comment(const bool verbose, const std::string& text) noexcept {
    if (verbose) shell.yellow(String::dateString() + " " + String::timeString() + ": " + text).endl();
}

int main(int argc, char **argv) {
    checkAsserts();

    CommandLineParser clp(argc, argv);
    if (clp.isSet("h")) usage(argv[0]);

    const bool verbose = clp.isSet("v");

    new Test(shell);
    new RealTimeGtfsToText(shell);
    new Download(shell);
    if (clp.isSet("interactive")) {
        shell.run();
    } else {
        const std::string gtfsURL(clp.value<std::string>("gtfs", ""));
        const std::string pbURL(clp.value<std::string>("pb", ""));
        const std::string dataFolder(clp.value<std::string>("data", ""));
        const std::string gtfsDownload(FileSystem::extendPath(dataFolder, "gtfs.zip"));
        const std::string pbDownload(FileSystem::extendPath(dataFolder, "tripUpdates.pb"));
        const std::string gtfsDataDirectory(FileSystem::extendPath(dataFolder, "gtfs/"));
        const std::string gtfsDataBaseFile(FileSystem::extendPath(gtfsDataDirectory, "data.binary"));
        const std::string intermediateDataDirectory(FileSystem::extendPath(dataFolder, "intermediate/"));
        const std::string intermediateDataBaseFile(FileSystem::extendPath(intermediateDataDirectory, "data"));
        const std::string intermediateDataBackupBaseFile(FileSystem::extendPath(intermediateDataDirectory, "backup"));
        std::string realTimeDataFolder(FileSystem::extendPath(intermediateDataDirectory, "00000000/"));

        FileSystem::makeDirectory(gtfsDataDirectory);
        FileSystem::makeDirectory(intermediateDataDirectory);

        bool running = true;
        std::thread keyListener(terminate, std::ref(running));
        while (running) {
            std::string newRealTimeDataFolder(FileSystem::extendPath(intermediateDataDirectory, String::dateString() + "/"));
            if (realTimeDataFolder != newRealTimeDataFolder) { // New day or program start
                realTimeDataFolder = newRealTimeDataFolder;
                FileSystem::makeDirectory(realTimeDataFolder);
                comment(verbose, "Downloading GTFS data from " + gtfsURL + " to " + gtfsDownload + ".");
                Curl::downloadFile(gtfsURL, gtfsDownload);
                comment(verbose, "Unzipping GTFS data from " + gtfsDownload + " to " + gtfsDataDirectory + ".");
                FileSystem::unzip(gtfsDownload, gtfsDataDirectory, verbose);
                comment(verbose, "Loading GTFS data from " + gtfsDownload + " to " + gtfsDataDirectory + ".");
                GTFS::Data gtfsData(GTFS::Data::FromGTFS(gtfsDataDirectory, verbose));
                gtfsData.write(gtfsDataBaseFile);

            }

            std::this_thread::sleep_for(std::chrono::milliseconds(500));
        }
        keyListener.join();
    }

    return 0;
}