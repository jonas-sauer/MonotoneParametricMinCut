
#include <iostream>
#include <string>
#include <vector>
#include <thread>

#include "../Helpers/Meta.h"
#include "../Helpers/Timer.h"
#include "../Helpers/Types.h"
#include "../Helpers/Debug.h"
#include "../Helpers/PStack.h"
#include "../Helpers/Helpers.h"
#include "../Helpers/HighlightText.h"
#include "../Helpers/String/String.h"
#include "../Helpers/Console/CommandLineParser.h"

#include "../Shell/BasicShell.h"
#include "../Shell/LineBuffer.h"

#include "../DataStructures/Container/Map.h"
#include "../DataStructures/Graph/Classes/GraphInterface.h"

void usage(const char *program) {
    cout << "Usage: " << program << " <options>" << endl
         << endl
         << "-p   <pid>       -- PID of the program to analyze" << endl
         << "-n   <name>      -- Name of the program to analyze" << endl
         << "-si  <integer>   -- Sampling Interval in ms (default: 10)" << endl
         << "-pi  <integer>   -- Print Interval in ms (default: 1000)" << endl
         << "-h               -- print help" << endl
         << endl;
    exit(0);
}

size_t firstRow = 0;
size_t lastRow = 0;
size_t maxTreeDepth = 5;
bool running = true;
bool printTree = true;
size_t currentRow = 0;

static char cmd[128];

static char *cmdLine(int pid) {
    int fd, len = -1, i;

    sprintf(cmd, "/proc/%d/cmdline", pid);
    if ((fd = open(cmd, O_RDONLY)) >= 0 && (len = read(fd, cmd, sizeof(cmd))) > 0) {
        for (i = 0; i < len; i++) if (!cmd[i]) cmd[i] = ' ';
        for ( ; len > 0 && cmd[len - 1] <= ' '; len--);
        cmd[len] = 0;
        if ((unsigned int) len >= sizeof(cmd) - 4) {
            strcpy(&cmd[sizeof(cmd) - 4], "...");
        }
    } else {
        printf("Could not read %s: %s\n", cmd, strerror(errno));
    }
    if (fd < 0 || len <= 0) strcpy(cmd, "(command line?)");
    if (fd >= 0) close(fd);

    return cmd;
}

struct StackFrame {
    StackFrame() :
        count(0) {
    }

    void print(Map<std::string, size_t>& callCounts, const double callCount, const int depth, const int maxDepth, const std::string& prefix = " ") {
        size_t numberOfChildren = children.size();
        for (auto& [name, stackFrame] : children) {
            if (currentRow >= lastRow) return;
            if (numberOfChildren == 1) {
                if (currentRow >= firstRow) std::cout << "\x1b[2K" << prefix << "└─";
                stackFrame.printData(name, callCounts[name], callCount);
                if (depth < maxDepth) {
                    stackFrame.print(callCounts, callCount, depth + 1, maxDepth, prefix + "  ");
                }
            } else {
                if (currentRow >= firstRow) std::cout << "\x1b[2K" << prefix << "├─";
                stackFrame.printData(name, callCounts[name], callCount);
                if (depth < maxDepth) {
                    stackFrame.print(callCounts, callCount, depth + 1, maxDepth, prefix + "│ ");
                }
            }
            numberOfChildren--;
        }
    }

    void printData(const std::string& name, const int totalCount, const double callCount) {
        if (currentRow >= firstRow) {
            if (children.empty()) {
                std::cout << "──";
            } else {
                std::cout << "┬─";
            }
            std::cout << "[" << String::percent(count / callCount) << "] " << shortType(name);
            if (count != totalCount) std::cout << " " << grey("[", String::percent(totalCount / callCount), "]");
            std::cout << "\n";
        }
        currentRow++;
    }

    size_t count;
    Map<std::string, StackFrame> children;
};

void input() {
    while (running) {
        uint64_t c = Shell::getch();
        if (c == 'q' || c == 'Q') { // Q
            running = false;
        } else if (c == 27) { // Escape sequence
            uint64_t c1 = Shell::getch();
            uint64_t c2 = Shell::getch();
            c = (c1 * 128) + c2;
            if (c2 >= 50 && c2 <= 55) {
                uint64_t c3 = Shell::getch();
                c = (c * 128) + c3;
            }
            if (c == 11713) { // Up Key
                if (firstRow > 0) firstRow--;
                printTree = true;
            } else if (c == 11714) { // Down Key
                firstRow++;
                printTree = true;
            } else if (c == 11715) { // Right Key
                maxTreeDepth++;
                printTree = true;
            } else if (c == 11716) { // Left Key
                if (maxTreeDepth > 1) maxTreeDepth--;
                printTree = true;
            }
        } else if (c == 127) { // Backspace
            firstRow = 0;
        }
    }
}

int main(int argc, char **argv) {
    CommandLineParser clp(argc, argv);
    if (clp.isSet("h")) usage(argv[0]);

    int pid(clp.value<int>("p", -1));
    if (pid < 0) pid = getProcIdByName(clp.value<std::string>("n", ""));
    if (pid < 0) usage(argv[0]);
    const int si(clp.value<int>("si", 50));
    const int pi(clp.value<int>("pi", 1000));

    thePid = pid;
    if (attach(thePid) != 0) {
        fprintf(stderr, "Could not attach to target %d: %s.\n", thePid, strerror(errno));
    } else {
        loadSymbols(thePid, false);
        detachall();

        for (size_t i = 2; i <= Shell::getScreenHeight(); i++) {
            std::cout << "\n";
        }

        Map<std::string, size_t> callCounts;
        StackFrame rootFrame;

        Timer printTimer;

        std::thread inoutThread(input);

        while (running) {
            Timer sampleTimer;
            attach(thePid);

            std::vector<std::string> stack = getStack(thePid);

            StackFrame* currentFrame = &rootFrame;
            while (!stack.empty()) {
                if (stack.back() == "_fini") {
                    stack.pop_back();
                    continue;
                }
                if (!currentFrame->children.contains(stack.back())) {
                    currentFrame->children.insert(stack.back(), StackFrame());
                }
                currentFrame = &((currentFrame->children)[stack.back()]);
                currentFrame->count++;
                if (!callCounts.contains(stack.back())) {
                    callCounts.insert(stack.back(), 0);
                }
                callCounts[stack.back()]++;
                stack.pop_back();
            }
            rootFrame.count++;

            if (printTimer.elapsedMilliseconds() >= pi || printTree) {
                std::cout << '\r';
                for (size_t i = 0; i <= Shell::getScreenHeight(); i++) {
                    std::cout << "\x1b[1A";
                }
                currentRow = 0;
                lastRow = firstRow + Shell::getScreenHeight() - 2;
                rootFrame.print(callCounts, rootFrame.count, 0, maxTreeDepth);
                while (currentRow < lastRow) {
                    std::cout << "\x1b[2K\n";
                    currentRow++;
                }
                std::cout << "\x1b[2K" << std::endl;
                printTimer.restart();
                printTree = false;
            }

            detachall();

            sleep(std::max(1, int(si - sampleTimer.elapsedMilliseconds())));
        }

        inoutThread.join();

        std::cout << '\r';
        for (size_t i = 0; i <= Shell::getScreenHeight(); i++) {
            std::cout << "\x1b[1A";
        }
        firstRow = 0;
        currentRow = 0;
        lastRow = 500;
        rootFrame.print(callCounts, rootFrame.count, 0, maxTreeDepth);

    }

    resetData();

    exit(0);
}
