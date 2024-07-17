#include "../Commands/Assignment.h"
using UsedCommand = CapacityAssignment;

#include <string>
#include <sstream>

#include "../../Helpers/String/String.h"
#include "../../Helpers/String/Enumeration.h"
#include "../../Helpers/MultiThreading.h"

int main(int argc, char** argv) {
    int index = 1;
    int core = 1;
    if ((argc >= 2) && ((std::string(argv[1]) == "-h") || (std::string(argv[1]) == "-help"))) {
        Shell::Shell shell;
        UsedCommand command(shell);
        std::cout << command.helpText() << std::endl;
        exit(0);
    }
    if ((argc >= 3) && (std::string(argv[1]) == "-core")) {
        core = String::lexicalCast<int>(argv[2]);
        index = 3;
    }
    Enumeration parameters;
    while (index < argc) {
        parameters << argv[index] << Sep(" ");
        index++;
    }

    pinThreadToCoreId(core);
    checkAsserts();

    Shell::Shell shell;
    shell.setReportParameters(true);
    UsedCommand* command = new UsedCommand(shell);
    Timer timer;
    static_cast<Command*>(command)->execute(parameters.str());
    shell << grey("[Finished in ", String::msToString(timer.elapsedMilliseconds()), "]") << newLine;
}
