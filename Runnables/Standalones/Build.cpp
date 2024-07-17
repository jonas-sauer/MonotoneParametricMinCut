#include <string>
#include <cstdlib>

#include "../../Shell/Shell.h"
#include "../../Shell/ParameterizedCommand.h"

#include "../../DataStructures/Container/Map.h"

#include "../../Helpers/String/String.h"
#include "../../Helpers/String/TextFileUtils.h"

class Build : public Shell::ParameterizedCommand {

public:
    const static inline std::string path = "../Commands/";

public:
    Build(Shell::BasicShell& shell) :
        ParameterizedCommand(shell, "build", "Displays all files in the current working directory.") {
        addParameter("Command");
        DIR* dir;
        struct dirent* ent;
        if ((dir = opendir(path.c_str())) != nullptr) {
            std::vector<std::string> dirs;
            while ((ent = readdir(dir)) != nullptr) {
                dirs.push_back(ent->d_name);
            }
            closedir(dir);
            std::sort(dirs.begin(), dirs.end());
            for (const std::string& name : dirs) {
                if (!String::endsWith(name, ".h")) continue;
                std::vector<std::string> lines = String::split(TextFile::read(path + name), '\n');
                for (const std::string& line : lines) {
                    std::vector<std::string> tokens = String::split(String::trim(line), ' ');
                    for (size_t i = 0; i + 1 < tokens.size(); i++) {
                        if ((tokens[i] != "class") && (tokens[i] != "struct")) continue;
                        commands.insert(tokens[i + 1], name);
                        names.emplace_back(tokens[i + 1]);
                    }
                }
            }
        } else {
            std::cout << "Could not open directory: \"" << path << "\"." << std::endl;
        }
    }

    virtual void execute() noexcept {
        const std::string name = getParameter("Command");
        if (!commands.contains(name)) {
            shell.error("Unknown command name: ", name, "!");
            return;
        }
        std::string makefile = TextFile::read("Makefile");
        replace(makefile, "NAME=", "\nCC", name);
        TextFile::write("Makefile", makefile);
        std::string runCommand = TextFile::read("RunCommand.cpp");
        replace(runCommand, "#include \"", "\"\nusing", path + commands[name]);
        replace(runCommand, "UsedCommand = ", ";\n\n", name);
        TextFile::write("RunCommand.cpp", runCommand);
        int result = std::system("make -B RunCommand");
        if (result != 0) {
            shell.error("Failure code: ", result, "!");
        }
    }

    virtual std::vector<std::string> parameterSuggestions() const {
        return names;
    }

private:
    inline void replace(std::string& s, const std::string& beginString, const std::string& endString, const std::string& replacement) const noexcept {
        const size_t begin = s.find(beginString) + beginString.size();
        const size_t end = s.find(endString, begin);
        s.replace(begin, end - begin, replacement);
    }

private:
    Map<std::string, std::string> commands;
    std::vector<std::string> names;

};

int main(int argc, char** argv) {
    Shell::Shell shell;
    Build* command = new Build(shell);
    static_cast<Shell::Command*>(command)->execute((argc > 1) ? argv[1] : "");
}
