#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
using std::cout;
using std::flush;
using std::endl;
using std::vector;
using std::string;

#include "../Helpers/Types.h"
#include "../Helpers/Vector/Vector.h"
#include "../Helpers/String/String.h"
#include "../Helpers/String/Enumeration.h"
#include "../Helpers/String/TextFileUtils.h"
#include "../Helpers/FileSystem/FileSystem.h"
#include "../Helpers/Console/CommandLineParser.h"

void usage(char *program) {
    cout << "Usage: " << program << " <Pdf_A> <Pdf_B> <options>" << endl
         << "   -w <integer>  -- minimal Word length (default: 3)" << endl
         << "   -m <integer>  -- minimal Match length (default: 50)" << endl
         << "   -a <double>   -- minimal Accordance between matched strings (default: 1.0)" << endl
         << "   -pa           -- print text a in complete" << endl
         << "   -pb           -- print text b in complete" << endl
         << "   -h            -- print help" << endl
         << endl
         << "The program 'pdftotext' is required for " << program << "!" << endl
         << "   -ptt <file>   -- location of the pdftotext executable (default: /usr/bin/pdftotext)" << endl
         << endl;
    exit(0);
}

struct TextPos {
    inline bool operator<(const TextPos& o) const {return (line < o.line) || ((line == o.line) && (index < o.index));}
    size_t line;
    size_t index;
};

struct Text {
    Text(const string& fromFile, int minWordLength = 1) :
        filename(fromFile),
        trimText("") {
        addFile(fromFile, minWordLength);
    }
    Text(const vector<string>& fromFiles, int minWordLength = 1) {
        Enumeration filenames;
        for (const string& fromFile : fromFiles) {
            filenames << fromFile << sep;
            addFile(fromFile, minWordLength);
        }
        filename = filenames.str();
    }
    inline void addFile(const string& fromFile, int minWordLength = 1) {
        std::ifstream from(fromFile);
        Assert(from.is_open());
        std::stringstream to;
        to << trimText;
        while (!from.eof()) {
            std::string line;
            getline(from, line);
            lines.push_back(line);
            line = String::toLower(line);
            int wordLength = 0;
            for (size_t i = 0; i < line.size(); i++) {
                char c = line[i];
                if (c >= 'a' && c <= 'z') {
                    wordLength++;
                } else {
                    if (wordLength >= minWordLength) {
                        for (size_t j = i - wordLength; j < i; j++) {
                            to << line[j];
                            trimTextToPos.push_back(TextPos{lines.size() - 1, j});
                        }
                    }
                    wordLength = 0;
                }
            }
            if (wordLength >= minWordLength) {
                for (size_t j = line.size() - wordLength; j < line.size(); j++) {
                    to << line[j];
                    trimTextToPos.push_back(TextPos{lines.size() - 1, j});
                }
            }
        }
        trimText = to.str();
    }
    inline void print() {
        for (string& line : lines) cout << line << endl;
        cout << endl << trimText << endl;
    }
    string filename;
    inline int size() const {return trimText.size();}
    vector<string> lines;
    string trimText;
    vector<TextPos> trimTextToPos;
};

struct TextSection {
    TextSection(const Text& text, const int begin, const int end) {
        Assert(begin >= 0);
        Assert(end < text.size());
        TextPosBegin = text.trimTextToPos[begin];
        TextPosEnd = text.trimTextToPos[end];
    }
    TextPos TextPosBegin;
    TextPos TextPosEnd;
};

struct WeakMatch {
    WeakMatch(const int i, const int j, const int size) :
        i(i),
        j(j),
        size(size),
        sizeI(size),
        sizeJ(size),
        accordance(size) {
    }
    inline bool operator<(const WeakMatch& m) const {return (m.i - m.sizeI >= i) && (m.j - m.sizeJ >= j);}
    inline int dist(const WeakMatch& m) const {return std::max(m.i - i, m.j - j) - m.size;}
    inline double accordancePercentage() const {return accordance / ((double)size);}
    int i;
    int j;
    int size;
    int sizeI;
    int sizeJ;
    int accordance;
};

struct Match {
    Match(const Text& a, const Text& b, const int aIndex, const int bIndex, const int size, double accordance = 1) :
        size(size),
        aSection(a, aIndex - size, aIndex - 1),
        bSection(b, bIndex - size, bIndex - 1),
        accordance(accordance) {
    }
    Match(const Text& a, const Text& b, const WeakMatch& wm) :
        size(wm.size),
        aSection(a, wm.i - wm.sizeI, wm.i - 1),
        bSection(b, wm.j - wm.sizeJ, wm.j - 1),
        accordance(wm.accordancePercentage()) {
    }
    inline bool operator<(const Match& other) const {return size > other.size;}
    int size;
    TextSection aSection;
    TextSection bSection;
    double accordance;
};

inline vector<Match> listMatches(const Text& a, const Text& b, const int minSize = 100) {
    Assert(minSize > 1);
    vector<Match> matches;
    if(a.trimText.empty() || b.trimText.empty()) return matches;

    int* curr = new int[b.trimText.size()];
    int* prev = new int[b.trimText.size()];

    for (size_t i = 0; i < a.trimText.size(); ++i) {
        for (size_t j = 0; j < b.trimText.size(); ++j) {
            if (a.trimText[i] != b.trimText[j]) {
                curr[j] = 0;
                if ((j != 0) && (i != 0) && (prev[j - 1] >= minSize)) {
                    matches.push_back(Match(a, b, i, j, prev[j - 1]));
                }
            } else {
                if(i == 0 || j == 0) {
                    curr[j] = 1;
                } else {
                    curr[j] = prev[j - 1] + 1;
                }
            }
        }
        if ((i != 0) && (prev[b.trimText.size() - 1] >= minSize)) {
            matches.push_back(Match(a, b, i, b.trimText.size(), prev[b.trimText.size() - 1]));
        }
        std::swap(curr, prev);
    }
    for (size_t j = minSize; j <= b.trimText.size(); ++j) {
        if (prev[j - 1] >= minSize) {
            matches.push_back(Match(a, b, a.trimText.size(), j, prev[j - 1]));
        }
    }

    delete[] curr;
    delete[] prev;

    return matches;
}

inline void combineWeakMatches(vector<WeakMatch>& matches, WeakMatch match, const double minAccordance = 0.8) {
    for (int i = matches.size() - 1; i >= 0; i--) {
        WeakMatch& m = matches[i];
        if (m < match) {
            int combineSize = m.size + m.dist(match) + match.size;
            int combineAccordance = m.accordance + match.accordance;
            if (combineSize * minAccordance < combineAccordance) {
                m.size = combineSize;
                m.sizeI = m.sizeI + match.i - m.i;
                m.sizeJ = m.sizeJ + match.j - m.j;
                m.accordance = combineAccordance;
                m.i = match.i;
                m.j = match.j;
                return;
            }
        }
    }
    matches.push_back(match);
}

inline vector<Match> listWeakMatches(const Text& a, const Text& b, const int minSize = 100, const double minAccordance = 0.8, const size_t minWordSize = 10) {
    Assert(minSize > 1);
    Assert(minWordSize > 1);
    vector<WeakMatch> weakMatches;
    if(a.trimText.empty() || b.trimText.empty()) return vector<Match>();

    size_t* curr = new size_t[b.trimText.size()];
    size_t* prev = new size_t[b.trimText.size()];

    for (size_t i = 0; i < a.trimText.size(); ++i) {
        for (size_t j = 0; j < b.trimText.size(); ++j) {
            if(a.trimText[i] != b.trimText[j]) {
                curr[j] = 0;
                if ((j != 0) && (i != 0) && (prev[j - 1] >= minWordSize)) {
                    combineWeakMatches(weakMatches, WeakMatch(i, j, prev[j - 1]), minAccordance);
                }
            } else {
                if(i == 0 || j == 0) {
                    curr[j] = 1;
                } else {
                    curr[j] = prev[j - 1] + 1;
                }
            }
        }
        if ((i >= minWordSize) && (prev[b.trimText.size() - 1] >= minWordSize)) {
            combineWeakMatches(weakMatches, WeakMatch(i, b.trimText.size(), prev[b.trimText.size() - 1]), minAccordance);
        }
        std::swap(curr, prev);
    }
    for (size_t j = minWordSize; j <= b.trimText.size(); ++j) {
        if (prev[j - 1] >= minWordSize) {
            combineWeakMatches(weakMatches, WeakMatch(a.trimText.size(), j, prev[j - 1]), minAccordance);
        }
    }

    delete[] curr;
    delete[] prev;

    vector<Match> matches;
    for (const WeakMatch& wm : weakMatches) {
        if (wm.size >= minSize) {
            matches.push_back(Match(a, b, wm));
        }
    }

    return matches;
}

inline double matchPercentage(const Text& text, const vector<Match>& matches) {
    int matchSize = 0;
    for (const Match& m : matches) matchSize += m.size;
    return matchSize / ((double) text.size());
}

template<typename T>
inline std::ostream& blue(const T& c) {return cout << "\033[1;34m" << c << "\033[0m";}

inline void printSection(const Text& text, const TextSection& section) {
    size_t lineBegin = size_t(std::max(int(0), int(section.TextPosBegin.line) - 1));
    size_t lineEnd = std::min(text.lines.size(), section.TextPosEnd.line + 1);
    cout << text.filename << " Line " << lineBegin << " - " << lineEnd << ":" << endl;
    cout << text.lines[lineBegin] << endl;
    if (section.TextPosBegin.line == section.TextPosEnd.line) {
        string sub = text.lines[section.TextPosBegin.line].substr(0, section.TextPosBegin.index);
        cout << sub;
        sub = text.lines[section.TextPosBegin.line].substr(section.TextPosBegin.index, section.TextPosEnd.index - section.TextPosBegin.index + 1);
        blue(sub);
        sub = text.lines[section.TextPosEnd.line].substr(section.TextPosEnd.index + 1);
        cout << sub << endl;
    } else {
        string sub = text.lines[section.TextPosBegin.line].substr(0, section.TextPosBegin.index);
        cout << sub;
        sub = text.lines[section.TextPosBegin.line].substr(section.TextPosBegin.index);
        blue(sub) << endl;
        for (size_t i = section.TextPosBegin.line + 1; i < section.TextPosEnd.line; i++) {
            blue(text.lines[i]) << endl;
        }
        sub = text.lines[section.TextPosEnd.line].substr(0, section.TextPosEnd.index + 1);
        blue(sub);
        sub = text.lines[section.TextPosEnd.line].substr(section.TextPosEnd.index + 1);
        cout << sub << endl;
    }
    cout << text.lines[lineEnd] << endl;
    cout << endl;
}

inline void printMatch(const Text& a, const Text& b, const Match& match) {
    if (match.accordance < 1) cout << "accordance: " << String::percent(match.accordance) << endl;
    printSection(a, match.aSection);
    printSection(b, match.bSection);
    cout << "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
    cout << endl;
}

template<typename SELECT_SECTION>
inline void printTextWithMatches(const Text& text, vector<Match>& matches, const SELECT_SECTION& selectSection) noexcept {
    cout << "\n\n\nFile: " << text.filename << ":\n" << endl;
    std::sort(matches.begin(), matches.end(), [&](const Match& a, const Match& b){
        return (selectSection(a).TextPosBegin < selectSection(b).TextPosBegin);
    });
    size_t matchI = 0;
    for (size_t textI = 0; textI < text.lines.size(); textI++) {
        string lineNumber = std::to_string(textI);
        for (size_t i = lineNumber.size(); i < 6; i++) {
            cout << " ";
        }
        cout << lineNumber << ": ";
        while ((matchI + 1 < matches.size()) && (selectSection(matches[matchI]).TextPosEnd.line < textI)) {
            matchI++;
        }
        TextSection section = selectSection(matches[matchI]);
        if ((section.TextPosBegin.line > textI) || (section.TextPosEnd.line < textI)) {
            cout << text.lines[textI] << endl;
        } else {
            size_t lineI = 0;
            while ((section.TextPosBegin.line <= textI) && (section.TextPosEnd.line >= textI)) {
                if (section.TextPosBegin.line == textI) {
                    string sub = text.lines[textI].substr(lineI, section.TextPosBegin.index);
                    cout << sub;
                    lineI = section.TextPosBegin.index;
                }
                if (section.TextPosEnd.line == textI) {
                    string sub = text.lines[textI].substr(lineI, section.TextPosEnd.index - lineI + 1);
                    blue(sub);
                    lineI = section.TextPosEnd.index + 1;
                    if (matchI + 1 < matches.size()) {
                        matchI++;
                        section = selectSection(matches[matchI]);
                    } else {
                        break;
                    }
                } else {
                    string sub = text.lines[textI].substr(lineI);
                    blue(sub);
                    lineI = text.lines[textI].size();
                    break;
                }
            }
            if (lineI < text.lines[textI].size()) {
                string sub = text.lines[textI].substr(lineI);
                cout << sub;
            }
            cout << endl;
        }
    }
}

int main(int argc, char **argv) {
    checkAsserts();

    CommandLineParser clp(argc, argv);
    if (clp.isSet("h")) usage(argv[0]);
    if (argc < 3) usage(argv[0]);

    int minWordLength = clp.value("w", 3);
    int minMatchLength = clp.value("m", 50);
    double minAccordance = clp.value("a", 1.0);
    bool printA = clp.isSet("pa");
    bool printB = clp.isSet("pb");

    string pdfFileA = argv[1];
    string pdfFileB = argv[2];
    cout << "Comparing " << pdfFileA << " to " << pdfFileB << endl;

    string textFileA = pdfFileA + ".txt";
    string textFileB = pdfFileB + ".txt";

    TextFile::pdfToText(pdfFileA, textFileA);
    TextFile::pdfToText(pdfFileB, textFileB);
    sleep(500);

    Text textA(textFileA, minWordLength);
    Text textB(textFileB, minWordLength);
    cout << textA.filename << ": " << textA.lines.size() << " lines" << endl;
    cout << textB.filename << ": " << textB.lines.size() << " lines" << endl;

    vector<Match> matches;
    if (minAccordance < 1) {
        matches = listWeakMatches(textA, textB, minMatchLength, minAccordance);
    } else {
        matches = listMatches(textA, textB, minMatchLength);
    }
    cout << "number of matches: " << matches.size() << endl << endl << endl;

    std::sort(matches.begin(), matches.end());

    for (const Match& m : matches) {
        printMatch(textA, textB, m);
    }

    if (printA) {
        printTextWithMatches(textA, matches, [](const Match& m){return m.aSection;});
    }

    if (printB) {
        printTextWithMatches(textB, matches, [](const Match& m){return m.bSection;});
    }

    cout << "accordance: \n    Text A: " << String::percent(matchPercentage(textA, matches)) << "\n    Text B: " << String::percent(matchPercentage(textB, matches)) << endl;

    return 0;
}
