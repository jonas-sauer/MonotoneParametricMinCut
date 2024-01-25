#pragma once

#include <iostream>
#include <algorithm>
#include <random>
#include <vector>
#include <string>

#include "../../Shell/Shell.h"

#include "../../Helpers/Assert.h"
#include "../../Helpers/Debug.h"
#include "../../Helpers/Timer.h"
#include "../../Helpers/Calendar.h"
#include "../../Helpers/MultiThreading.h"
#include "../../Helpers/IO/File.h"
#include "../../Helpers/IO/CSVData.h"
#include "../../Helpers/IO/ParserCSV.h"
#include "../../Helpers/Vector/Vector.h"
#include "../../Helpers/Ranges/Range.h"

#include "../../Helpers/String/String.h"
#include "../../Helpers/String/Enumeration.h"
#include "../../Helpers/String/TextFileUtils.h"
#include "../../Helpers/Vector/Vector.h"
#include "../../Helpers/Vector/Permutation.h"
#include "../../Helpers/FileSystem/FileSystem.h"

#include "../../DataStructures/Graph/Graph.h"
#include "../../DataStructures/GTFS/Entities/Calendar.h"
#include "../../DataStructures/GTFS/Data.h"
#include "../../DataStructures/Intermediate/Data.h"
#include "../../DataStructures/RAPTOR/Data.h"
#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/Demand/IdVertexDemand.h"
#include "../../DataStructures/Demand/PassengerData.h"

using namespace Shell;

class TikzPlot : public Command {

public:
    TikzPlot(BasicShell& shell) : shell(shell) {
        shell.addCommand(this);
    }

    virtual std::string name() const noexcept {
        return "tikzPlot";
    }

    virtual std::string helpText() const noexcept {
        return name() + " <input file> <output file>: generates a Tikz plot from a csv file.";
    }

    virtual void execute(const std::string& parameter) noexcept {
        std::vector<std::string> tokens = String::split(String::trim(parameter), ' ');
        if (tokens.size() < 1) tokens.emplace_back(shell.ask("CSV file> "));
        if (tokens[0] == "") return;
        IO::CSVData<double, IO::TrimChars<>, IO::DoubleQuoteEscape<',','"'>> data(tokens[0]);
        if (tokens.size() < 2) tokens.emplace_back(shell.ask("tex output file> "));
        if (tokens[1] == "") return;
        std::ofstream os(tokens[1]);
        Assert(os, "cannot open file: " << tokens[1]);
        Assert(os.is_open(), "cannot open file: " << tokens[1]);
        int minX = intMax;
        int maxX = -intMax;
        int minY = 0;
        int maxY = -intMax;
        for (const size_t r : range(data.numRows())) {
            minX = std::min(minX, static_cast<int>(data.columnData[0][r]));
            maxX = std::max(maxX, static_cast<int>(data.columnData[0][r] + 0.5));
            for (size_t c = 1; c < data.numColumns(); c++) {
                minY = std::min(minY, static_cast<int>(data.columnData[c][r]));
                maxY = std::max(maxY, static_cast<int>(data.columnData[c][r]));
            }
        }
        maxY = maxY * 1.1;
        maxY = maxY - (maxY % (maxY / 20));
        int stepX = std::max(1, static_cast<int>((maxX - minX) / 10 + 0.8));
        Enumeration xTicks;
        for (int x = minX; x <= maxX; x+= stepX) {
            xTicks << x << Sep(", ");
        }
        Enumeration legendEntries;
        for (size_t c = 1; c < data.numColumns(); c++) {
            legendEntries << data.columnNames[c] << Sep(", ");
        }
        os << "\\tikzstyle{markSign} = [mark=o]" << std::endl
           << "\\tikzstyle{shortenLines} = [shorten <= 3.5pt,shorten >= 3.5pt]" << std::endl
           << "" << std::endl
           << "\\begin{tikzpicture}" << std::endl
           << "\\pgfplotsset{" << std::endl
           << "   grid style = {dash pattern = on 1pt off 1pt, KITblack25,line width = 0.5pt  }" << std::endl
           << "}" << std::endl
           << "" << std::endl
           << "\\begin{axis}[" << std::endl
           << "   height=6.7cm," << std::endl
           << "   width=0.9\\textwidth," << std::endl
           << "   xmin=" << (minX - 0.5) << "," << std::endl
           << "   xmax=" << (maxX + 0.5) << "," << std::endl
           << "   ymin=" << minY << "," << std::endl
           << "   ymax=" << maxY << "," << std::endl
           << "   xlabel={" << data.columnNames[0] << "}," << std::endl
           << "   ylabel={Time [ms]}," << std::endl
           << "   xtick={" << xTicks.str() << "}," << std::endl
           << "%  xticklabel=\\pgfmathparse{\\tick}$2^{\\pgfmathprintnumber{\\pgfmathresult}}$," << std::endl
           << "   grid=major," << std::endl
           << "   legend entries={" << legendEntries.str() << "}," << std::endl
           << "   legend cell align=left," << std::endl
           << "   legend style={at={(0.05,0.92)}," << std::endl
           << "   anchor=north west," << std::endl
           << "   font=\\scriptsize}" << std::endl
           << "]" << std::endl
           << "" << std::endl;
        for (size_t c = 1; c < data.numColumns(); c++) {
            os << "\\addplot [color=plotColor" << c << ",line width=1.5pt] table {" << std::endl
               << "   " << (maxX + 100) << " 0" << std::endl
               << "   " << (maxX + 100) << " 1" << std::endl
               << "};" << std::endl;
        }
        os << "" << std::endl;
        for (size_t c = 1; c < data.numColumns(); c++) {
            for (size_t r = 1; r < data.numRows(); r++) {
                if (data.columnData[c][r - 1] != data.columnData[c][r - 1] || data.columnData[c][r] != data.columnData[c][r]) continue;
                os << "\\addplot [color=plotColor" << c << ",markSign,shortenLines,line width=1pt] table {" << std::endl
                   << "   " << data.columnData[0][r - 1] << " " << data.columnData[c][r - 1] << std::endl
                   << "   " << data.columnData[0][r] << " " << data.columnData[c][r] << std::endl
                   << "};" << std::endl;
            }
            os << "" << std::endl;
        }
        os << "\\end{axis}" << std::endl
           << "\\end{tikzpicture}" << std::endl
           << "" << std::endl;
    }

private:
    BasicShell& shell;

};

class AccumulateDijkstraRank : public Command {

public:
    AccumulateDijkstraRank(BasicShell& shell) {
        shell.addCommand(this);
    }

    virtual std::string name() const noexcept {
        return "accumulateDijkstraRank";
    }

    virtual std::string helpText() const noexcept {
        std::stringstream text;
        text << name() + " <rank column> <value column> <function> <output file> <input files>...: accumulates Dijkstra Rank raw data.\n"
             << "function:\n"
             << "   min\n"
             << "   mean\n"
             << "   median\n"
             << "   max\n"
             << "   sum";
        return text.str();
    }

    virtual void execute(const std::string& parameter) noexcept {
        std::vector<std::string> tokens = String::split(String::trim(parameter), ' ');
        if (tokens.size() < 5) return;
        if (tokens[2] == "min") {
            accumulateData(tokens, [](const std::vector<double>& data){
                if (data.empty()) return std::numeric_limits<double>::quiet_NaN();
                return Vector::min(data);
            });
        } else if (tokens[2] == "mean") {
            accumulateData(tokens, [](const std::vector<double>& data){
                if (data.empty()) return std::numeric_limits<double>::quiet_NaN();
                long double sum = 0;
                for (const double d : data) {
                    sum += d;
                }
                return static_cast<double>(sum / data.size());
            });
        } else if (tokens[2] == "median") {
            accumulateData(tokens, [](const std::vector<double>& data){
                if (data.empty()) return std::numeric_limits<double>::quiet_NaN();
                std::vector<double> sortedData = data;
                std::sort(sortedData.begin(), sortedData.end());
                if (sortedData.size() % 2 != 0) {
                    return sortedData[sortedData.size() / 2];
                } else {
                    return 0.5 * (sortedData[sortedData.size() / 2] + sortedData[(sortedData.size() / 2) - 1]);
                }
            });

        } else if (tokens[2] == "max") {
            accumulateData(tokens, [](const std::vector<double>& data){
                if (data.empty()) return std::numeric_limits<double>::quiet_NaN();
                return Vector::max(data);
            });
        } else if (tokens[2] == "sum") {
            accumulateData(tokens, [](const std::vector<double>& data){
                if (data.empty()) return std::numeric_limits<double>::quiet_NaN();
                double sum = 0;
                for (const double d : data) {
                    sum += d;
                }
                return sum;
            });
        } else {
            error("The function " + tokens[2] + " is not supported!");
        }
    }

private:
    template<typename ACCUMULATE>
    inline void accumulateData(const std::vector<std::string>& tokens, const ACCUMULATE& accumulate) const noexcept {
        std::vector<std::string> fileNames = getFileNames(tokens);
        std::vector<double> ranks;
        Map<double, std::vector<std::vector<double>>> data;
        for (size_t i = 0; i < fileNames.size(); i++) {
            IO::CSVReader<2, IO::TrimChars<>, IO::DoubleQuoteEscape<'\t','"'>> in(fileNames[i]);
            in.readHeader(tokens[0], tokens[1]);
            double rank;
            double value;
            while (in.readRow(rank, value)) {
                if (!data.contains(rank)) {
                    data.insert(rank, std::vector<std::vector<double>>(fileNames.size(), std::vector<double>()));
                    ranks.emplace_back(rank);
                }
                data[rank][i].emplace_back(value);
            }
        }
        std::sort(ranks.begin(), ranks.end());
        IO::CSVData<double, IO::TrimChars<>, IO::DoubleQuoteEscape<',','"'>> result;
        result.columnNames = getColumnNames(tokens);
        result.columnData = std::vector<std::vector<double>>(result.columnNames.size(), std::vector<double>());
        result.columnData[0] = ranks;
        for (const double rank : ranks) {
            for (size_t i = 0; i < fileNames.size(); i++) {
                result.columnData[i + 1].emplace_back(accumulate(data[rank][i]));
            }
        }
        result.write(tokens[3]);
    }

    inline std::vector<std::string> getFileNames(const std::vector<std::string>& tokens) const noexcept {
        std::vector<std::string> result;
        for (size_t i = 4; i < tokens.size(); i++) {
            result.emplace_back(tokens[i]);
        }
        return result;
    }

    inline std::vector<std::string> getColumnNames(const std::vector<std::string>& tokens) const noexcept {
        std::vector<std::string> result;
        result.emplace_back(tokens[0]);
        for (size_t i = 4; i < tokens.size(); i++) {
            result.emplace_back(FileSystem::getFileNameWithoutExtension(tokens[i]));
        }
        return result;
    }

};

class DijkstraRankPlot : public Command {

public:
    DijkstraRankPlot(BasicShell& shell) {
        shell.addCommand(this);
    }

    virtual std::string name() const noexcept {
        return "dijkstraRankPlot";
    }

    virtual std::string helpText() const noexcept {
        return name() + " <input file> <output file> <rank column> <value columns>...: Creates a Dijkstra Rank Plot in tikz/pgf format.\n";
    }

    virtual void execute(const std::string& parameter) noexcept {
        std::vector<std::string> tokens = String::split(String::trim(parameter), ' ');
        if (tokens.size() < 4) return;
        IO::CSVData<double, IO::TrimChars<>, IO::DoubleQuoteEscape<',','"'>> data;
        data.read(tokens[0]);
        std::vector<std::vector<std::vector<double>>> dataByColumnByRank;
        std::vector<double>& ranks = data.getColumn(tokens[2]);
        std::vector<std::string> columnNames;
        std::vector<std::vector<double>> columns;
        size_t minX = intMax;
        size_t maxX = 0;
        int minY = intMax;
        int maxY = -intMax;
        for (size_t i = 3; i < tokens.size(); i++) {
            columnNames.emplace_back(tokens[i]);
            columns.emplace_back(data.getColumn(tokens[i]));
        }
        for (size_t row = 0; row < data.numRows(); row++) {
            const size_t rank = ranks[row];
            minX = std::min(minX, rank);
            maxX = std::max(maxX, rank);
            if (dataByColumnByRank.size() <= rank) dataByColumnByRank.resize(rank + 1, std::vector<std::vector<double>>(columns.size()));
            for (size_t column = 0; column < columns.size(); column++) {
                dataByColumnByRank[rank][column].emplace_back(columns[column][row]);
                const int value = columns[column][row];
                minY = std::min(minY, value);
                maxY = std::max(maxY, value);
            }
        }
        Enumeration xTicks;
        for (size_t i = minX; i <= maxX; i++) {
            xTicks << i << sep;
        }
        const double boxWidth = 1 / ((1.7 * columns.size()) + 0.3);
        std::ofstream os(tokens[1]);
        Assert(os, "cannot open file: " << tokens[1]);
        Assert(os.is_open(), "cannot open file: " << tokens[1]);
        os << "\\newcommand{\\plotHeight}{6.5cm}\n"
           << "\\newcommand{\\plotWidth}{0.951\\textwidth}\n"
           << "\\newcommand{\\boxWidth}{" << boxWidth << "}\n"
           << "\\newcommand{\\minValueX}{" << minX << "}\n"
           << "\\newcommand{\\maxValueX}{" << maxX << "}\n"
           << "\\newcommand{\\minValueY}{" << minY << "}\n"
           << "\\newcommand{\\maxValueY}{" << maxY << "}\n"
           << "\\newcommand{\\markSign}{x}\n"
           << "\\newcommand{\\markSize}{1.5}\n"
           << "\\newcommand{\\lineWidth}{1.0pt}\n"
           << "\n"
           << "\\begin{tikzpicture}\n"
           << "\\colorlet{lineColor1}{KITseablue}\n"
           << "\\colorlet{lineColor2}{KITgreen}\n"
           << "\\colorlet{lineColor3}{KITred}\n"
           << "\\colorlet{lineColor4}{KITorange}\n"
           << "\\colorlet{lineColor5}{KITcyanblue}\n"
           << "\\colorlet{lineColor6}{KITpalegreen}\n"
           << "\\colorlet{lineColor7}{KITlilac}\n"
           << "\\colorlet{lineColor8}{KITblue}\n"
           << "\\colorlet{fillColor1}{KITseablue30}\n"
           << "\\colorlet{fillColor2}{KITgreen30}\n"
           << "\\colorlet{fillColor3}{KITred30}\n"
           << "\\colorlet{fillColor4}{KITorange30}\n"
           << "\\colorlet{fillColor5}{KITcyanblue30}\n"
           << "\\colorlet{fillColor6}{KITpalegreen30}\n"
           << "\\colorlet{fillColor7}{KITlilac30}\n"
           << "\\colorlet{fillColor8}{KITblue30}\n"
           << "\\colorlet{markColor1}{KITseablue70}\n"
           << "\\colorlet{markColor2}{KITgreen70}\n"
           << "\\colorlet{markColor3}{KITred70}\n"
           << "\\colorlet{markColor4}{KITorange70}\n"
           << "\\colorlet{markColor5}{KITcyanblue70}\n"
           << "\\colorlet{markColor6}{KITpalegreen70}\n"
           << "\\colorlet{markColor7}{KITlilac70}\n"
           << "\\colorlet{markColor8}{KITblue70}\n"
           << "\n"
           << "\n"
           << "\\pgfplotsset{\n"
           << "   grid style={KITblack20,line width = 0.2pt,dash pattern = on 2pt off 1pt},\n"
           << "   xtick style={line width=0.4pt,black},\n"
           << "   ytick style={line width=0.4pt,black},\n"
           << "   boxplot/every median/.style={black, line width = 1.5pt},\n"
           << "}\n"
           << "\n"
           << "\\begin{axis}[\n"
           << "   height=\\plotHeight,\n"
           << "   width=\\plotWidth,\n"
           << "   grid=both,\n"
           << "   xmin=\\minValueX,\n"
           << "   xmax=\\maxValueX + 1,\n"
           << "   xtick align=outside,\n"
           << "   xtick pos=left,\n"
           << "   xtick={" << xTicks << ", " << (maxX + 1) << "},\n"
           << "   xticklabel=\\pgfmathparse{(\\tick)}${}$,\n"
           << "   ymin=\\minValueY,\n"
           << "   ymax=\\maxValueY,\n"
           << "   ymode=log,\n"
           << "   yticklabel=\\pgfmathparse{(\\tick)}${}$\n"
           << "]\n"
           << "\\end{axis}\n"
           << "\n"
           << "\\begin{axis}[\n"
           << "   height=\\plotHeight,\n"
           << "   width=\\plotWidth,\n"
           << "   ticklabel pos=right,\n"
           << "   xmin=\\minValueX,\n"
           << "   xmax=\\maxValueX + 1,\n"
           << "   xtick style={draw=none},\n"
           << "   xticklabel=\\pgfmathparse{(\\tick)}${}$,\n"
           << "   ymin=\\minValueY,\n"
           << "   ymax=\\maxValueY,\n"
           << "   ymode=log,\n"
           << "   yminorticks=false,\n"
           << "   ytick pos=right,\n"
           << "   ytick align=outside\n"
           << "]\n"
           << "\\end{axis}\n"
           << "\n"
           << "\\begin{axis}[\n"
           << "   height=\\plotHeight,\n"
           << "   width=\\plotWidth,\n"
           << "   boxplot/draw direction=y,\n"
           << "   xmin=\\minValueX - 0.5,\n"
           << "   xmax=\\maxValueX + 0.5,\n"
           << "   xtick style={draw=none},\n"
           << "   xtick={" << xTicks << "},\n"
           << "   ymin=\\minValueY,\n"
           << "   ymax=\\maxValueY,\n"
           << "   ymode=log,\n"
           << "   yminorticks=false,\n"
           << "   ytick pos=left,\n"
           << "   ytick align=outside\n"
           << "]\n"
           << "\n";
        for (size_t rank = minX; rank <= maxX; rank++) {
            for (size_t column = 0; column < columns.size(); column++) {
                std::sort(dataByColumnByRank[rank][column].begin(), dataByColumnByRank[rank][column].end());
                const double x = rank + (((1.7 * column) + 1) * boxWidth) - 0.5;
                const double median = percentile(dataByColumnByRank[rank][column], 0.5);
                const double lowerQuartile = percentile(dataByColumnByRank[rank][column], 0.25);
                const double upperQuartile = percentile(dataByColumnByRank[rank][column], 0.75);
                const double iqr = 1.5 * (upperQuartile - lowerQuartile);
                double lowerWhisker = lowerQuartile;
                double upperWhisker = upperQuartile;
                std::vector<double> points;
                for (const double d : dataByColumnByRank[rank][column]) {
                    if (d < lowerWhisker) {
                        if (d >= lowerQuartile - iqr) {
                            lowerWhisker = d;
                        } else {
                            points.emplace_back(d);
                        }
                    }
                }
                for (const double d : dataByColumnByRank[rank][column]) {
                    if (d > upperWhisker) {
                        if (d <= upperQuartile + iqr) {
                            upperWhisker = d;
                        } else {
                            points.emplace_back(d);
                        }
                    }
                }
                os << "\\addplot[\n"
                   << "   mark=\\markSign,\n"
                   << "   mark size=\\markSize,\n"
                   << "   mark options={markColor" << (column + 1) << "},\n"
                   << "   boxplot prepared={\n"
                   << "      draw position=" << x << ",\n"
                   << "      lower whisker=" << lowerWhisker << ",\n"
                   << "      lower quartile=" << lowerQuartile << ",\n"
                   << "      median=" << median << ",\n"
                   << "      upper quartile=" << upperQuartile << ",\n"
                   << "      upper whisker=" << upperWhisker << ",\n"
                   << "      box extend=\\boxWidth,\n"
                   << "      every box/.style={lineColor" << (column + 1) << ",fill=fillColor" << (column + 1) << ",line width = \\lineWidth},\n"
                   << "      every whisker/.style={lineColor" << (column + 1) << ",line width = \\lineWidth},\n"
                   << "   },\n"
                   << "] coordinates {\n"
                   << "   ";
                for (const double p : points) {
                    os << "(0," << p << ") ";
                }
                os << "\n"
                   << "};\n"
                   << "\n";
            }
        }
        os << "\\end{axis}\n"
           << "\n"
           << "\\end{tikzpicture}";
    }

private:
    inline double percentile(const std::vector<double>& sortedData, const double p) const noexcept {
        Assert(!sortedData.empty(), "Percentile is not defined for empty data sets!");
        Assert(p >= 0, "Percentile cannot be negative!");
        Assert(p <= 1, "Percentile cannot be greater than one!");
        if (sortedData.size() == 1) return sortedData.front();
        const double index = (sortedData.size() - 1) * p;
        const size_t lowerIndex = index;
        const size_t higherIndex = lowerIndex + 1;
        if (higherIndex == sortedData.size()) return sortedData.back();
        const double lambda = higherIndex - index;
        return (lambda * sortedData[lowerIndex]) + ((1 - lambda) * sortedData[higherIndex]);
    }

};

class BoxCode : public ParameterizedCommand {

public:
    BoxCode(BasicShell& shell) :
        ParameterizedCommand(shell, "boxCode", "Generates box and whiskers code for PGF from CSV.") {
        addParameter("CSV file");
        addParameter("X");
        addParameter("Color");
    }

    virtual void execute() noexcept {
        const std::string csvFile = getParameter("CSV file");
        const std::string x = getParameter("X");
        const std::string color = getParameter("Color");

        IO::CSVReader<1, IO::TrimChars<>, IO::DoubleQuoteEscape<',','"'>> in(csvFile);
        in.readHeader(IO::IGNORE_EXTRA_COLUMN, "QueryTime");
        std::vector<double> data;
        double value;
        while (in.readRow(value)) {
            data.emplace_back(value);
        }

        std::sort(data.begin(), data.end());
        const double median = Vector::percentile(data, 0.5);
        const double lowerQuartile = Vector::percentile(data, 0.25);
        const double upperQuartile = Vector::percentile(data, 0.75);
        const double iqr = 1.5 * (upperQuartile - lowerQuartile);
        double lowerWhisker = lowerQuartile;
        double upperWhisker = upperQuartile;
        std::vector<double> points;
        for (const double d : data) {
            if (d < lowerWhisker) {
                if (d >= lowerQuartile - iqr) {
                    lowerWhisker = d;
                } else {
                    points.emplace_back(d);
                }
            }
        }
        for (const double d : data) {
            if (d > upperWhisker) {
                if (d <= upperQuartile + iqr) {
                    upperWhisker = d;
                } else {
                    points.emplace_back(d);
                }
            }
        }
        std::cout << "\\addplot[\n"
                  << "   mark=\\markSign,\n"
                  << "   mark size=\\markSize,\n"
                  << "   mark options={markColor" << color << "},\n"
                  << "   boxplot prepared={\n"
                  << "      draw position=" << x << ",\n"
                  << "      lower whisker=" << lowerWhisker << ",\n"
                  << "      lower quartile=" << lowerQuartile << ",\n"
                  << "      median=" << median << ",\n"
                  << "      upper quartile=" << upperQuartile << ",\n"
                  << "      upper whisker=" << upperWhisker << ",\n"
                  << "      box extend=\\boxWidth,\n"
                  << "      every box/.style={lineColor" << color << ",fill=fillColor" << color << ",line width = \\lineWidth},\n"
                  << "      every whisker/.style={lineColor" << color << ",line width = \\lineWidth},\n"
                  << "   },\n"
                  << "] coordinates {\n"
                  << "   ";
        for (const double p : points) {
            std::cout << "(0," << p << ") ";
        }
        std::cout << "\n"
                  << "};\n"
                  << "\n";

    }

};

class AnalyzeData : public Command {

public:
    AnalyzeData(BasicShell& shell) {
        shell.addCommand(this);
    }

    virtual std::string name() const noexcept {
        return "analyzeData";
    }

    virtual std::string helpText() const noexcept {
        return name() + " <value column> <output file> <input files>...: accumulates Dijkstra Rank raw data.\n";
    }

    virtual void execute(const std::string& parameter) noexcept {
        std::vector<std::string> tokens = String::split(String::trim(parameter), ' ');
        if (tokens.size() < 3) return;
        IO::CSVData<std::string, IO::TrimChars<>, IO::DoubleQuoteEscape<',','"'>> result;
        result.columnNames = std::vector<std::string>{"file", "min", "mean", "median", "max", "sum"};
        result.columnData = std::vector<std::vector<std::string>>(result.columnNames.size(), std::vector<std::string>());
        std::vector<std::string> fileNames = getFileNames(tokens);
        for (size_t i = 0; i < fileNames.size(); i++) {
            double sum = 0;
            std::vector<double> data;
            IO::CSVReader<1, IO::TrimChars<>, IO::DoubleQuoteEscape<'\t','"'>> in(fileNames[i]);
            in.readHeader(tokens[0]);
            double value;
            while (in.readRow(value)) {
                data.emplace_back(value);
                sum += value;
            }
            std::sort(data.begin(), data.end());
            const double median = (data.size() % 2 != 0) ? (data[data.size() / 2]) : (0.5 * (data[data.size() / 2] + data[(data.size() / 2) - 1]));
            result.columnData[0].emplace_back(FileSystem::getFileNameWithoutExtension(fileNames[i]));
            result.columnData[1].emplace_back(std::to_string(data.front()));
            result.columnData[2].emplace_back(std::to_string(sum / data.size()));
            result.columnData[3].emplace_back(std::to_string(median));
            result.columnData[4].emplace_back(std::to_string(data.back()));
            result.columnData[5].emplace_back(std::to_string(sum));
        }
        result.write(tokens[1]);
    }

private:
    inline std::vector<std::string> getFileNames(const std::vector<std::string>& tokens) const noexcept {
        std::vector<std::string> result;
        for (size_t i = 2; i < tokens.size(); i++) {
            result.emplace_back(tokens[i]);
        }
        return result;
    }

};

class AnalyzeHeuristics : public Command {

public:
    AnalyzeHeuristics(BasicShell& shell) {
        shell.addCommand(this);
    }

    virtual std::string name() const noexcept {
        return "analyzeHeuristics";
    }

    virtual std::string helpText() const noexcept {
        return name() + " <value column> <output file> <input files>...: accumulates Dijkstra Rank raw data.\n";
    }

    virtual void execute(const std::string& parameter) noexcept {
        std::vector<std::string> tokens = String::split(String::trim(parameter), ' ');
        if (tokens.size() < 2) return;
        std::vector<std::string> fileNames = getFileNames(tokens);
        std::vector<IO::CSVData<double, IO::TrimChars<>, IO::DoubleQuoteEscape<'\t','"'>>> files;
        for (const std::string& fileName : fileNames) {
            files.emplace_back(fileName);
        }
        IO::OFStream result = IO::OFStream(tokens[0]);
        const size_t rows = files[0].numRows();
        const std::vector<double>& optimum = files[0].getColumn("travelTime");
        for (size_t i = 0; i < fileNames.size(); i++) {
            double feasible = 0;
            double optimal = 0;
            std::vector<double> deviation;
            const std::vector<double>& labels = files[i].getColumn("q_relaxed_edges");
            const std::vector<double>& domination = files[i].getColumn("q_domination");
            const std::vector<double>& travelTime = files[i].getColumn("travelTime");
            const std::vector<double>& time = files[i].getColumn("totalTime");
            for (size_t j = 0; j < rows; j++) {
                if (travelTime[j] >= 0) {
                    feasible++;
                    if (travelTime[j] != optimum[j]) {
                        deviation.emplace_back(std::max(travelTime[j] / optimum[j], optimum[j] / travelTime[j]));
                    } else {
                        deviation.emplace_back(1);
                        optimal++;
                    }
                }
            }
            result << FileSystem::getFileNameWithoutExtension(fileNames[i]) << " ";
            result << Vector::mean(labels) << " ";
            result << Vector::mean(domination) << " ";
            result << Vector::mean(time) << " ";
            result << (feasible / rows) << " ";
            result << (optimal / feasible) << " ";
            result << Vector::mean(deviation) << " ";
            result << Vector::max(deviation) << "\n";
        }
    }

private:
    inline std::vector<std::string> getFileNames(const std::vector<std::string>& tokens) const noexcept {
        std::vector<std::string> result;
        for (size_t i = 1; i < tokens.size(); i++) {
            result.emplace_back(tokens[i]);
        }
        return result;
    }

};

class BuildDfgScript : public Command {

public:
    BuildDfgScript(BasicShell& shell) : shell(shell) {
        shell.addCommand(this);
    }

    virtual std::string name() const noexcept {
        return "buildDfgScript";
    }

    virtual std::string helpText() const noexcept {
        return name() + " <input folder> <settings file>: Generates an Group assignment script for all files in the folder.\n";
    }

    virtual void execute(const std::string& parameter) noexcept {
        std::vector<std::string> tokens = String::split(String::trim(parameter), ' ');
        if (tokens.size() < 1) tokens.emplace_back(shell.ask("input folder> "));
        if (tokens[0] == "") return;
        if (tokens.size() < 2) tokens.emplace_back(shell.ask("settings file> "));
        if (tokens[1] == "") return;
        const std::string baseFolder = tokens[0];
        const std::string binaryDataFolder = FileSystem::extendPath(baseFolder, "CSA_Assignment_Binary_Data/");
        const std::string resultsFolder = FileSystem::extendPath(baseFolder, "CSA_Assignment_Results/");
        const std::string settingsFolder = FileSystem::extendPath(baseFolder, "CSA_Assignment_Settings/");
        std::ofstream createBinaryNetworksScript(FileSystem::extendPath(settingsFolder, "createBinaryNetworks.script"));
        std::ofstream computeAssignmentsScript(FileSystem::extendPath(settingsFolder, "computeAssignments.script"));
        std::vector<std::string> files = FileSystem::getFiles(baseFolder);
        for (const std::string& file : files) {
            if (!FileSystem::isFile(FileSystem::extendPath(baseFolder, file))) continue;
            if (!(String::endsWith(file, "connection.csv") || String::endsWith(file, "connections.csv"))) continue;
            const std::string identifier = file.substr(0, String::lastIndexOf(file, "connection"));
            const std::string dataFiles = FileSystem::extendPath(baseFolder, identifier);
            const std::string demandFiles = FileSystem::extendPath(baseFolder, identifier + "demand.csv");
            const std::string binaryFiles = FileSystem::extendPath(binaryDataFolder, identifier + "csa.binary");
            const std::string resultFiles = FileSystem::extendPath(resultsFolder, identifier);
            const std::string resultAggregateFile = FileSystem::extendPath(resultsFolder, "aggregate.csv");
            const std::string defaultSettingsFiles = FileSystem::extendPath(settingsFolder, tokens[1]);
            const std::string parameters = String::replaceAll(String::replaceAll(String::replaceAll(String::replaceAll(identifier, "FPZ=", ""), "__ZI=", ","), "__dAbf+", ","), "min__", "");
            createBinaryNetworksScript << "parseCSV " << dataFiles << " " << binaryFiles << " true false true" << std::endl;
            computeAssignmentsScript << "groupAssignment " << defaultSettingsFiles << " " << binaryFiles << " " << demandFiles << " " << resultFiles << " 0 1 " << resultAggregateFile << " " << parameters << std::endl;
        }
    }

private:
    BasicShell& shell;

};
