#pragma once

#include <vector>
#include <string>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "ParserCSV.h"

#include "../Assert.h"
#include "../String/String.h"

namespace IO {

template<class DATA_TYPE = std::string, class TRIM_POLICY = TrimChars<>, class QUOTE_POLICY = NoQuoteEscape<','>, class OVERFLOW_POLICY = ThrowOnOverflow, class COMMENT_POLICY = EmptyLineComment>
class CSVData {

public:
    using DataType = DATA_TYPE;
    using TrimPolics = TRIM_POLICY;
    using QuotePolicy = QUOTE_POLICY;
    using OverflowPolicy = OVERFLOW_POLICY;
    using CommentPolicy = COMMENT_POLICY;
    using Type = CSVData<DataType, TrimPolics, QuotePolicy, OverflowPolicy, CommentPolicy>;

public:
    CSVData() {}
    CSVData(const std::string& fileName) {
        read(fileName);
    }
    CSVData(const std::vector<std::string>& columnNames) : columnNames(columnNames) {
        columnData = std::vector<std::vector<DATA_TYPE>>(columnNames.size(), std::vector<DATA_TYPE>());
    }

    inline size_t numRows() const noexcept {
        return (columnData.empty()) ? (0) : (columnData[0].size());
    }

    inline size_t numColumns() const noexcept {
        return columnNames.size();
    }

    inline size_t getColumnIndex(const std::string& columnName) const noexcept {
        for (size_t i = 0; i < columnNames.size(); i++) {
            if (columnNames[i] == columnName) return i;
        }
        error("Data does not contain a column named " + columnName + "!");
        exit(1);
    }

    inline std::vector<DataType>& getColumn(const std::string& columnName) noexcept {
        return columnData[getColumnIndex(columnName)];
    }

    inline const std::vector<DataType>& getColumn(const std::string& columnName) const noexcept {
        return columnData[getColumnIndex(columnName)];
    }

    inline void appendRow(const std::vector<DataType>& row) noexcept {
        Assert(row.size() == numColumns(), "Cannot append a row of length " << row.size() << " to a table of width " << numColumns() << "!");
        for (size_t i = 0; i < row.size(); i++) {
            columnData[i].emplace_back(row[i]);
        }
    }

    inline void read(const std::string& fileName) noexcept {
        LineReader in(fileName);
        try {
            char* line = nullptr;
            do {
                line = in.nextLine();
                if (!line) throw Error::HeaderMissing();
            } while(COMMENT_POLICY::isComment(line));
            columnNames = Detail::parseLine<TRIM_POLICY, QUOTE_POLICY>(line);
            columnData.resize(columnNames.size());
            while (true) {
                line = in.nextLine();
                if (!line) break;
                if (COMMENT_POLICY::isComment(line)) continue;
                std::vector<std::string> rowData = Detail::parseLine<TRIM_POLICY, QUOTE_POLICY>(line);
                for (size_t i = 0; i < rowData.size() && i < columnData.size(); i++) {
                    DataType value;
                    Detail::parse<OVERFLOW_POLICY>(rowData[i].c_str(), value);
                    columnData[i].emplace_back(value);
                }
            }
        } catch(Error::WithFileName& error) {
            error.setFileName(in.getTruncatedFileName());
        }
    }

    inline void write(std::ostream& os) const noexcept {
        Assert(os, "cannot write to stream");
        if (!columnNames.empty()) {
            os << columnNames.front();
            for (size_t i = 1; i < columnNames.size(); i++) {
                os << QUOTE_POLICY::Sep << columnNames[i];
            }
            os << std::endl;
        }
        if (!columnData.empty()) {
            for (size_t j = 0; j < columnData[0].size(); j++) {
                os << columnData[0][j];
                for (size_t i = 1; i < columnData.size(); i++) {
                    if (columnData[i].size() > j) {
                        os << QUOTE_POLICY::Sep << columnData[i][j];
                    }
                }
                os << std::endl;
            }
        }
    }

    inline void write(const std::string& fileName) const noexcept {
        std::ofstream os(fileName);
        Assert(os, "cannot open file: " << fileName);
        Assert(os.is_open(), "cannot open file: " << fileName);
        write(os);
    }

    inline size_t columnWidth(const size_t index) const noexcept {
        size_t result = columnNames[index].size();
        for (const DataType& d : columnData[index]) {
            std::stringstream ss;
            ss << d;
            const size_t size = ss.str().size();
            if (size > result) result = size;
        }
        return result;
    }

    inline size_t columnWidth(const std::string& columnName) const noexcept {
        return columnWidth(getColumnIndex(columnName));
    }

    inline std::vector<size_t> columnWidth() const noexcept {
        std::vector<size_t> result;
        for (size_t i = 0; i < columnNames.size(); i++) {
            result.emplace_back(columnWidth(i));
        }
        return result;
    }

    inline void print() const noexcept {
        std::vector<size_t> width = columnWidth();
        if (!columnNames.empty()) {
            std::cout << std::setw(width[0] + 2) << columnNames.front();
            for (size_t i = 1; i < columnNames.size(); i++) {
                std::cout << std::setw(width[i] + 2) << columnNames[i];
            }
            std::cout << std::endl;
        }
        if (!columnData.empty()) {
            for (size_t j = 0; j < columnData[0].size(); j++) {
                std::cout << std::setw(width[0] + 2) << columnData[0][j];
                for (size_t i = 1; i < columnData.size(); i++) {
                    if (columnData[i].size() > j) {
                        std::cout << std::setw(width[i] + 2) << columnData[i][j];
                    }
                }
                std::cout << std::endl;
            }
        }
    }

public:
    std::vector<std::string> columnNames;
    std::vector<std::vector<DataType>> columnData;

};

}
