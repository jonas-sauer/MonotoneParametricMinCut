#pragma once

#include <string>
#include <sstream>
#include <iostream>

class JSONString {

public:

    JSONString() :
        json(""),
        indentation(0),
        emptyElement(true) {
        json.precision(10);
    }

    void startNewObject() {
        addNewElement();
        indentation++;
        json << "{";
        newLine();
    }

    void endObject() {
        indentation--;
        newLine();
        json << "}";
        emptyElement = false;
    }

    void startNewArray() {
        addNewElement();
        indentation++;
        json << "[";
        newLine();
    }

    void endArray() {
        indentation--;
        newLine();
        json << "]";
        emptyElement = false;
    }

    void addQualifier(std::string name) {
        addNewElement();
        json << "\"" << name << "\":";
    }

    void addValue(std::string value) {
        json << "\"" << value << "\"";
        emptyElement = false;
    }

    void addValue(int value) {
        json << value;
        emptyElement = false;
    }

    void addValue(double value) {
        json << value;
        emptyElement = false;
    }

    void addQualifierWithValue(std::string name, std::string value) {
        addQualifier(name);
        addValue(value);
    }

    void addQualifierWithValue(std::string name, int value) {
        addQualifier(name);
        addValue(value);
    }

    void addQualifierWithValue(std::string name, double value) {
        addQualifier(name);
        addValue(value);
    }

    void newLine() {
        json << std::endl;
        addIndentation(indentation);
    }

    std::string str() {
        json << std::flush;
        return json.str();
    }

protected:

    void addNewElement() {
        if (!emptyElement) {
            json << ",";
            newLine();
            emptyElement = true;
        }
    }

    void addIndentation(int depth = 1) {
        for (int i = 0; i < depth; i++) {
            json << "    ";
        }
    }

    std::stringstream json;
    int indentation;
    bool emptyElement;

};
