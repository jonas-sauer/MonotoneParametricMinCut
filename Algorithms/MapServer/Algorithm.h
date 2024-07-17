#pragma once

#include <map>
#include <vector>

#include "../../Helpers/Types.h"

#include "../../DataStructures/Result.h"
#include "../../DataStructures/Parameter.h"

class Algorithm {

public:
    Algorithm() {}

    virtual ~Algorithm() {}

    virtual std::vector<Result> run(const Vertex sourceVertexId, const Vertex targetNodeId) = 0;

    virtual std::vector<Result> runSourceOnly(const Vertex sourceVertexId) noexcept {
        std::vector<Result> resultSet;
        resultSet.push_back(Result("Source Node", KIT::green));
        resultSet.back().nodes.push_back(Result::Node(sourceVertexId, "S"));
        return resultSet;
    }

    virtual std::vector<Result> runTargetOnly(const Vertex targetNodeId) noexcept {
        std::vector<Result> resultSet;
        resultSet.push_back(Result("Target Node", KIT::green));
        resultSet.back().nodes.push_back(Result::Node(targetNodeId, "T"));
        return resultSet;
    }

    inline const std::map<std::string, Parameter>& getParameters() const noexcept {
        return parameters;
    }

    inline std::map<std::string, Parameter>& getParameters() noexcept {
        return parameters;
    }

protected:
    inline void addParameter(Parameter p) noexcept {
        parameters[p.name] = p;
    }

    inline void constParameter(const std::string& name, const std::string& value) noexcept {
        std::map<std::string, Parameter>::iterator position = parameters.find(name);
        if (position != parameters.end()) {
            constParameters[position->first] = position->second;
            constParameters[position->first].defaultValue = value;
            parameters.erase(position);
        }
    }

    template <typename T>
    inline T getParameter(const std::string& name) noexcept {
        std::map<std::string, Parameter>::iterator position = parameters.find(name);
        if (position != parameters.end()) {
            return position->second.getValue<T>();
        } else {
            position = constParameters.find(name);
            if (position != constParameters.end()) {
                return position->second.getValue<T>();
            } else {
                std::cout << "WARNING! Parameter " << name << " not found, returning: " << T() << std::endl;
                return T();
            }
        }
    }

    template <typename T>
    inline void setParameter(const std::string& name, const T& value) noexcept {
        std::map<std::string, Parameter>::iterator position = parameters.find(name);
        if (position != parameters.end()) {
            position->second.value = std::to_string(value);
        }
    }

    std::map<std::string, Parameter> parameters;

private:
    std::map<std::string, Parameter> constParameters;

};
