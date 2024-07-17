#pragma once

#include <string>
#include <vector>
#include <iostream>

#include "../Container/Map.h"
#include "../Geometry/Point.h"
#include "../Geometry/Rectangle.h"
#include "../Geometry/CoordinateTree.h"

#include "../../Helpers/Vector/Vector.h"
#include "../../Helpers/String/String.h"
#include "../../Helpers/IO/Serialization.h"
#include "../../Visualization/MapVisualization.h"

namespace BikeSharing {

class Stations {

public:
    Stations() {}

public:
    static Stations FromBinary(const std::string& fileName) {
        Stations result;
        result.deserialize(fileName);
        return result;
    }

    static Stations FromIntermediateFormat(const std::string& fileName, const double coordinateFactor = 0.000001) {
        Stations result;
        result.readIntermediateFile(fileName, coordinateFactor);
        return result;
    }

public:
    inline void deleteOperator(const std::string& operatorName) noexcept {
        if (!map.contains(operatorName)) return;
        map.erase(operatorName);
    }

    inline void deleteSmallOperators(const size_t maximumDeletedStations = 1, const std::vector<std::string>& exceptions = std::vector<std::string>()) noexcept {
        for (const std::string& operatorName : getOperators()) {
            if (Vector::contains(exceptions, operatorName)) continue;
            if (map[operatorName].size() > maximumDeletedStations) continue;
            map.erase(operatorName);
        }
    }

    inline void reassignSmallOperators(const size_t maximumReassignedStations = 1, const std::vector<std::string>& exceptions = std::vector<std::string>()) noexcept {
        std::vector<std::string> operatorNames;
        std::vector<Geometry::Point> coordinates;
        for (auto [key, value] : map) {
            if (Vector::contains(exceptions, key)) continue;
            if (value.size() <= maximumReassignedStations) continue;
            for (const Geometry::Point& coordinate : value) {
                operatorNames.emplace_back(key);
                coordinates.emplace_back(coordinate);
            }
        }
        Geometry::Rectangle boundingBox = Geometry::Rectangle::BoundingBox(coordinates);
        Geometry::GeoMetricAproximation metric = Geometry::GeoMetricAproximation::ComputeCorrection(boundingBox.center());
        CoordinateTree<Geometry::GeoMetricAproximation> ct(metric, coordinates);
        for (const std::string& operatorName : getOperators()) {
            if (Vector::contains(exceptions, operatorName)) continue;
            if (map[operatorName].size() > maximumReassignedStations) continue;
            for (const Geometry::Point& coordinate : map[operatorName]) {
                map[operatorNames[ct.getNearestNeighbor(coordinate)]].emplace_back(coordinate);
            }
            map.erase(operatorName);
        }
    }

    inline void reassignStationsOfOperator(const std::string& operatorName, const std::vector<std::string>& exceptions = std::vector<std::string>()) noexcept {
        if (!map.contains(operatorName)) return;
        std::vector<std::string> operatorNames;
        std::vector<Geometry::Point> coordinates;
        for (auto [key, value] : map) {
            if (key == operatorName) continue;
            if (Vector::contains(exceptions, key)) continue;
            for (const Geometry::Point& coordinate : value) {
                operatorNames.emplace_back(key);
                coordinates.emplace_back(coordinate);
            }
        }
        Geometry::Rectangle boundingBox = Geometry::Rectangle::BoundingBox(coordinates);
        Geometry::GeoMetricAproximation metric = Geometry::GeoMetricAproximation::ComputeCorrection(boundingBox.center());
        CoordinateTree<Geometry::GeoMetricAproximation> ct(metric, coordinates);
        for (const Geometry::Point& coordinate : map[operatorName]) {
            map[operatorNames[ct.getNearestNeighbor(coordinate)]].emplace_back(coordinate);
        }
        map.erase(operatorName);
    }

    inline void reassignOutlierStations(const double minimalOutlierDistanceinCM = 2000000, const std::vector<std::string>& exceptions = std::vector<std::string>()) noexcept {
        std::vector<std::string> operatorNames;
        std::vector<Geometry::Point> coordinates;
        for (auto [key, value] : map) {
            if (Vector::contains(exceptions, key)) continue;
            for (const Geometry::Point& coordinate : value) {
                operatorNames.emplace_back(key);
                coordinates.emplace_back(coordinate);
            }
        }
        Geometry::Rectangle boundingBox = Geometry::Rectangle::BoundingBox(coordinates);
        Geometry::GeoMetricAproximation metric = Geometry::GeoMetricAproximation::ComputeCorrection(boundingBox.center());
        CoordinateTree<Geometry::GeoMetricAproximation> ct(metric, coordinates);
        Map<std::string, std::vector<Geometry::Point>> newMap;
        for (auto [key, value] : map) {
            Geometry::Rectangle operatorBoundingBox = Geometry::Rectangle::BoundingBox(value);
            Geometry::GeoMetricAproximation operatorMetric = Geometry::GeoMetricAproximation::ComputeCorrection(operatorBoundingBox.center());
            CoordinateTree<Geometry::GeoMetricAproximation> operatorCT(operatorMetric, value);
            for (const Geometry::Point& coordinate : value) {
                const size_t i = operatorCT.getNearestNeighbor(coordinate, 0.001);
                if (i >= value.size()) {
                    warning("No nearest neighbor with the same operator for the station ", coordinate, " was found!");
                    continue;
                }
                const double distance = geoDistanceInCM(coordinate, value[i]);
                if (distance > minimalOutlierDistanceinCM) {
                    const size_t j = ct.getNearestNeighbor(coordinate, 0.001);
                    if (j >= coordinates.size()) {
                        warning("No nearest neighbor for the station ", coordinate, " was found!");
                        continue;
                    }
                    const std::string newKey = operatorNames[j];
                    if (!newMap.contains(newKey)) {
                        newMap.insert(newKey, std::vector<Geometry::Point>());
                    }
                    newMap[newKey].emplace_back(coordinate);
                } else {
                    if (!newMap.contains(key)) {
                        newMap.insert(key, std::vector<Geometry::Point>());
                    }
                    newMap[key].emplace_back(coordinate);
                }
            }
        }
        map = std::move(newMap);
    }

    inline const Map<std::string, std::vector<Geometry::Point>>& getOperatorMap() const noexcept {
        return map;
    }

    inline Map<std::string, std::vector<Geometry::Point>> getOperatorMap() noexcept {
        return map;
    }

    inline std::vector<std::string> getOperators() const noexcept {
        std::vector<std::string> result;
        for (auto [key, value] : map) {
            suppressUnusedParameterWarning(value);
            result.emplace_back(key);
        }
        return result;
    }

    inline std::vector<Geometry::Point> getStations() const noexcept {
        std::vector<Geometry::Point> result;
        for (auto [key, value] : map) {
            suppressUnusedParameterWarning(key);
            for (const Geometry::Point& coordinate : value) {
                result.emplace_back(coordinate);
            }
        }
        return result;
    }

    inline std::vector<Geometry::Point> getStations(const std::string& operatorName) const noexcept {
        if (!map.contains(operatorName)) return std::vector<Geometry::Point>();
        std::vector<Geometry::Point> result;
        for (const Geometry::Point& coordinate : map.at(operatorName)) {
            result.emplace_back(coordinate);
        }
        return result;
    }

    inline std::vector<Geometry::Point> getStations(const std::vector<std::string>& exceptions) const noexcept {
        std::vector<Geometry::Point> result;
        for (auto [key, value] : map) {
            if (Vector::contains(exceptions, key)) continue;
            for (const Geometry::Point& coordinate : value) {
                result.emplace_back(coordinate);
            }
        }
        return result;
    }

    inline size_t numberOfOperators() const noexcept {
        return map.size();
    }

    inline size_t numberOfStations() const noexcept {
        size_t result = 0;
        for (auto [key, value] : map) {
            suppressUnusedParameterWarning(key);
            result += value.size();
        }
        return result;
    }

    inline void printInfo() const noexcept {
        std::cout << "Bike sharing station data:" << std::endl;
        std::cout << "   Number of Operators:      " << std::setw(12) << String::prettyInt(numberOfOperators()) << std::endl;
        std::cout << "   Number of Stations:       " << std::setw(12) << String::prettyInt(numberOfStations()) << std::endl;
    }

    inline void serialize(const std::string& fileName) const noexcept {
        IO::serialize(fileName, map);
    }

    inline void deserialize(const std::string& fileName) noexcept {
        IO::deserialize(fileName, map);
    }

    inline void readIntermediateFile(const std::string& fileName, const double coordinateFactor = 0.000001) noexcept {
        map.clear();
        std::ifstream bssFile(fileName);
        Assert(bssFile.is_open(), "Could not open file " << fileName << "!");
        std::string line = "";
        while (!bssFile.eof()) {
            getline(bssFile, line);
            std::vector<std::string> tokens = String::split(line, ' ');
            if (tokens.size() < 3) continue;
            if (!map.contains(tokens[2])) {
                map.insert(tokens[2], std::vector<Geometry::Point>());
            }
            map[tokens[2]].emplace_back(Construct::XY, String::lexicalCast<double>(tokens[0]) * coordinateFactor, String::lexicalCast<double>(tokens[1]) * coordinateFactor);
        }
        bssFile.close();
    }

    template<typename DOCUMENT>
    inline void drawTogether(DOCUMENT& doc, const std::string& country) const noexcept {
        doc.drawBordes(country);
        int i = 0;
        for (auto [key, value] : map) {
            doc.drawPoints(value, cyclicColor(i), cyclicIcon(i / 8), 30, key);
            i++;
        }
    }

    template<typename DOCUMENT>
    inline void drawSeparate(DOCUMENT& doc, const std::string& country) const noexcept {
        int i = 0;
        for (auto [key, value] : map) {
            Geometry::Rectangle operatorBoundingBox = Geometry::Rectangle::BoundingBox(value);
            Geometry::GeoMetricAproximation operatorMetric = Geometry::GeoMetricAproximation::ComputeCorrection(operatorBoundingBox.center());
            CoordinateTree<Geometry::GeoMetricAproximation> operatorCT(operatorMetric, value);
            double maxDistance = 0;
            Geometry::Point maxPoint;
            for (const Geometry::Point& coordinate : value) {
                const size_t j = operatorCT.getNearestNeighbor(coordinate, 0.001);
                if (j >= value.size()) continue;
                const double distance = geoDistanceInCM(coordinate, value[j]);
                if (maxDistance >= distance) continue;
                maxDistance = distance;
                maxPoint = coordinate;
            }
            doc.newPage();
            doc.drawBordes(country);
            doc.drawPoint(maxPoint, cyclicColor(i), Icon::MalteseCross, 80);
            doc.drawPoints(value, cyclicColor(i), cyclicIcon(i / 8 + 1), 30, key);
            i++;
        }
    }

private:
    Map<std::string, std::vector<Geometry::Point>> map;

};

}
