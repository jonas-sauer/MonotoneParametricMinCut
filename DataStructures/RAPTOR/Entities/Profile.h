#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "ArrivalLabel.h"

#include "../../../Helpers/Assert.h"
#include "../../../Helpers/Vector/Vector.h"
#include "../../../Helpers/String/String.h"

#include "../../Geometry/Point.h"

namespace RAPTOR {

struct ProfileEntry {

    ProfileEntry(const int departureTime = never, const int arrivalTime = never, const size_t numberOfTrips = -1) :
        departureTime(departureTime),
        arrivalTime(arrivalTime),
        numberOfTrips(numberOfTrips) {
    }

    ProfileEntry(const ArrivalLabel& other, const int departureTime) :
        departureTime(departureTime),
        arrivalTime(other.arrivalTime),
        numberOfTrips(other.numberOfTrips) {
    }

    inline int travelTime() const noexcept {
        return arrivalTime - departureTime;
    }

    inline bool operator<(const ProfileEntry& other) const noexcept {
        return (departureTime < other.departureTime) || ((departureTime == other.departureTime) && (arrivalTime < other.arrivalTime));
    }

    inline bool operator==(const ProfileEntry& other) const noexcept {
        return (departureTime == other.departureTime) && (arrivalTime == other.arrivalTime) && (numberOfTrips == other.numberOfTrips);
    }

    inline bool dominates(const ProfileEntry& other) const noexcept {
        return (departureTime >= other.departureTime) && (arrivalTime <= other.arrivalTime);
    }

    inline friend std::ostream& operator<<(std::ostream& out, const ProfileEntry& entry) noexcept {
        return out << "Departure: " << String::secToTime(entry.departureTime) << ", Arrival: " << String::secToTime(entry.arrivalTime) << ", trips: " << entry.numberOfTrips;
    }

    int departureTime;
    int arrivalTime;
    size_t numberOfTrips;
};

using Profile = std::vector<ProfileEntry>;

class ProfileHandle : public Profile {

private:
    struct Function {
        inline void emplace_back(const int x, const int y) noexcept {
            if (points.empty() || points.back().x != x || points.back().y < y) {
                points.emplace_back(Geometry::Point(Construct::XY, x, y));
            } else {
                points.back().y = y;
            }
        }
        inline const Geometry::Point& back() const noexcept {
            return points.back();
        }
        std::vector<Geometry::Point> points;
    };

public:
    ProfileHandle(const Profile& profile, const int minDepartureTime, const int maxDepartureTime, const int walkingTime = intMax) :
        Profile(profile),
        minDepartureTime(minDepartureTime),
        maxDepartureTime(maxDepartureTime),
        walkingTime(std::min(walkingTime, intMax / 2)) {
    }

    inline void swap(ProfileHandle& other) noexcept {
        Profile::swap(*this);
        std::swap(minDepartureTime, other.minDepartureTime);
        std::swap(maxDepartureTime, other.maxDepartureTime);
        std::swap(walkingTime, other.walkingTime);
    }

    inline void sort() noexcept {
        std::stable_sort(begin(), end());
    }

    inline int getMaxTrips() const noexcept {
        return getMaxTrips(*this);
    }

    inline std::vector<Geometry::Point> getArrivalTimeFunction() const noexcept {
        return constructArrivalTimeFunction(*this);
    }

    inline std::vector<Geometry::Point> getTravelTimeFunction() const noexcept {
        return constructTravelTimeFunction(*this);
    }

    inline std::vector<Geometry::Point> getWalkingArrivalTimeFunction() const noexcept {
        std::vector<Geometry::Point> result;
        if (walkingTime < intMax / 2) {
            result.emplace_back(Geometry::Point(Construct::XY, minDepartureTime, minDepartureTime + walkingTime));
            result.emplace_back(Geometry::Point(Construct::XY, maxDepartureTime, maxDepartureTime + walkingTime));
        }
        return result;
    }

    inline std::vector<Geometry::Point> getWalkingTravelTimeFunction() const noexcept {
        std::vector<Geometry::Point> result;
        if (walkingTime < intMax / 2) {
            result.emplace_back(Geometry::Point(Construct::XY, minDepartureTime, walkingTime));
            result.emplace_back(Geometry::Point(Construct::XY, maxDepartureTime, walkingTime));
        }
        return result;
    }

    inline std::vector<std::vector<Geometry::Point>> getArrivalTimeFunctions() const noexcept {
        std::vector<std::vector<Geometry::Point>> result{getWalkingArrivalTimeFunction()};
        Profile profile = reduceMaxTripsTo(getMaxTrips());
        while (!profile.empty()) {
            result.emplace_back(constructArrivalTimeFunction(profile));
            profile = reduceMaxTrips(profile);
        }
        return result;
    }

    inline std::vector<std::vector<Geometry::Point>> getTravelTimeFunctions() const noexcept {
        std::vector<std::vector<Geometry::Point>> result{getWalkingTravelTimeFunction()};
        Profile profile = reduceMaxTripsTo(getMaxTrips());
        while (!profile.empty()) {
            result.emplace_back(constructTravelTimeFunction(profile));
            profile = reduceMaxTrips(profile);
        }
        return result;
    }

    inline ProfileHandle reduceToEarliestArrivalTime() const noexcept {
        Profile profile;
        ProfileEntry bestAfterDepartureTime;
        size_t i = size() - 1;
        for (; i != size_t(-1) && (*this)[i].departureTime >= maxDepartureTime; i--) {
            const ProfileEntry& entry = (*this)[i];
            if (entry.arrivalTime < bestAfterDepartureTime.arrivalTime) {
                bestAfterDepartureTime = entry;
            }
        }
        if (bestAfterDepartureTime.arrivalTime < never) {
            profile.emplace_back(bestAfterDepartureTime);
        }
        for (; i != size_t(-1); i--) {
            const ProfileEntry& entry = (*this)[i];
            if (entry.arrivalTime - entry.departureTime >= walkingTime) continue;
            if (!profile.empty() && entry.arrivalTime >= profile.back().arrivalTime) continue;
            if (!profile.empty() && entry.departureTime == profile.back().departureTime) {
                profile.back() = entry;
            } else {
                profile.emplace_back(entry);
            }
        }
        Vector::reverse(profile);
        return ProfileHandle(profile, minDepartureTime, maxDepartureTime, walkingTime);
    }

    inline friend std::ostream& operator<<(std::ostream& out, const ProfileHandle& profile) noexcept {
        out << std::endl << "Profile: " << std::endl;
        out << "    Walking time = " << String::secToString(profile.walkingTime) << std::endl;
        if (profile.empty()) {
            out << "    -- no entries --" << std::endl;
        }
        for (const auto& entry : profile) {
            out << "    " << entry << std::endl;
        }
        return out;
    }

private:
    inline int getMaxTrips(const Profile& profile) const noexcept {
        size_t maxTrips = 0;
        for (const ProfileEntry& entry : profile) {
            if (maxTrips < entry.numberOfTrips) {
                maxTrips = entry.numberOfTrips;
            }
        }
        return maxTrips;
    }

    inline Profile reduceMaxTripsTo(const size_t maxTrips) const noexcept {
        Profile result;
        for (const ProfileEntry& entry : *this) {
            if (entry.numberOfTrips > maxTrips) continue;
            if ((!result.empty()) && (result.back().dominates(entry))) continue;
            while ((!result.empty()) && (entry.dominates(result.back()))) result.pop_back();
            result.emplace_back(entry);
        }
        return result;
    }

    inline Profile reduceMaxTrips(const Profile& profile) const noexcept {
        return reduceMaxTripsTo(getMaxTrips(profile) - 1);
    }

    inline std::vector<Geometry::Point> constructArrivalTimeFunction(const Profile& profile) const noexcept {
        AssertMsg(!profile.empty(), "Cannot construct Function from empty profile!");
        Function function;
        function.emplace_back(minDepartureTime, std::min(profile.front().arrivalTime, minDepartureTime + walkingTime));
        for (size_t i = 0; i < profile.size(); i++) {
            if (profile[i].arrivalTime - profile[i].departureTime >= walkingTime) continue;
            if (profile[i].arrivalTime - function.back().x >= walkingTime) {
                function.emplace_back(function.back().x, function.back().x + walkingTime);
                function.emplace_back(profile[i].arrivalTime - walkingTime, profile[i].arrivalTime);
            } else {
                function.emplace_back(function.back().x, profile[i].arrivalTime);
            }
            function.emplace_back(profile[i].departureTime, profile[i].arrivalTime);
        }
        if (function.back().x < maxDepartureTime) {
            function.emplace_back(function.back().x, function.back().x + walkingTime);
            function.emplace_back(maxDepartureTime, maxDepartureTime + walkingTime);
        }
        return function.points;
    }

    inline std::vector<Geometry::Point> constructTravelTimeFunction(const Profile& profile) const noexcept {
        AssertMsg(!profile.empty(), "Cannot construct Function from empty profile!");
        Function function;
        function.emplace_back(minDepartureTime, std::min(profile.front().arrivalTime - minDepartureTime, walkingTime));
        for (size_t i = 0; i < profile.size(); i++) {
            if (profile[i].arrivalTime - profile[i].departureTime >= walkingTime) continue;
            if (profile[i].arrivalTime - function.back().x >= walkingTime) {
                function.emplace_back(function.back().x, walkingTime);
                function.emplace_back(profile[i].arrivalTime - walkingTime, walkingTime);
            } else {
                function.emplace_back(function.back().x, profile[i].arrivalTime - function.back().x);
            }
            function.emplace_back(profile[i].departureTime, profile[i].arrivalTime - profile[i].departureTime);
        }
        if (function.back().x < maxDepartureTime) {
            function.emplace_back(function.back().x, walkingTime);
            function.emplace_back(maxDepartureTime, walkingTime);
        }
        return function.points;
    }

private:
    int minDepartureTime;
    int maxDepartureTime;
    int walkingTime;

};

}
