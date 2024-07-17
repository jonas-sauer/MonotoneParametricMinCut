#pragma once

#include "../../../Helpers/Types.h"

namespace CSA {

struct ProfileEntry {
    ProfileEntry(const int departureTime = never, const int arrivalTime = never) :
        departureTime(departureTime),
        arrivalTime(arrivalTime) {
    }

    int departureTime;
    int arrivalTime;

    inline friend std::ostream& operator<<(std::ostream& out, const ProfileEntry& entry) noexcept {
        return out << "Departure time: " << entry.departureTime << ", Arrival time: " << entry.arrivalTime;
    }
};

}
