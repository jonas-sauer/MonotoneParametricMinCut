#pragma once

#include <iostream>
#include <vector>
#include <string>

using DestinationSpecificPassengerId = u_int32_t;
using GlobalPassengerId = u_int64_t;

using DestinationSpecificPassengerList = std::vector<DestinationSpecificPassengerId>;
using GlobalPassengerList = std::vector<GlobalPassengerId>;

inline constexpr DestinationSpecificPassengerId getDestinationSpecificPassengerId(const GlobalPassengerId id) noexcept {
    return static_cast<DestinationSpecificPassengerId>(id);
}

inline constexpr GlobalPassengerId getGlobalPassengerId(const u_int32_t destination, const DestinationSpecificPassengerId id) noexcept {
    return (static_cast<GlobalPassengerId>(destination) << 32) + (static_cast<u_int32_t>(id));
}

inline constexpr u_int32_t getDestination(const GlobalPassengerId id) noexcept {
    return static_cast<u_int32_t>(id >> 32);
}
