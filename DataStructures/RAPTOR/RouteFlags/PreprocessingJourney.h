#pragma once

#include "../Data.h"

namespace RAPTOR::RouteFlags {
    struct Shortcut {
        Vertex origin;
        Vertex destination;
        int travelTime;
    };

    struct PreprocessingJourney {
        std::vector<const StopEvent*> departureStopEvents;
        std::vector<const StopEvent*> arrivalStopEvents;
        std::vector<Shortcut> shortcuts;

        inline bool empty() const noexcept {
            return departureStopEvents.empty() && arrivalStopEvents.empty() && shortcuts.empty();
        }

        inline void reverse() noexcept {
            Vector::reverse(departureStopEvents);
            Vector::reverse(arrivalStopEvents);
            Vector::reverse(shortcuts);
        }

        inline PreprocessingJourney& operator+=(const PreprocessingJourney& otherJourney) noexcept {
            departureStopEvents += otherJourney.departureStopEvents;
            arrivalStopEvents += otherJourney.arrivalStopEvents;
            shortcuts += otherJourney.shortcuts;
            return *this;
        }
    };
}
