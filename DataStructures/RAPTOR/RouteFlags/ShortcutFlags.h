#pragma once

#include "../Data.h"
#include "../../../Helpers/Types.h"

#include <vector>

namespace RAPTOR::RouteFlags {

struct ShortcutFlags {
    ShortcutFlags(const size_t fromCells, const size_t toCells, const TransferGraph& shortcuts) :
        shortcuts(shortcuts),
        shortcutFlags(fromCells, std::vector<std::vector<bool>>(toCells, std::vector<bool>(shortcuts.numEdges(), false))) {
    }

    inline void markShortcut(const size_t fromCell, const std::vector<size_t>& toCells, const Vertex fromVertex, const Vertex toVertex) noexcept {
        const Edge edge = shortcuts.findEdge(fromVertex, toVertex);
        if (edge == noEdge) return;
        for (const size_t toCell : toCells) {
            markShortcut(fromCell, toCell, edge);
        }
    }

    inline void markShortcut(const std::vector<size_t>& fromCells, const std::vector<size_t>& toCells, const Vertex fromVertex, const Vertex toVertex) noexcept {
        const Edge edge = shortcuts.findEdge(fromVertex, toVertex);
        if (edge == noEdge) return;
        for (const size_t fromCell : fromCells) {
            for (const size_t toCell : toCells) {
                markShortcut(fromCell, toCell, edge);
            }
        }
    }

    inline void markShortcut(const std::vector<size_t>& fromCells, const size_t toCell, const Vertex fromVertex, const Vertex toVertex) noexcept {
        const Edge edge = shortcuts.findEdge(fromVertex, toVertex);
        if (edge == noEdge) return;
        for (const size_t fromCell : fromCells) {
            markShortcut(fromCell, toCell, edge);
        }
    }

    inline void markShortcut(const size_t fromCell, const size_t toCell, const Vertex fromVertex, const Vertex toVertex) noexcept {
        const Edge edge = shortcuts.findEdge(fromVertex, toVertex);
        if (edge == noEdge) return;
        markShortcut(fromCell, toCell, edge);

    }

    inline void markShortcut(const size_t fromCell, const size_t toCell, const Edge edge) noexcept {
        shortcutFlags[fromCell][toCell][edge] = true;
    }

    inline void incorporate(const ShortcutFlags& other) noexcept {
        for (size_t i = 0; i < shortcutFlags.size(); i++) {
            for (size_t j = 0; j < shortcutFlags[i].size(); j++) {
                for (size_t k = 0; k < shortcutFlags[i][j].size(); k++) {
                    shortcutFlags[i][j][k] = shortcutFlags[i][j][k] | other.shortcutFlags[i][j][k];
                }
            }
        }
    }

    const TransferGraph& shortcuts;
    std::vector<std::vector<std::vector<bool>>> shortcutFlags;
};

}
