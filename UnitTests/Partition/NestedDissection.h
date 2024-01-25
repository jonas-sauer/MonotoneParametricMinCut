#pragma once

#include "../UnitTests.h"

#include "../../DataStructures/Partition/NestedDissection.h"
#include "../../DataStructures/Container/Set.h"

namespace UnitTests {

class NestedDissection {

public:
    inline void check() {
        ::NestedDissection nd = buildNestedDissection();
        UnitTests::check(nd.numberOfVertices() == 30, "Number of vertices should be 30 but is ", nd.numberOfVertices());
        UnitTests::check(nd.numberOfCells() == 31, "Number of cells should be 31 but is ", nd.numberOfCells());
        UnitTests::check(nd.numberOfLevels() == 5, "Number of levels should be 5 but is ", nd.numberOfLevels());

        UnitTests::check(nd.levelOfCell(0) == 0, "Level of cell 0 should be 0 but is ", nd.levelOfCell(0));
        UnitTests::check(nd.levelOfCell(1) == 1, "Level of cell 1 should be 1 but is ", nd.levelOfCell(1));
        UnitTests::check(nd.levelOfCell(2) == 1, "Level of cell 2 should be 1 but is ", nd.levelOfCell(2));
        UnitTests::check(nd.levelOfCell(3) == 2, "Level of cell 3 should be 2 but is ", nd.levelOfCell(3));
        UnitTests::check(nd.levelOfCell(4) == 2, "Level of cell 4 should be 2 but is ", nd.levelOfCell(4));
        UnitTests::check(nd.levelOfCell(5) == 2, "Level of cell 5 should be 2 but is ", nd.levelOfCell(5));
        UnitTests::check(nd.levelOfCell(6) == 2, "Level of cell 6 should be 2 but is ", nd.levelOfCell(6));
        UnitTests::check(nd.levelOfCell(7) == 3, "Level of cell 7 should be 3 but is ", nd.levelOfCell(7));
        UnitTests::check(nd.levelOfCell(8) == 3, "Level of cell 8 should be 3 but is ", nd.levelOfCell(8));
        UnitTests::check(nd.levelOfCell(9) == 3, "Level of cell 9 should be 3 but is ", nd.levelOfCell(9));

        UnitTests::check(nd.parentOfCell(0) == size_t(-1), "Parent of cell 0 should be -1 but is ", nd.parentOfCell(0));
        UnitTests::check(nd.parentOfCell(1) == 0, "Parent of cell 1 should be 0 but is ", nd.parentOfCell(1));
        UnitTests::check(nd.parentOfCell(2) == 0, "Parent of cell 2 should be 0 but is ", nd.parentOfCell(2));
        UnitTests::check(nd.parentOfCell(3) == 1, "Parent of cell 3 should be 1 but is ", nd.parentOfCell(3));
        UnitTests::check(nd.parentOfCell(4) == 1, "Parent of cell 4 should be 1 but is ", nd.parentOfCell(4));
        UnitTests::check(nd.parentOfCell(5) == 2, "Parent of cell 5 should be 2 but is ", nd.parentOfCell(5));
        UnitTests::check(nd.parentOfCell(6) == 2, "Parent of cell 6 should be 2 but is ", nd.parentOfCell(6));
        UnitTests::check(nd.parentOfCell(7) == 3, "Parent of cell 7 should be 3 but is ", nd.parentOfCell(7));
        UnitTests::check(nd.parentOfCell(8) == 3, "Parent of cell 8 should be 3 but is ", nd.parentOfCell(8));
        UnitTests::check(nd.parentOfCell(9) == 4, "Parent of cell 9 should be 4 but is ", nd.parentOfCell(9));

        UnitTests::check(nd.firstChildOfCell(0) == 1, "Left child of cell 0 should be 1 but is ", nd.firstChildOfCell(0));
        UnitTests::check(nd.firstChildOfCell(1) == 3, "Left child of cell 1 should be 3 but is ", nd.firstChildOfCell(1));
        UnitTests::check(nd.firstChildOfCell(2) == 5, "Left child of cell 2 should be 5 but is ", nd.firstChildOfCell(2));
        UnitTests::check(nd.firstChildOfCell(3) == 7, "Left child of cell 3 should be 7 but is ", nd.firstChildOfCell(3));
        UnitTests::check(nd.firstChildOfCell(4) == 9, "Left child of cell 4 should be 9 but is ", nd.firstChildOfCell(4));
        UnitTests::check(nd.firstChildOfCell(5) == 11, "Left child of cell 5 should be 11 but is ", nd.firstChildOfCell(5));
        UnitTests::check(nd.firstChildOfCell(6) == 13, "Left child of cell 6 should be 13 but is ", nd.firstChildOfCell(6));
        UnitTests::check(nd.firstChildOfCell(7) == 15, "Left child of cell 7 should be 15 but is ", nd.firstChildOfCell(7));
        UnitTests::check(nd.firstChildOfCell(8) == 17, "Left child of cell 8 should be 17 but is ", nd.firstChildOfCell(8));
        UnitTests::check(nd.firstChildOfCell(9) == 19, "Left child of cell 9 should be 19 but is ", nd.firstChildOfCell(9));

        UnitTests::check(nd.secondChildOfCell(0) == 2, "Right child of cell 0 should be 2 but is ", nd.secondChildOfCell(0));
        UnitTests::check(nd.secondChildOfCell(1) == 4, "Right child of cell 1 should be 4 but is ", nd.secondChildOfCell(1));
        UnitTests::check(nd.secondChildOfCell(2) == 6, "Right child of cell 2 should be 6 but is ", nd.secondChildOfCell(2));
        UnitTests::check(nd.secondChildOfCell(3) == 8, "Right child of cell 3 should be 8 but is ", nd.secondChildOfCell(3));
        UnitTests::check(nd.secondChildOfCell(4) == 10, "Right child of cell 4 should be 10 but is ", nd.secondChildOfCell(4));
        UnitTests::check(nd.secondChildOfCell(5) == 12, "Right child of cell 5 should be 12 but is ", nd.secondChildOfCell(5));
        UnitTests::check(nd.secondChildOfCell(6) == 14, "Right child of cell 6 should be 14 but is ", nd.secondChildOfCell(6));
        UnitTests::check(nd.secondChildOfCell(7) == 16, "Right child of cell 7 should be 16 but is ", nd.secondChildOfCell(7));
        UnitTests::check(nd.secondChildOfCell(8) == 18, "Right child of cell 8 should be 18 but is ", nd.secondChildOfCell(8));
        UnitTests::check(nd.secondChildOfCell(9) == 20, "Right child of cell 9 should be 20 but is ", nd.secondChildOfCell(9));

        UnitTests::check(nd.firstCellOfLevel(0) == 0, "First Cell of level 0 should be 0 but is ", nd.firstCellOfLevel(0));
        UnitTests::check(nd.firstCellOfLevel(1) == 1, "First Cell of level 1 should be 1 but is ", nd.firstCellOfLevel(1));
        UnitTests::check(nd.firstCellOfLevel(2) == 3, "First Cell of level 2 should be 3 but is ", nd.firstCellOfLevel(2));
        UnitTests::check(nd.firstCellOfLevel(3) == 7, "First Cell of level 3 should be 7 but is ", nd.firstCellOfLevel(3));
        UnitTests::check(nd.firstCellOfLevel(4) == 15, "First Cell of level 4 should be 15 but is ", nd.firstCellOfLevel(4));

        UnitTests::check(nd.lastCellOfLevel(0) == 0, "First Cell of level 0 should be 0 but is ", nd.firstCellOfLevel(0));
        UnitTests::check(nd.lastCellOfLevel(1) == 2, "First Cell of level 1 should be 1 but is ", nd.firstCellOfLevel(1));
        UnitTests::check(nd.lastCellOfLevel(2) == 6, "First Cell of level 2 should be 3 but is ", nd.firstCellOfLevel(2));
        UnitTests::check(nd.lastCellOfLevel(3) == 14, "First Cell of level 3 should be 7 but is ", nd.firstCellOfLevel(3));
        UnitTests::check(nd.lastCellOfLevel(4) == 30, "First Cell of level 4 should be 15 but is ", nd.firstCellOfLevel(4));

        UnitTests::check(nd.sizeOfCell(0) == 30, "Level of cell 0 should be 30 but is ", nd.sizeOfCell(0));
        UnitTests::check(nd.sizeOfCell(1) == 12, "Level of cell 1 should be 12 but is ", nd.sizeOfCell(1));
        UnitTests::check(nd.sizeOfCell(2) == 12, "Level of cell 2 should be 12 but is ", nd.sizeOfCell(2));
        UnitTests::check(nd.sizeOfCell(3) == 4, "Level of cell 3 should be 4 but is ", nd.sizeOfCell(3));
        UnitTests::check(nd.sizeOfCell(4) == 6, "Level of cell 4 should be 6 but is ", nd.sizeOfCell(4));
        UnitTests::check(nd.sizeOfCell(5) == 6, "Level of cell 5 should be 6 but is ", nd.sizeOfCell(5));
        UnitTests::check(nd.sizeOfCell(6) == 5, "Level of cell 6 should be 5 but is ", nd.sizeOfCell(6));
        UnitTests::check(nd.sizeOfCell(7) == 1, "Level of cell 7 should be 1 but is ", nd.sizeOfCell(7));
        UnitTests::check(nd.sizeOfCell(8) == 1, "Level of cell 8 should be 1 but is ", nd.sizeOfCell(8));
        UnitTests::check(nd.sizeOfCell(9) == 2, "Level of cell 9 should be 2 but is ", nd.sizeOfCell(9));

        IndexedSet<true> set;
        for (const Vertex vertex : nd.getCell(4)) {
            set.insert(vertex);
        }
        UnitTests::check(set.size() == 6, "Cell 4 should have size 2 but has size ", set.size());
        UnitTests::check(set.contains(15), "Cell 4 should contain vertex 15 but does not");
        UnitTests::check(set.contains(16), "Cell 4 should contain vertex 16 but does not");
        UnitTests::check(set.contains(20), "Cell 4 should contain vertex 20 but does not");
        UnitTests::check(set.contains(21), "Cell 4 should contain vertex 21 but does not");
        UnitTests::check(set.contains(25), "Cell 4 should contain vertex 25 but does not");
        UnitTests::check(set.contains(26), "Cell 4 should contain vertex 26 but does not");

        set.clear();
        for (const Vertex vertex : nd.getSeparatorOfCell(1)) {
            set.insert(vertex);
        }
        UnitTests::check(set.size() == 2, "Separator of cell 1 should have size 2 but has size ", set.size());
        UnitTests::check(set.contains(10), "Separator of cell 1 should contain vertex 10 but does not");
        UnitTests::check(set.contains(11), "Separator of cell 1 should contain vertex 11 but does not");

        set.clear();
        for (const Vertex vertex : nd.getSeparatorOfLevel(2)) {
            set.insert(vertex);
        }
        UnitTests::check(set.size() == 9, "Separator of level 2 should have size 2 but has size ", set.size());
        UnitTests::check(set.contains(2), "Separator of level 2 should contain vertex 2 but does not");
        UnitTests::check(set.contains(7), "Separator of level 2 should contain vertex 7 but does not");
        UnitTests::check(set.contains(12), "Separator of level 2 should contain vertex 12 but does not");
        UnitTests::check(set.contains(17), "Separator of level 2 should contain vertex 17 but does not");
        UnitTests::check(set.contains(22), "Separator of level 2 should contain vertex 22 but does not");
        UnitTests::check(set.contains(27), "Separator of level 2 should contain vertex 27 but does not");
        UnitTests::check(set.contains(10), "Separator of level 2 should contain vertex 10 but does not");
        UnitTests::check(set.contains(11), "Separator of level 2 should contain vertex 11 but does not");
        UnitTests::check(set.contains(18), "Separator of level 2 should contain vertex 18 but does not");

        UnitTests::check(nd.getCellIdOfVertex(Vertex(0)) == 3, "Cell of vertex 0 should be 3 but is ", nd.getCellIdOfVertex(Vertex(0)));
        UnitTests::check(nd.getCellIdOfVertex(Vertex(1)) == 7, "Cell of vertex 1 should be 7 but is ", nd.getCellIdOfVertex(Vertex(1)));
        UnitTests::check(nd.getCellIdOfVertex(Vertex(2)) == 0, "Cell of vertex 2 should be 0 but is ", nd.getCellIdOfVertex(Vertex(2)));
        UnitTests::check(nd.getCellIdOfVertex(Vertex(3)) == 11, "Cell of vertex 3 should be 11 but is ", nd.getCellIdOfVertex(Vertex(3)));
        UnitTests::check(nd.getCellIdOfVertex(Vertex(4)) == 24, "Cell of vertex 4 should be 24 but is ", nd.getCellIdOfVertex(Vertex(4)));
        UnitTests::check(nd.getCellIdOfVertex(Vertex(5)) == 8, "Cell of vertex 5 should be 8 but is ", nd.getCellIdOfVertex(Vertex(5)));
        UnitTests::check(nd.getCellIdOfVertex(Vertex(6)) == 3, "Cell of vertex 6 should be 3 but is ", nd.getCellIdOfVertex(Vertex(6)));
        UnitTests::check(nd.getCellIdOfVertex(Vertex(7)) == 0, "Cell of vertex 7 should be 0 but is ", nd.getCellIdOfVertex(Vertex(7)));
        UnitTests::check(nd.getCellIdOfVertex(Vertex(8)) == 23, "Cell of vertex 8 should be 23 but is ", nd.getCellIdOfVertex(Vertex(8)));
        UnitTests::check(nd.getCellIdOfVertex(Vertex(9)) == 11, "Cell of vertex 9 should be 11 but is ", nd.getCellIdOfVertex(Vertex(9)));

        UnitTests::check(nd.getCellIdOfVertex(Vertex(0), 2) == 3, "Cell of vertex 0 should be 3 but is ", nd.getCellIdOfVertex(Vertex(0), 2));
        UnitTests::check(nd.getCellIdOfVertex(Vertex(1), 2) == 3, "Cell of vertex 1 should be 3 but is ", nd.getCellIdOfVertex(Vertex(1), 2));
        UnitTests::check(nd.getCellIdOfVertex(Vertex(2), 2) == 0, "Cell of vertex 2 should be 0 but is ", nd.getCellIdOfVertex(Vertex(2), 2));
        UnitTests::check(nd.getCellIdOfVertex(Vertex(3), 2) == 5, "Cell of vertex 3 should be 5 but is ", nd.getCellIdOfVertex(Vertex(3), 2));
        UnitTests::check(nd.getCellIdOfVertex(Vertex(4), 2) == 5, "Cell of vertex 4 should be 5 but is ", nd.getCellIdOfVertex(Vertex(4), 2));
        UnitTests::check(nd.getCellIdOfVertex(Vertex(5), 2) == 3, "Cell of vertex 5 should be 3 but is ", nd.getCellIdOfVertex(Vertex(5), 2));
        UnitTests::check(nd.getCellIdOfVertex(Vertex(6), 2) == 3, "Cell of vertex 6 should be 3 but is ", nd.getCellIdOfVertex(Vertex(6), 2));
        UnitTests::check(nd.getCellIdOfVertex(Vertex(7), 2) == 0, "Cell of vertex 7 should be 0 but is ", nd.getCellIdOfVertex(Vertex(7), 2));
        UnitTests::check(nd.getCellIdOfVertex(Vertex(8), 2) == 5, "Cell of vertex 8 should be 5 but is ", nd.getCellIdOfVertex(Vertex(8), 2));
        UnitTests::check(nd.getCellIdOfVertex(Vertex(9), 2) == 5, "Cell of vertex 9 should be 5 but is ", nd.getCellIdOfVertex(Vertex(9), 2));
    }

protected:
    inline ::NestedDissection buildNestedDissection() const noexcept {
        ::NestedDissection nd(30);

        UnitTests::check(nd.numberOfCells() == 1, "Number of cells should be 1 but is ", nd.numberOfCells());
        UnitTests::check(nd.numberOfLevels() == 1, "Number of levels should be 1 but is ", nd.numberOfLevels());
        nd.divideCell(0, {Vertex(2), Vertex(7), Vertex(12), Vertex(17), Vertex(22), Vertex(27)}, {Vertex(0), Vertex(1), Vertex(5), Vertex(6), Vertex(10), Vertex(11), Vertex(15), Vertex(16), Vertex(20), Vertex(21), Vertex(25), Vertex(26)}, {Vertex(3), Vertex(4), Vertex(8), Vertex(9), Vertex(13), Vertex(14), Vertex(18), Vertex(19), Vertex(23), Vertex(24), Vertex(28), Vertex(29)});

        UnitTests::check(nd.numberOfCells() == 3, "Number of cells should be 3 but is ", nd.numberOfCells());
        UnitTests::check(nd.numberOfLevels() == 2, "Number of levels should be 2 but is ", nd.numberOfLevels());
        nd.divideCell(1, {Vertex(10), Vertex(11)}, {Vertex(0), Vertex(1), Vertex(5), Vertex(6)}, {Vertex(15), Vertex(16), Vertex(20), Vertex(21), Vertex(25), Vertex(26)});
        UnitTests::check(nd.numberOfCells() == 7, "Number of cells should be 7 but is ", nd.numberOfCells());
        UnitTests::check(nd.numberOfLevels() == 3, "Number of levels should be 3 but is ", nd.numberOfLevels());
        nd.divideCell(2, {Vertex(18)}, {Vertex(3), Vertex(4), Vertex(8), Vertex(9), Vertex(13), Vertex(14)}, {Vertex(19), Vertex(23), Vertex(24), Vertex(28), Vertex(29)});

        UnitTests::check(nd.numberOfCells() == 7, "Number of cells should be 7 but is ", nd.numberOfCells());
        UnitTests::check(nd.numberOfLevels() == 3, "Number of levels should be 3 but is ", nd.numberOfLevels());
        nd.divideCell(3, {Vertex(0), Vertex(6)}, {Vertex(1)}, {Vertex(5)});
        UnitTests::check(nd.numberOfCells() == 15, "Number of cells should be 15 but is ", nd.numberOfCells());
        UnitTests::check(nd.numberOfLevels() == 4, "Number of levels should be 4 but is ", nd.numberOfLevels());
        nd.divideCell(4, {Vertex(20), Vertex(21)}, {Vertex(25), Vertex(26)}, {Vertex(16), Vertex(15)});
        UnitTests::check(nd.numberOfCells() == 15, "Number of cells should be 15 but is ", nd.numberOfCells());
        UnitTests::check(nd.numberOfLevels() == 4, "Number of levels should be 4 but is ", nd.numberOfLevels());
        nd.divideCell(5, {Vertex(13)}, {Vertex(3), Vertex(9), Vertex(4), Vertex(8)}, {Vertex(14)});
        UnitTests::check(nd.numberOfCells() == 15, "Number of cells should be 15 but is ", nd.numberOfCells());
        UnitTests::check(nd.numberOfLevels() == 4, "Number of levels should be 4 but is ", nd.numberOfLevels());
        nd.divideCell(6, {Vertex(24)}, {Vertex(19)}, {Vertex(23), Vertex(28), Vertex(29)});

        UnitTests::check(nd.numberOfCells() == 15, "Number of cells should be 15 but is ", nd.numberOfCells());
        UnitTests::check(nd.numberOfLevels() == 4, "Number of levels should be 4 but is ", nd.numberOfLevels());
        nd.divideCell(11, {Vertex(3), Vertex(9)}, {Vertex(8)}, {Vertex(4)});
        UnitTests::check(nd.numberOfCells() == 31, "Number of cells should be 31 but is ", nd.numberOfCells());
        UnitTests::check(nd.numberOfLevels() == 5, "Number of levels should be 5 but is ", nd.numberOfLevels());
        nd.divideCell(14, {Vertex(28)}, {Vertex(23)}, {Vertex(29)});

        UnitTests::check(nd.numberOfCells() == 31, "Number of cells should be 31 but is ", nd.numberOfCells());
        UnitTests::check(nd.numberOfLevels() == 5, "Number of levels should be 5 but is ", nd.numberOfLevels());

        return nd;
    }

};

}
