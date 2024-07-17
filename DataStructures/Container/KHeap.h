#pragma once

#include <limits>
#include <iostream>
#include <cstring>
#include <cassert>
#include <vector>

#include "../../Helpers/Assert.h"

template<int LOG_K, typename KEY_TYPE>
class KHeap {

public:
    using KeyType = KEY_TYPE;
    static const int NULLINDEX = 0xFFFFFFFF;
    static const int K = 1 << LOG_K;

    struct HeapElement {
        KeyType key;
        int element;
        HeapElement() : key(0), element(0) {}
        HeapElement(const KeyType k, const int e) : key(k), element(e) {}
    };

public:
    KHeap(int n) : maxNumberOfElements(n), heap(NULL), position(NULL) {
        position = new int[n];
        heap = new HeapElement[n];
        initialize();
    }

    ~KHeap() {
        if (position != NULL) delete[] position;
        if (heap != NULL) delete[] heap;
    }

    inline int size() const noexcept {return numberOfElements;}
    inline bool empty() const noexcept {return size() == 0;}

    inline void extractMin(int &element, KeyType &key) noexcept {
        Assert(!empty(), "Trying to extract element from an empty heap!");
        HeapElement &front = heap[0];
        element = front.element;
        key = front.key;
        position[element] = NULLINDEX;
        numberOfElements--;
        if (!empty()) {
            front = heap[numberOfElements];
            position[front.element] = 0;
            siftDown(0);
        }
    }

    inline void popMin() noexcept {
        Assert(!empty(), "Trying to extract element from an empty heap!");
        HeapElement &front = heap[0];
        int element = front.element;
        position[element] = NULLINDEX;
        numberOfElements--;
        if (!empty()) {
            front = heap[numberOfElements];
            position[front.element] = 0;
            siftDown(0);
        }
    }

    inline void deleteElement(const int element) noexcept {
        if (!contains(element)) return;
        int elementPosition = position[element];
        HeapElement &heapElement = heap[elementPosition];
        KeyType previousKey = heapElement.key;
        position[element] = NULLINDEX;
        numberOfElements--;
        if (elementPosition < numberOfElements) {
            heapElement = heap[numberOfElements];
            position[heapElement.element] = elementPosition;
            if (heapElement.key > previousKey) {
                siftDown(elementPosition);
            } else {
                siftUp(elementPosition);
            }
        }
    }

    inline void update(const int element, const KeyType key) noexcept {
        if (position[element] == NULLINDEX) {
            HeapElement &back = heap[numberOfElements];
            back.key = key;
            back.element = element;
            position[element] = numberOfElements;
            siftUp(numberOfElements++);
        } else {
            int el_pos = position[element];
            HeapElement &el = heap[el_pos];
            if (key > el.key) {
                el.key = key;
                siftDown(el_pos);
            } else {
                el.key = key;
                siftUp(el_pos);
            }
        }
    }

    inline KeyType minKey() const noexcept {
        Assert(!empty(), "Heap is empty!");
        return heap[0].key;
    }

    inline int minElement() const noexcept {
        Assert(!empty(), "Heap is empty!");
        return heap[0].element;
    }

    inline void min(int &element, KeyType &key) noexcept {
        Assert(!empty(), "Heap is empty!");
        key = heap[0].key;
        element = heap[0].element;
    }

    inline KeyType key(const int element) const noexcept {
        Assert(contains(element), "Element " << element << " is not contained!");
        return heap[position[element]].key;
    }

    inline void reset() noexcept {initialize();}

    inline void clear() noexcept {
        for (int i = 0; i < numberOfElements; ++i) {
            position[heap[i].element] = NULLINDEX;
        }
        numberOfElements = 0;
    }

    inline bool contains(const int element) const noexcept {
        return position[element] != NULLINDEX;
    }

    void dump(std::ostream& os) const noexcept {
        os << "Heap: ";
        for (int i = 0; i < size(); ++i) {
            std::cout << "[" << i << ":" << heap[i].key << "] " << std::flush;
        }
        os << std::endl;
        os << "Pos:  ";
        for (int e = 0; e < maxNumberOfElements; ++e) {
            if (position[e] != NULLINDEX) {
                os << "[" << e << ":" << position[e] << "] " << std::flush;
            }
        }
        os << std::endl;
    }

protected:
    inline void initialize() noexcept {
        numberOfElements = 0;
        memset(position, 0xFF, sizeof(int) * maxNumberOfElements);
        for (int i = 0; i < maxNumberOfElements; ++i) {
            assert(position[i] == NULLINDEX);
        }
    }

    inline void siftUp(int i) noexcept {
        assert(i < numberOfElements);
        int cur_i = i;
        while (cur_i > 0) {
            int parent_i = (cur_i-1) >> LOG_K;
            if (heap[parent_i].key > heap[cur_i].key) {
                swap(cur_i, parent_i);
            } else {
                break;
            }
            cur_i = parent_i;
        }
    }

    inline void siftDown(int i) noexcept {
        assert(i < numberOfElements);
        while (true) {
            int min_ind = i;
            KeyType min_key = heap[i].key;
            int child_ind_l = (i << LOG_K) + 1;
            int child_ind_u = std::min(child_ind_l + K, numberOfElements);
            for (int j = child_ind_l; j < child_ind_u; ++j) {
                if (heap[j].key < min_key) {
                    min_ind = j;
                    min_key = heap[j].key;
                }
            }
            if (min_ind != i) {
                swap(i, min_ind);
                i = min_ind;
            } else {
                break;
            }
        }
    }

    inline void swap(const int i, const int j) noexcept {
        HeapElement &el_i = heap[i];
        HeapElement &el_j = heap[j];
        position[el_i.element] = j;
        position[el_j.element] = i;
        HeapElement temp = el_i;
        el_i = el_j;
        el_j = temp;
    }

private:
    int numberOfElements;
    int maxNumberOfElements;
    HeapElement *heap;
    int *position;

};
