#pragma once

#include <iostream>
#include <vector>

#include "../../Helpers/Assert.h"

template<typename ELEMENT_TYPE>
struct Less {
    using ElementType = ELEMENT_TYPE;
    inline bool operator()(const ElementType& a, const ElementType& b) const {return a < b;}
};

template<typename ELEMENT_TYPE, typename LESS_TYPE = Less<ELEMENT_TYPE>>
class Heap {

public:
    using ElementType = ELEMENT_TYPE;
    using LessType = LESS_TYPE;
    using Type = Heap<ElementType, LessType>;

public:
    Heap(const LessType& less) : less(less) {}
    Heap(const size_t size = 0, const LessType& less = LessType()) :
        less(less) {
        elements.reserve(size);
    }

    inline void push_back(const ElementType& element) {
        elements.emplace_back(element);
        siftUp(size() - 1);
    }

    template<typename... ARGS>
    inline void emplace_back(ARGS&&... args) {
        elements.emplace_back(args...);
        siftUp(size() - 1);
    }

    inline void remove_min() {
        Assert(!empty());
        siftDownHole(0);
    }

    inline ElementType pop_min() {
        Assert(!empty());
        ElementType result = std::move(elements[0]);
        siftDownHole(0);
        return result;
    }

    inline void decreaseKey(const size_t i) {
        Assert(i >= 0);
        Assert(i < elements.size());
        Assert((left(i) >= elements.size()) || (!less(elements[left(i)], elements[i])));
        Assert((right(i) >= elements.size()) || (!less(elements[right(i)], elements[i])));
        siftUp(i);
    }

    template<typename Range>
    inline void build(const Range& range) {
        clear();
        for (const ElementType& element : range) {
            elements.emplace_back(element);
        }
        for (size_t i = parent(size() - 1); i < size(); i--) {
            siftDown(i);
        }
    }

    inline void clear() {elements.clear();}
    inline size_t size() const {return elements.size();}
    inline bool empty() const {return elements.empty();}
    inline ElementType& min() {return elements[0];}
    inline const ElementType& min() const {return elements[0];}
    inline ElementType& front() {return elements[0];}
    inline const ElementType& front() const {return elements[0];}
    inline ElementType& operator[](const size_t i) {return elements[i];}
    inline const ElementType& operator[](const size_t i) const {return elements[i];}

    inline void printErrors(const size_t i = 0) {
        const size_t l = left(i);
        if ((l < size()) && (less(elements[l], elements[i]))) {
            std::cout << "Heap is broken! (" << i << ", " << l << ")" << std::endl;
            printErrors(l);
        }
        const size_t r = right(i);
        if ((r < size()) && (less(elements[r], elements[i]))) {
            std::cout << "Heap is broken! (" << i << ", " << r << ")" << std::endl;
            printErrors(r);
        }
    }

private:
    inline size_t left(const size_t i) const {return (i * 2) + 1;}
    inline size_t right(const size_t i) const {return (i * 2) + 2;}
    inline size_t parent(const size_t i) const {return (i - 1) / 2;}
    inline bool isLeaf(const size_t i) const {return left(i) >= size();}

    inline void siftDown(size_t i = 0) {
        using std::swap;
        Assert(i >= 0);
        Assert(i < size());
        while (true) {
            size_t minChild = left(i);
            if (minChild >= size()) return;
            if (minChild + 1 < size() && less(elements[minChild + 1], elements[minChild])) minChild++;
            if (less(elements[i], elements[minChild])) return;
            swap(elements[minChild], elements[i]);
            i = minChild;
        }
    }

    inline void siftDownHole(size_t i = 0) {
        AssertMsg(i >= 0, "Cannot siftDownHole " << i << " on heap of size " << size());
        AssertMsg(i < size(), "Cannot siftDownHole " << i << " on heap of size " << size());
        while (true) {
            const size_t r = right(i);
            if (r < size()) {
                const size_t l = left(i);
                if (less(elements[l], elements[r])) {
                    elements[i] = std::move(elements[l]);
                    i = l;
                } else {
                    elements[i] = std::move(elements[r]);
                    i = r;
                }
            } else {
                if (i + 1 < size()) {
                    elements[i] = std::move(elements.back());
                    elements.pop_back();
                    siftUp(i);
                    return;
                } else {
                    elements.pop_back();
                    return;
                }
            }
        }
    }

    inline void siftUp(size_t i) {
        using std::swap;
        Assert(i >= 0);
        Assert(i < size());
        while (i > 0) {
            const size_t p = parent(i);
            if (less(elements[i], elements[p])) {
                swap(elements[p], elements[i]);
                i = p;
            } else {
                return;
            }
        }
    }

private:
    std::vector<ElementType> elements;
    LessType less;

};
