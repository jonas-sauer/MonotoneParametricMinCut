#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "Color.h"
#include "Visualization.h"

#include "../Helpers/Helpers.h"

template<typename DOCUMENT_TYPE>
class FunctionVisualization : public Visualization<DOCUMENT_TYPE> {

public:
    using DocumentType = DOCUMENT_TYPE;
    using Type = FunctionVisualization<DocumentType>;
    using Point = Geometry::Point;
    using Rectangle = Geometry::Rectangle;

private:
    using Super = Visualization<DocumentType>;
    constexpr static int Height = 2000;
    constexpr static int MarginLeft = 400;
    constexpr static int MarginRight = 50;
    constexpr static int MarginTop = 50;
    constexpr static int MarginBottom = 100;

private:
    inline static Rectangle extendedBoundingBox(const Rectangle& boundingBox, const int width) noexcept {
        const double ux = boundingBox.dx() / (width - MarginLeft - MarginRight);
        const double uy = boundingBox.dy() / (Height - MarginTop - MarginBottom);
        return Rectangle::BoundingBox(boundingBox.min - Point(Construct::XY, MarginLeft*ux, MarginBottom*uy), boundingBox.max + Point(Construct::XY, MarginRight*ux, MarginTop*uy));
    }

    FunctionVisualization(const std::string& filename, const int width, const int height, const Rectangle& boundingBox) :
        Super(filename, width, height, extendedBoundingBox(boundingBox, width)),
        boundingBox(boundingBox),
        ux(Super::ex(10)),
        uy(Super::ey(10)) {
        Super::setDefaultSize(32);
        Super::setDefaultStroke(5);
    }

public:
    inline static Type FromBoundingBox(const std::string& filename, const Rectangle& boundingBox, const double aspectRatio) noexcept {
        const int width = ((Height - MarginTop - MarginBottom) * aspectRatio) + MarginLeft + MarginRight;
        return Type(filename, width, Height, boundingBox);
    }

    inline static Type FromBoundingBox(const std::string& filename, const Rectangle& boundingBox) noexcept {
        const int width = ((Height - MarginTop - MarginBottom) * boundingBox.dx()) / (boundingBox.dy()) + MarginLeft + MarginRight;
        return Type(filename, width, Height, boundingBox);
    }

    inline static Type FromBoundingBox(const std::string& filename, const Rectangle& boundingBox, const Point& canvasSize) noexcept {
        return Type(filename, canvasSize.x + MarginLeft + MarginRight, canvasSize.y + MarginBottom + MarginTop, boundingBox);
    }

    inline static Type FromFunction(const std::string& filename, const std::vector<Point>& function, const double aspectRatio) noexcept {
        return FromBoundingBox(filename, Rectangle::BoundingBox(function), aspectRatio);
    }

    inline static Type FromFunction(const std::string& filename, const std::vector<Point>& function) noexcept {
        return FromBoundingBox(filename, Rectangle::BoundingBox(function));
    }

    inline static Type FromFunction(const std::string& filename, const std::vector<Point>& function, const double dx, const double dy, const double aspectRatio) noexcept {
        Rectangle boundingBox = Rectangle::BoundingBox(function);
        boundingBox.min.x = floor(boundingBox.min.x, dx);
        boundingBox.min.y = floor(boundingBox.min.y, dy);
        boundingBox.max.x = ceil(boundingBox.max.x, dx);
        boundingBox.max.y = ceil(boundingBox.max.y, dy);
        return FromBoundingBox(filename, boundingBox, aspectRatio);
    }

    inline static Type FromFunction(const std::string& filename, const std::vector<Point>& function, const double dx, const double dy) noexcept {
        Rectangle boundingBox = Rectangle::BoundingBox(function);
        boundingBox.min.x = floor(boundingBox.min.x, dx);
        boundingBox.min.y = floor(boundingBox.min.y, dy);
        boundingBox.max.x = ceil(boundingBox.max.x, dx);
        boundingBox.max.y = ceil(boundingBox.max.y, dy);
        return FromBoundingBox(filename, boundingBox);
    }

    inline static Type FromFunctionAndBoundingBox(const std::string& filename, const std::vector<Point>& function, const double dx, const double dy, const double aspectRatio) noexcept {
        Rectangle boundingBox = Rectangle::BoundingBox(function);
        boundingBox.min.x = floor(boundingBox.min.x, dx);
        boundingBox.min.y = floor(boundingBox.min.y, dy);
        boundingBox.max.x = ceil(boundingBox.max.x, dx);
        boundingBox.max.y = ceil(boundingBox.max.y, dy);
        return FromBoundingBox(filename, boundingBox, aspectRatio);
    }

    inline static Type FromFunctionAndBoundingBox(const std::string& filename, const std::vector<Point>& function, const double dx, const double dy) noexcept {
        Rectangle boundingBox = Rectangle::BoundingBox(function);
        boundingBox.min.x = floor(boundingBox.min.x, dx);
        boundingBox.min.y = floor(boundingBox.min.y, dy);
        boundingBox.max.x = ceil(boundingBox.max.x, dx);
        boundingBox.max.y = ceil(boundingBox.max.y, dy);
        return FromBoundingBox(filename, boundingBox);
    }

public:
    inline void clearOverflow() noexcept {
        Super::fillRectangle(Rectangle::BoundingBox(Super::boundingBox.min, Point(Construct::XY, boundingBox.min.x, Super::boundingBox.max.y)), Color::White);
        Super::fillRectangle(Rectangle::BoundingBox(Super::boundingBox.min, Point(Construct::XY, Super::boundingBox.min.x, boundingBox.max.y)), Color::White);
        Super::fillRectangle(Rectangle::BoundingBox(Point(Construct::XY, boundingBox.max.x, Super::boundingBox.min.y), Super::boundingBox.max), Color::White);
        Super::fillRectangle(Rectangle::BoundingBox(Point(Construct::XY, Super::boundingBox.max.x, boundingBox.min.y), Super::boundingBox.max), Color::White);
    }

    inline void drawBox() noexcept {
        Super::drawLine(boundingBox.min, Point(Construct::XY, boundingBox.min.x, boundingBox.max.y), 4);
        Super::drawLine(boundingBox.min, Point(Construct::XY, boundingBox.max.x, boundingBox.min.y), 4);
        Super::drawLine(boundingBox.max, Point(Construct::XY, boundingBox.min.x, boundingBox.max.y), 4);
        Super::drawLine(boundingBox.max, Point(Construct::XY, boundingBox.max.x, boundingBox.min.y), 4);
    }

    inline void drawGrid(const double dx, const double dy, const Color& color) noexcept {
        for (double x = boundingBox.min.x; x <= boundingBox.max.x; x += dx) {
            Super::drawLine(x, boundingBox.min.y, x, boundingBox.max.y, color, 1);
        }
        for (double y = boundingBox.min.y; y <= boundingBox.max.y; y += dy) {
            Super::drawLine(boundingBox.min.x, y, boundingBox.max.x, y, color, 1);
        }
    }
    inline void drawGrid(const double dx, const double dy) noexcept {drawGrid(dx, dy, Color::LightGrey);}

    template<typename LABEL>
    inline void drawXAxis(const double dx, const LABEL& label, const Color& color) noexcept {
        Super::drawLine(boundingBox.min, Point(Construct::XY, boundingBox.max.x, boundingBox.min.y), color, 6);
        for (double x = boundingBox.min.x; x <= boundingBox.max.x; x += dx) {
            Super::drawLine(x, boundingBox.min.y - uy, x, boundingBox.min.y + uy, color, 6);
            Super::write(label(x), x, boundingBox.min.y - 1.5*uy - Super::ey(), color, XAlignment::Center, YAlignment::Bottom);
        }
    }
    template<typename LABEL>
    inline void drawXAxis(const double dx, const LABEL& label) noexcept {drawXAxis(dx, label, Super::defaultColor);}
    inline void drawXAxis(const double dx) noexcept {drawXAxis(dx, [](const double x){return std::to_string(x);}, Super::defaultColor);}

    template<typename LABEL>
    inline void drawYAxis(const double dy, const LABEL& label, const Color& color) noexcept {
        Super::drawLine(boundingBox.min, Point(Construct::XY, boundingBox.min.x, boundingBox.max.y), color, 6);
        for (double y = boundingBox.min.y; y <= boundingBox.max.y; y += dy) {
            Super::drawLine(boundingBox.min.x - ux, y, boundingBox.min.x + ux, y, color, 6);
            Super::write(label(y), boundingBox.min.x - 2*ux, y - 0.4*Super::ey(), color, XAlignment::Right, YAlignment::Bottom);
        }
    }
    template<typename LABEL>
    inline void drawYAxis(const double dy, const LABEL& label) noexcept {drawYAxis(dy, label, Super::defaultColor);}
    inline void drawYAxis(const double dy) noexcept {drawYAxis(dy, [](const double y){return std::to_string(y);}, Super::defaultColor);}

    inline void drawFunction(const std::vector<Point>& function, const Color& color, const double stroke) noexcept {
        if (function.size() < 2) return;
        for (size_t i = 0; i < function.size() - 1; i++) {
            if (function[i].x == function[i + 1].x) {
                Super::setDashPattern({3 * stroke});
            } else {
                Super::clearDashPattern();
            }
            cairo_set_line_width(Super::context, stroke);
            cairo_set_source_rgba(Super::context, color.r, color.g, color.b, color.a);
            cairo_move_to(Super::context, Super::x(function[i]), Super::y(function[i]));
            cairo_line_to(Super::context, Super::x(function[i + 1]), Super::y(function[i + 1]));
            cairo_stroke(Super::context);
        }
        Super::clearDashPattern();
    }
    inline void drawFunction(const std::vector<Point>& function, const Color& color) noexcept {drawFunction(function, color, Super::defaultStroke);}
    inline void drawFunction(const std::vector<Point>& function, const double stroke) noexcept {drawFunction(function, Super::defaultColor, stroke);}
    inline void drawFunction(const std::vector<Point>& function) noexcept {drawFunction(function, Super::defaultColor, Super::defaultStroke);}

private:
    const Rectangle boundingBox;
    const double ux;
    const double uy;

};
