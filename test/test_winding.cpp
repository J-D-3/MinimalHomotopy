// Copyright Ingo Proff 2016-2026.
// https://github.com/CrikeeIP/MinimalHomotopy
// Distributed under the MIT Software License (X11 license).

#include <doctest/doctest.h>

#include "minimal_homotopy/winding.hpp"

using namespace mh;

namespace {
Point_2 P(double x, double y) { return Point_2(x, y); }
}  // namespace

TEST_CASE("unit square loop has winding area 1 and is consistent") {
  // CCW unit square, returned to its start.
  Curve square = {P(0, 0), P(1, 0), P(1, 1), P(0, 1), P(0, 0)};
  auto a = analyze_closed_curve(square);
  CHECK(a.consistent);
  CHECK(a.total_area == doctest::Approx(1.0));
}

TEST_CASE("figure-eight has inconsistent winding numbers") {
  // Two lobes traced in opposite senses: one face wn = +1, the other wn = -1.
  Curve eight = {P(0, 0), P(1, 1), P(2, 0), P(1, -1), P(0, 0),
                 P(-1, 1), P(-2, 0), P(-1, -1), P(0, 0)};
  auto a = analyze_closed_curve(eight);
  CHECK_FALSE(a.consistent);
}

TEST_CASE("a triangle loop reports its area") {
  Curve tri = {P(0, 0), P(4, 0), P(0, 3), P(0, 0)};  // legs 4 and 3 -> area 6
  auto a = analyze_closed_curve(tri);
  CHECK(a.consistent);
  CHECK(a.total_area == doctest::Approx(6.0));
}
