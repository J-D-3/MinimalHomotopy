// Copyright Ingo Proff 2016-2026.
// https://github.com/CrikeeIP/MinimalHomotopy
// Distributed under the MIT Software License (X11 license).
//
// Heavier cases: many anchors, non-convex consistent regions, and invariants
// (symmetry, translation, scaling) that must hold for any valid input.

#include <doctest/doctest.h>

#include "minimal_homotopy/homotopy_area.hpp"

#include <vector>

using namespace mh;

namespace {

Point_2 P(double x, double y) { return Point_2(x, y); }

double area_of(const Curve& p, const Curve& q) {
  auto r = calculate_min_homotopy_area(p, q);
  REQUIRE(r.ok());
  return *r;
}

// Canonical zig-zag: P is the segment (0,0)->(N,0); Q rises/falls to unit peaks
// at every integer x, crossing P transversally at each half-integer. Every lobe
// has opposite sign to its neighbour, so the optimal homotopy must split at
// every crossing, and the total area is exactly N/2 (two 3/4 end-triangles plus
// (N-3) half-triangles in the middle).
std::pair<Curve, Curve> zigzag(int N) {
  Curve p = {P(0, 0), P(static_cast<double>(N), 0)};
  Curve q;
  q.push_back(P(0, 0));
  for (int i = 1; i < N; ++i) q.push_back(P(i, (i % 2 == 1) ? 1.0 : -1.0));
  q.push_back(P(static_cast<double>(N), 0));
  return {p, q};
}

Curve transformed(const Curve& c, double s, double dx, double dy) {
  Curve out;
  for (const auto& v : c)
    out.push_back(P(s * CGAL::to_double(v.x()) + dx,
                    s * CGAL::to_double(v.y()) + dy));
  return out;
}

}  // namespace

TEST_CASE("two anchors: +,-,+ lobes sum to 2.0") {
  auto [p, q] = zigzag(4);
  CHECK(area_of(p, q) == doctest::Approx(2.0));
}

TEST_CASE("three anchors: four alternating lobes sum to 2.5") {
  auto [p, q] = zigzag(5);
  CHECK(area_of(p, q) == doctest::Approx(2.5));
}

TEST_CASE("many anchors follow the N/2 closed form") {
  for (int N : {3, 4, 5, 6, 8, 10}) {
    auto [p, q] = zigzag(N);
    CHECK(area_of(p, q) == doctest::Approx(0.5 * N));
  }
}

TEST_CASE("non-convex consistent region (L-shape) has its polygon area") {
  // P along the x-axis; Q traces the top of an L. P ∘ rev(Q) is a simple
  // notched rectangle of area 4 with no interior crossings -> consistent.
  Curve p = {P(0, 0), P(3, 0)};
  Curve q = {P(0, 0), P(0, 1), P(2, 1), P(2, 2), P(3, 2), P(3, 0)};
  CHECK(area_of(p, q) == doctest::Approx(4.0));
}

TEST_CASE("symmetry holds on the stress cases") {
  std::vector<std::pair<Curve, Curve>> cases = {
      zigzag(4), zigzag(7),
      {{P(0, 0), P(3, 0)},
       {P(0, 0), P(0, 1), P(2, 1), P(2, 2), P(3, 2), P(3, 0)}}};
  for (auto& [p, q] : cases) CHECK(area_of(p, q) == doctest::Approx(area_of(q, p)));
}

TEST_CASE("translation invariance") {
  auto [p, q] = zigzag(6);
  CHECK(area_of(transformed(p, 1, 10, -7), transformed(q, 1, 10, -7)) ==
        doctest::Approx(area_of(p, q)));
}

TEST_CASE("scaling by s multiplies area by s^2") {
  auto [p, q] = zigzag(4);  // area 2.0
  const double s = 3.0;
  CHECK(area_of(transformed(p, s, 0, 0), transformed(q, s, 0, 0)) ==
        doctest::Approx(2.0 * s * s));  // 18.0
}
