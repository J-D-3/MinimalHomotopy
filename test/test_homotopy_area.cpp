// Copyright Ingo Proff 2016-2026.
// https://github.com/CrikeeIP/MinimalHomotopy
// Distributed under the MIT Software License (X11 license).

#include <doctest/doctest.h>

#include "minimal_homotopy/homotopy_area.hpp"

using namespace mh;

namespace {
Point_2 P(double x, double y) { return Point_2(x, y); }
}  // namespace

TEST_CASE("disjoint curves bounding a unit square -> area 1") {
  // P and Q meet only at the shared endpoints; P ∘ rev(Q) is a unit square.
  Curve p = {P(0, 0), P(1, 0), P(1, 1)};
  Curve q = {P(0, 0), P(0, 1), P(1, 1)};
  auto r = calculate_min_homotopy_area(p, q);
  REQUIRE(r.ok());
  CHECK(*r == doctest::Approx(1.0));
}

TEST_CASE("measure is symmetric in its arguments") {
  Curve p = {P(0, 0), P(1, 0), P(1, 1)};
  Curve q = {P(0, 0), P(0, 1), P(1, 1)};
  auto a = calculate_min_homotopy_area(p, q);
  auto b = calculate_min_homotopy_area(q, p);
  REQUIRE(a.ok());
  REQUIRE(b.ok());
  CHECK(*a == doctest::Approx(*b));
}

TEST_CASE("a straight P vs a tent Q -> triangle area") {
  Curve p = {P(0, 0), P(2, 0)};
  Curve q = {P(0, 0), P(1, 1), P(2, 0)};  // tent of height 1 over base 2
  auto r = calculate_min_homotopy_area(p, q);
  REQUIRE(r.ok());
  CHECK(*r == doctest::Approx(1.0));
}

TEST_CASE("inconsistent winding forces an anchor split") {
  // Q zig-zags across P once at (1.5, 0). The whole loop has a +1 lobe and a -1
  // lobe (cancelling to Tw = 0), so it is inconsistent and must be split at the
  // crossing into two triangles of area 0.75 each -> total 1.5.
  Curve p = {P(0, 0), P(3, 0)};
  Curve q = {P(0, 0), P(1, 1), P(2, -1), P(3, 0)};
  auto r = calculate_min_homotopy_area(p, q);
  REQUIRE(r.ok());
  CHECK(*r == doctest::Approx(1.5));
}

TEST_CASE("identical curves -> zero area") {
  Curve p = {P(0, 0), P(1, 1), P(2, 0)};
  Curve q = p;
  auto r = calculate_min_homotopy_area(p, q, /*validate=*/false);
  REQUIRE(r.ok());
  CHECK(*r == doctest::Approx(0.0));
}

TEST_CASE("validation rejects mismatched endpoints") {
  Curve p = {P(0, 0), P(1, 0), P(1, 1)};
  Curve q = {P(0, 0), P(0, 1), P(2, 2)};  // different end point
  auto r = calculate_min_homotopy_area(p, q);
  CHECK_FALSE(r.ok());
}

TEST_CASE("validation rejects a self-intersecting input") {
  Curve p = {P(0, 0), P(2, 2), P(2, 0), P(0, 2), P(3, 3)};  // self-crossing P
  Curve q = {P(0, 0), P(1, 0), P(3, 3)};
  auto r = calculate_min_homotopy_area(p, q);
  CHECK_FALSE(r.ok());
}
