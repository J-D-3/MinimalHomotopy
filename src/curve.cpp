// Copyright Ingo Proff 2016-2026.
// https://github.com/CrikeeIP/MinimalHomotopy
// Distributed under the MIT Software License (X11 license).

#include "minimal_homotopy/curve.hpp"

#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/intersections.h>
#include <CGAL/version.h>

#include <variant>
#include <vector>

namespace mh {
namespace {

using SegTraits = CGAL::Arr_segment_traits_2<Kernel>;
using Arr       = CGAL::Arrangement_2<SegTraits>;

// Pull the held alternative out of the result of CGAL::intersection, regardless
// of whether this CGAL version returns a std::variant (>= 6.0) or a
// boost::variant (5.x).
template <class T, class Variant>
const T* held(const Variant* v) {
#if CGAL_VERSION_MAJOR >= 6
  return std::get_if<T>(v);
#else
  return boost::get<T>(v);
#endif
}

} // namespace

bool is_simple_polyline(const Curve& poly) {
  if (poly.size() < 2) return false;

  std::vector<Segment_2> segs;
  segs.reserve(poly.size() - 1);
  for (std::size_t i = 0; i + 1 < poly.size(); ++i) {
    if (poly[i] == poly[i + 1]) return false;  // degenerate edge
    segs.emplace_back(poly[i], poly[i + 1]);
  }

  // For a simple path the arrangement is exactly that path: one vertex per
  // point and one edge per segment. Any self-crossing or self-touch introduces
  // extra vertices and splits edges, so the counts no longer match.
  Arr arr;
  CGAL::insert(arr, segs.begin(), segs.end());
  return arr.number_of_vertices() == poly.size() &&
         arr.number_of_edges() == poly.size() - 1;
}

Result<bool> validate_inputs(const Curve& P, const Curve& Q) {
  if (P.size() < 2 || Q.size() < 2)
    return Result<bool>::failure("Polylines must have at least 2 points.");

  if (!(P.front() == Q.front()) || !(P.back() == Q.back()))
    return Result<bool>::failure("P and Q must share start and end points.");

  if (!is_simple_polyline(P))
    return Result<bool>::failure("P is not a simple polyline.");
  if (!is_simple_polyline(Q))
    return Result<bool>::failure("Q is not a simple polyline.");

  // Reject collinear overlaps between P and Q: an overlap means a piece of the
  // closed curve is traced twice, which breaks the winding-number propagation.
  for (std::size_t a = 0; a + 1 < P.size(); ++a) {
    Segment_2 sp(P[a], P[a + 1]);
    for (std::size_t b = 0; b + 1 < Q.size(); ++b) {
      Segment_2 sq(Q[b], Q[b + 1]);
      auto res = CGAL::intersection(sp, sq);
      if (res && held<Segment_2>(&*res))
        return Result<bool>::failure(
            "P and Q overlap collinearly (intersections must be transversal).");
    }
  }

  return Result<bool>::success(true);
}

} // namespace mh
