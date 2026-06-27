// Copyright Ingo Proff 2016-2026.
// https://github.com/CrikeeIP/MinimalHomotopy
// Distributed under the MIT Software License (X11 license).
//
// Arrangement construction + winding numbers. This file is the CGAL-heavy core
// that replaces the original graphs.{h,cpp} (Boost.Graph planar_face_traversal)
// and the hand-rolled intersection/refinement code.
//
// Idea: build the arrangement of the closed curve C. Each arrangement edge is
// traced by exactly one directed arc of C (true because intersections are
// transversal — no overlaps). The winding number of two faces sharing an edge
// differs by exactly 1: the face to the LEFT of C's direction is the higher
// one. So we seed the unbounded face with wn = 0 and flood-fill outward.

#include "minimal_homotopy/winding.hpp"

#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

#include <cmath>
#include <cstddef>
#include <limits>
#include <queue>
#include <vector>

namespace mh {
namespace {

// Segment traits carrying, per curve, the index of its directed C-segment.
// The data rides along when CGAL splits edges at intersection points, so every
// resulting half-edge can be traced back to the arc of C that induced it.
using SegTraits = CGAL::Arr_segment_traits_2<Kernel>;
using Traits    = CGAL::Arr_curve_data_traits_2<SegTraits, std::size_t>;
using Dcel      = CGAL::Arr_face_extended_dcel<Traits, int>;  // int = winding nr
using Arr       = CGAL::Arrangement_2<Traits, Dcel>;

constexpr int UNVISITED = std::numeric_limits<int>::min();

// Shoelace signed area of one connected boundary component (CCB), walking the
// CCB via next() from `start`. A face's outer CCB is counter-clockwise
// (positive); its holes are clockwise (negative), so summing all CCBs of a face
// yields its true positive area.
double signed_area_ccb(Arr::Halfedge_handle start) {
  double a = 0.0;
  Arr::Halfedge_handle h = start;
  do {
    const auto& p = h->source()->point();
    const auto& q = h->target()->point();
    a += CGAL::to_double(p.x()) * CGAL::to_double(q.y()) -
         CGAL::to_double(q.x()) * CGAL::to_double(p.y());
    h = h->next();
  } while (h != start);
  return 0.5 * a;
}

// Collect every half-edge of the CCB starting at `start`, walking next()
// pointers and comparing handle identity (robust, unlike circulator !=).
void collect_ccb(Arr::Halfedge_handle start,
                 std::vector<Arr::Halfedge_handle>& out) {
  Arr::Halfedge_handle h = start;
  do {
    out.push_back(h);
    h = h->next();
  } while (h != start);
}

double face_area(Arr::Face_handle f) {
  if (f->is_unbounded()) return 0.0;
  double a = signed_area_ccb(f->outer_ccb());
  for (auto h = f->holes_begin(); h != f->holes_end(); ++h)
    a += signed_area_ccb(*h);  // *h converts to a half-edge handle on the hole
  return a;
}

// The outer ring of a bounded face as a list of its vertices (no closing dup).
Curve outer_ring(Arr::Halfedge_handle start) {
  Curve ring;
  Arr::Halfedge_handle h = start;
  do {
    ring.push_back(h->target()->point());
    h = h->next();
  } while (h != start);
  return ring;
}

// Does half-edge `h` run in the same direction as the C-arc that induced it?
// The two vectors are parallel by construction, so the sign of their dot
// product is unambiguous (no floating-point knife's edge).
bool aligned_with_curve(Arr::Halfedge_const_handle h,
                        const std::vector<Segment_2>& segs) {
  const std::size_t k = h->curve().data();
  const auto        d = segs[k].to_vector();
  const auto&       s = h->source()->point();
  const auto&       t = h->target()->point();
  const double dot = CGAL::to_double((t.x() - s.x()) * d.x() +
                                     (t.y() - s.y()) * d.y());
  return dot > 0.0;
}

} // namespace

WindingAnalysis analyze_closed_curve(const Curve& closed) {
  WindingAnalysis out;

  // Directed segments of C (skipping any zero-length steps), each tagged with
  // its own index so we can recover its direction after edge splitting.
  std::vector<Segment_2>       segs;
  std::vector<Traits::Curve_2> curves;
  for (std::size_t i = 0; i + 1 < closed.size(); ++i) {
    if (closed[i] == closed[i + 1]) continue;
    const std::size_t k = segs.size();
    segs.emplace_back(closed[i], closed[i + 1]);
    curves.emplace_back(SegTraits::Curve_2(closed[i], closed[i + 1]), k);
  }
  if (segs.empty()) return out;  // degenerate: zero area, trivially consistent

  Arr arr;
  CGAL::insert(arr, curves.begin(), curves.end());

  // Flood-fill winding numbers outward from the unbounded face.
  for (auto f = arr.faces_begin(); f != arr.faces_end(); ++f)
    f->set_data(UNVISITED);

  Arr::Face_handle uf = arr.unbounded_face();
  uf->set_data(0);
  std::queue<Arr::Face_handle> bfs;
  bfs.push(uf);

  while (!bfs.empty()) {
    Arr::Face_handle f = bfs.front();
    bfs.pop();

    // Collect every boundary half-edge of f: its outer CCB (if bounded) plus
    // the CCB of each hole. Crossing any of them reaches a neighbouring face.
    std::vector<Arr::Halfedge_handle> boundary;
    if (!f->is_unbounded()) collect_ccb(f->outer_ccb(), boundary);
    for (auto hit = f->holes_begin(); hit != f->holes_end(); ++hit)
      collect_ccb(*hit, boundary);

    for (Arr::Halfedge_handle h : boundary) {
      Arr::Face_handle g = h->twin()->face();
      if (g->data() != UNVISITED) continue;
      // h's incident face f lies to its LEFT. If h runs along C, the left face
      // has the higher winding number, so the neighbour g is one lower.
      const int delta = aligned_with_curve(h, segs) ? -1 : +1;
      g->set_data(f->data() + delta);
      bfs.push(g);
    }
  }

  // Total winding area Σ|wn|·area, and the consistency test.
  int sign = 0;
  for (auto f = arr.faces_begin(); f != arr.faces_end(); ++f) {
    if (f->is_unbounded()) continue;
    int wn = f->data();
    if (wn == UNVISITED) wn = 0;  // safeguard; shouldn't happen for a connected C
    if (wn != 0) {
      const int s = wn > 0 ? 1 : -1;
      if (sign == 0)
        sign = s;
      else if (sign != s)
        out.consistent = false;
      out.faces.push_back({outer_ring(f->outer_ccb()), wn});
    }
    out.total_area += std::abs(wn) * face_area(f);
  }

  return out;
}

} // namespace mh
