// Copyright Ingo Proff 2016-2026.
// https://github.com/CrikeeIP/MinimalHomotopy
// Distributed under the MIT Software License (X11 license).
//
// The dynamic program of Chambers/Wang (section 4.2). Replaces the recursive
// decompose_at_anchor_points of the original code with the paper's forward DP.

#include "minimal_homotopy/homotopy_area.hpp"

#include "minimal_homotopy/curve.hpp"
#include "minimal_homotopy/winding.hpp"

#include <CGAL/intersections.h>
#include <CGAL/version.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <tuple>
#include <variant>
#include <vector>

namespace mh {
namespace {

// One intersection point of P and Q, with where it sits along each curve.
struct Inter {
  Point_2     pt;
  std::size_t p_seg = 0;  // index of the P-segment it lies on, + param along it
  double      p_t   = 0.0;
  std::size_t q_seg = 0;  // index of the Q-segment it lies on, + param along it
  double      q_t   = 0.0;
  std::size_t q_rank = 0; // order of this point along Q (0 = Q's start)
};

template <class T, class Variant>
const T* held(const Variant* v) {
#if CGAL_VERSION_MAJOR >= 6
  return std::get_if<T>(v);
#else
  return boost::get<T>(v);
#endif
}

// Parameter of x along segment a->b, in [0,1] for points on the segment. Only
// used to order coincident-segment intersections, so inexact arithmetic is fine.
double seg_param(const Point_2& a, const Point_2& b, const Point_2& x) {
  const double dx = CGAL::to_double(b.x() - a.x());
  const double dy = CGAL::to_double(b.y() - a.y());
  const double len2 = dx * dx + dy * dy;
  if (len2 <= 0.0) return 0.0;
  const double px = CGAL::to_double(x.x() - a.x());
  const double py = CGAL::to_double(x.y() - a.y());
  return (px * dx + py * dy) / len2;
}

// Quantised coordinate key, used to deduplicate intersection points that were
// constructed (in floating point) from different segment pairs. This tolerance
// is the price of the EPICK kernel; switching to EPECK would make it exact.
std::pair<std::int64_t, std::int64_t> coord_key(const Point_2& p) {
  auto r = [](double v) { return static_cast<std::int64_t>(std::llround(v * 1e7)); };
  return {r(CGAL::to_double(p.x())), r(CGAL::to_double(p.y()))};
}

// Sub-polyline of `poly` from intersection A to intersection B (A before B).
Curve extract(const Curve& poly, std::size_t from_seg, const Point_2& from_pt,
              std::size_t to_seg, const Point_2& to_pt) {
  Curve out;
  out.push_back(from_pt);
  for (std::size_t s = from_seg; s < to_seg; ++s) out.push_back(poly[s + 1]);
  out.push_back(to_pt);
  return out;
}

// The closed loop C[A,B] = P[A,B] ∘ rev(Q[A,B]).
Curve build_closed(const Curve& P, const Curve& Q, const Inter& A, const Inter& B) {
  Curve sp = extract(P, A.p_seg, A.pt, B.p_seg, B.pt);  // A..B along P
  Curve sq = extract(Q, A.q_seg, A.pt, B.q_seg, B.pt);  // A..B along Q

  Curve closed = sp;
  for (auto it = sq.rbegin() + 1; it != sq.rend(); ++it)  // rev(Q), drop dup B
    closed.push_back(*it);

  // Collapse any consecutive duplicate points so we never emit zero-length edges.
  Curve dedup;
  for (const auto& p : closed)
    if (dedup.empty() || !(dedup.back() == p)) dedup.push_back(p);
  return dedup;
}

} // namespace

Result<double> calculate_min_homotopy_area(const Curve& P, const Curve& Q,
                                           bool validate) {
  if (P.size() < 2 || Q.size() < 2)
    return Result<double>::failure("Polylines must have at least 2 points.");
  // A curve is maximally similar to itself: no deformation, zero swept area.
  // (Handled up front because identical curves overlap everywhere, which the
  // transversal-intersection machinery below is not meant to process.)
  if (P == Q) return Result<double>::success(0.0);
  if (!(P.front() == Q.front()) || !(P.back() == Q.back()))
    return Result<double>::failure("P and Q must share start and end points.");
  if (validate) {
    auto v = validate_inputs(P, Q);
    if (!v.ok()) return Result<double>::failure(v.error);
  }

  // --- 1. Collect P∩Q intersections (brute force; predicates are exact). ----
  std::vector<Inter> raw;
  for (std::size_t a = 0; a + 1 < P.size(); ++a) {
    Segment_2 sp(P[a], P[a + 1]);
    for (std::size_t b = 0; b + 1 < Q.size(); ++b) {
      Segment_2 sq(Q[b], Q[b + 1]);
      auto res = CGAL::intersection(sp, sq);
      if (!res) continue;
      if (held<Segment_2>(&*res))
        return Result<double>::failure(
            "P and Q overlap collinearly (intersections must be transversal).");
      const Point_2* pt = held<Point_2>(&*res);
      if (!pt) continue;
      Inter it;
      it.pt = *pt;
      it.p_seg = a;
      it.q_seg = b;
      it.p_t = seg_param(P[a], P[a + 1], *pt);
      it.q_t = seg_param(Q[b], Q[b + 1], *pt);
      raw.push_back(it);
    }
  }

  // --- 2. Deduplicate, then order along P and along Q. ----------------------
  std::map<std::pair<std::int64_t, std::int64_t>, Inter> uniq;
  for (const auto& it : raw) {
    auto key = coord_key(it.pt);
    auto f = uniq.find(key);
    if (f == uniq.end())
      uniq.emplace(key, it);
    else if (std::tie(it.p_seg, it.p_t) <
             std::tie(f->second.p_seg, f->second.p_t))
      f->second = it;  // keep the lowest (p_seg, p_t) representative
  }

  std::vector<Inter> X;
  X.reserve(uniq.size());
  for (auto& kv : uniq) X.push_back(kv.second);

  // Order along P (this becomes the index of each point in the DP).
  std::sort(X.begin(), X.end(), [](const Inter& a, const Inter& b) {
    return std::tie(a.p_seg, a.p_t) < std::tie(b.p_seg, b.p_t);
  });

  // Order along Q -> q_rank, looked up by point key.
  std::vector<Inter> byq = X;
  std::sort(byq.begin(), byq.end(), [](const Inter& a, const Inter& b) {
    return std::tie(a.q_seg, a.q_t) < std::tie(b.q_seg, b.q_t);
  });
  std::map<std::pair<std::int64_t, std::int64_t>, std::size_t> qrank;
  for (std::size_t i = 0; i < byq.size(); ++i) qrank[coord_key(byq[i].pt)] = i;
  for (auto& x : X) x.q_rank = qrank[coord_key(x.pt)];

  if (X.size() < 2)
    return Result<double>::failure(
        "Expected the shared endpoints as intersections; found fewer.");

  // --- 3. Dynamic program (Chambers/Wang 4.2). ------------------------------
  // T[i] = optimal homotopy area between P[x0,xi] and Q[x0,xi].
  const double INF = std::numeric_limits<double>::infinity();
  const std::size_t last = X.size() - 1;
  std::vector<double> T(X.size(), INF);
  T[0] = 0.0;

  for (std::size_t i = 1; i < X.size(); ++i) {
    // If the whole prefix loop is consistent, its total winding is optimal.
    auto whole = analyze_closed_curve(build_closed(P, Q, X[0], X[i]));
    if (whole.consistent) {
      T[i] = whole.total_area;
      continue;
    }
    // Otherwise split at the best valid anchor x_j (same order along P and Q,
    // with a consistent loop C[j,i]).
    double best = INF;
    for (std::size_t j = 1; j < i; ++j) {
      if (X[j].q_rank >= X[i].q_rank) continue;  // must keep order along Q
      if (!std::isfinite(T[j])) continue;
      auto part = analyze_closed_curve(build_closed(P, Q, X[j], X[i]));
      if (!part.consistent) continue;
      best = std::min(best, part.total_area + T[j]);
    }
    T[i] = best;
  }

  if (!std::isfinite(T[last]))
    return Result<double>::failure("No valid anchor decomposition found.");
  return Result<double>::success(T[last]);
}

} // namespace mh
