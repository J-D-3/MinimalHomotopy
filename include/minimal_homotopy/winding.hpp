// Copyright Ingo Proff 2016-2026.
// https://github.com/CrikeeIP/MinimalHomotopy
// Distributed under the MIT Software License (X11 license).

#pragma once

#include "kernel.hpp"

namespace mh {

// Result of analysing one closed curve C (= a sub-loop P[a,b] ∘ rev(Q[a,b])).
struct WindingAnalysis {
  // Σ |wn(face)| · area(face) over all bounded faces of the arrangement of C.
  // When `consistent` is true this equals |Tw(C)|, the total winding number,
  // which (Chambers/Wang Lemma 4.3) is the optimal homotopy area of that loop.
  double total_area = 0.0;

  // True iff every face has winding number ≥ 0, or every face has wn ≤ 0.
  // Inconsistent loops cannot be swept directly and must be split at an anchor.
  bool consistent = true;
};

// Build the CGAL arrangement of the closed polyline `closed` (which must satisfy
// closed.front() == closed.back()), assign each face its winding number by
// propagating outward from the unbounded face, and return the totals above.
//
// This single function replaces both the hand-rolled intersection/refinement
// machinery and the Boost.Graph planar_face_traversal of the original code.
WindingAnalysis analyze_closed_curve(const Curve& closed);

} // namespace mh
