// Copyright Ingo Proff 2016-2026.
// https://github.com/CrikeeIP/MinimalHomotopy
// Distributed under the MIT Software License (X11 license).

#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <vector>

namespace mh {

// The one and only place the CGAL kernel is named.
//
// EPICK = Exact_predicates_inexact_constructions_kernel: geometric *decisions*
// (do two segments cross? which side of a line is a point on?) are exact, while
// *constructed* values (intersection points, areas) are ordinary doubles. That
// makes the arrangement topology robust while keeping arithmetic fast.
//
// If area sums ever look unstable on near-degenerate inputs, switch this single
// typedef to CGAL::Exact_predicates_exact_constructions_kernel (EPECK) — the
// rest of the code is written kernel-agnostic and needs no other change.
using Kernel    = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2   = Kernel::Point_2;
using Segment_2 = Kernel::Segment_2;

// An open polygonal curve (polyline): an ordered list of points.
using Curve = std::vector<Point_2>;

} // namespace mh
