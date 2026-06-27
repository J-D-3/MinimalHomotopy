// Copyright Ingo Proff 2016-2026.
// https://github.com/CrikeeIP/MinimalHomotopy
// Distributed under the MIT Software License (X11 license).

#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <vector>

namespace mh {

// The one and only place the CGAL kernel is named.
//
// EPECK = Exact_predicates_exact_constructions_kernel: both geometric
// *decisions* (do two segments cross? which side of a line is a point on?) and
// *constructed* values (intersection points, areas) are exact. Intersection
// points built from different segment pairs are therefore bit-identical, so we
// can compare and deduplicate them exactly — no tolerance fudge.
//
// EPECK trades speed for that exactness (its FT is a lazy exact number). If a
// future use needs raw throughput over exact areas, switch this single typedef
// to CGAL::Exact_predicates_inexact_constructions_kernel (EPICK) — the rest of
// the code is written kernel-agnostic. (Under EPICK, restore a tolerance in the
// intersection dedup in homotopy_area.cpp, which currently assumes exactness.)
using Kernel    = CGAL::Exact_predicates_exact_constructions_kernel;
using Point_2   = Kernel::Point_2;
using Segment_2 = Kernel::Segment_2;

// An open polygonal curve (polyline): an ordered list of points.
using Curve = std::vector<Point_2>;

} // namespace mh
