// Copyright Ingo Proff 2016-2026.
// https://github.com/CrikeeIP/MinimalHomotopy
// Distributed under the MIT Software License (X11 license).

#pragma once

#include "kernel.hpp"
#include "result.hpp"

namespace mh {

// True if `poly` is a simple (non-self-intersecting) open polyline with no
// degenerate (zero-length) edges. Reproduces the old `cw_admissible_polygon`.
bool is_simple_polyline(const Curve& poly);

// Checks the Chambers/Wang input restrictions on the pair (P, Q):
//   1. P and Q share their start and end points,
//   2. each of P, Q is a simple polyline,
//   3. P and Q do not overlap collinearly (intersections must be transversal).
// Returns ok() on success, otherwise an explanatory error.
//
// Note: full transversality (rejecting tangential touch-but-not-cross contacts)
// is only partially enforced — collinear overlaps, which would break the
// winding computation, are the case we reject hard. See README/CLAUDE.md.
Result<bool> validate_inputs(const Curve& P, const Curve& Q);

} // namespace mh
