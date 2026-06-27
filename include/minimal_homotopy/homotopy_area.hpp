// Copyright Ingo Proff 2016-2026.
// https://github.com/CrikeeIP/MinimalHomotopy
// Distributed under the MIT Software License (X11 license).

#pragma once

#include "kernel.hpp"
#include "result.hpp"

namespace mh {

// Minimal homotopy area between two open polylines P and Q that share their
// start and end points: the smallest total area swept while continuously
// deforming P into Q in the plane (Chambers/Wang 2013, planar case).
//
// Set `validate` to false to skip the input-restriction checks (e.g. when the
// caller already knows the curves are admissible).
Result<double> calculate_min_homotopy_area(const Curve& P, const Curve& Q,
                                           bool validate = true);

} // namespace mh
