// Copyright Ingo Proff 2016-2026.
// https://github.com/CrikeeIP/MinimalHomotopy
// Distributed under the MIT Software License (X11 license).

#pragma once

#include "homotopy_area.hpp"
#include "kernel.hpp"

#include <string>

namespace mh {

// Render a standalone SVG document showing P and Q in distinct colours with the
// swept minimal-homotopy area hatched. `decomp` is the result of
// decompose_min_homotopy_area(P, Q); its faces are drawn as the hatched regions.
std::string render_svg(const Curve& P, const Curve& Q,
                       const HomotopyDecomposition& decomp);

} // namespace mh
