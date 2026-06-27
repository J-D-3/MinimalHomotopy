// Copyright Ingo Proff 2016-2026.
// https://github.com/CrikeeIP/MinimalHomotopy
// Distributed under the MIT Software License (X11 license).
//
// Reads two polylines from stdin (same format as mh_cli) and writes an SVG of
// the two curves with the swept minimal-homotopy area hatched to stdout.
//
//   printf "3\n0 0\n1 0\n1 1\n3\n0 0\n0 1\n1 1\n" | mh_svg > out.svg

#include "minimal_homotopy/homotopy_area.hpp"
#include "minimal_homotopy/svg.hpp"

#include <iostream>

namespace {

mh::Curve read_curve(std::istream& in) {
  std::size_t n = 0;
  in >> n;
  mh::Curve c;
  c.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    double x = 0.0, y = 0.0;
    in >> x >> y;
    c.emplace_back(x, y);
  }
  return c;
}

}  // namespace

int main() {
  mh::Curve p = read_curve(std::cin);
  mh::Curve q = read_curve(std::cin);

  auto d = mh::decompose_min_homotopy_area(p, q);
  if (!d.ok()) {
    std::cerr << "error: " << d.error << '\n';
    return 1;
  }
  std::cout << mh::render_svg(p, q, *d.value);
  return 0;
}
