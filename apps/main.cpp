// Copyright Ingo Proff 2016-2026.
// https://github.com/CrikeeIP/MinimalHomotopy
// Distributed under the MIT Software License (X11 license).
//
// Tiny demo CLI. Reads two polylines from stdin and prints their minimal
// homotopy area. Input format:
//
//   n
//   x0 y0
//   ... (n lines)
//   m
//   x0 y0
//   ... (m lines)
//
// Example (unit square): printf "3\n0 0\n1 0\n1 1\n3\n0 0\n0 1\n1 1\n" | mh_cli

#include "minimal_homotopy/homotopy_area.hpp"

#include <iostream>
#include <vector>

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

  auto result = mh::calculate_min_homotopy_area(p, q);
  if (!result.ok()) {
    std::cerr << "error: " << result.error << '\n';
    return 1;
  }
  std::cout << *result << '\n';
  return 0;
}
