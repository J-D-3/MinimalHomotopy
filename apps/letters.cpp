// Copyright Ingo Proff 2016-2026.
// https://github.com/CrikeeIP/MinimalHomotopy
// Distributed under the MIT Software License (X11 license).
//
// Experiment: can the minimal homotopy area recognise a letter across fonts?
//
// Each glyph is a single pen stroke drawn between the same two anchor points
// (top-left -> bottom-right of a shared box), so any two glyphs share their
// endpoints and the homotopy area between them is defined. A "font" bends the
// interior of the stroke by an endpoint-preserving sine offset. We then ask
// whether same-letter-across-fonts is closer (smaller area) than different
// letters. Writes comparison SVGs and prints the area table.

#include "minimal_homotopy/homotopy_area.hpp"
#include "minimal_homotopy/svg.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace mh;

namespace {

const double PI = 3.14159265358979323846;
Point_2 Pt(double x, double y) { return Point_2(x, y); }

// A "font": offset every interior vertex by amp*sin(pi*t) along a fixed generic
// direction (20 deg), leaving the shared endpoints fixed. Different amp ->
// different font (a different glyph if amp differs a lot). The direction is
// chosen off the glyph segment angles (0/45/90 deg) so no segment stays
// collinear across fonts (which would be an inadmissible overlap).
Curve font(const Curve& base, double amp) {
  Curve out = base;
  const std::size_t n = base.size() - 1;
  const double dx = 0.93969, dy = 0.34202;  // (cos 20, sin 20)
  for (std::size_t i = 1; i < n; ++i) {
    const double t = static_cast<double>(i) / static_cast<double>(n);
    const double s = std::sin(PI * t);
    out[i] = Pt(CGAL::to_double(base[i].x()) + amp * s * dx,
                CGAL::to_double(base[i].y()) + amp * s * dy);
  }
  return out;
}

double area(const Curve& p, const Curve& q) {
  auto r = calculate_min_homotopy_area(p, q);
  return r.ok() ? *r : -1.0;  // -1 marks an inadmissible pair
}

void write_svg(const std::string& path, const Curve& p, const Curve& q) {
  auto d = decompose_min_homotopy_area(p, q);
  std::ofstream o(path);
  if (d.ok()) o << render_svg(p, q, *d.value);
  else        std::cerr << "  (skip " << path << ": " << d.error << ")\n";
}

}  // namespace

int main(int argc, char** argv) {
  const std::string dir = (argc > 1) ? argv[1] : "examples/letters";

  // Base glyphs: single strokes from TL (0,100) to BR (100,0).
  const Curve L = {Pt(0, 100), Pt(0, 0), Pt(100, 0)};            // down then right
  const Curve seven = {Pt(0, 100), Pt(100, 100), Pt(100, 0)};   // right then down
  const Curve Z = {Pt(0, 100), Pt(100, 100), Pt(0, 0), Pt(100, 0)};
  const Curve S = {Pt(0, 100), Pt(75, 80), Pt(25, 55), Pt(75, 30), Pt(100, 0)};

  const double A = 8, B = -8, BIG = 78;  // two normal fonts + one extreme font

  std::cout << "same-letter (font A vs font B):\n";
  std::cout << "  L: " << area(font(L, A), font(L, B)) << "\n";
  std::cout << "  7: " << area(font(seven, A), font(seven, B)) << "\n";
  std::cout << "  Z: " << area(font(Z, A), font(Z, B)) << "\n";
  std::cout << "  S: " << area(font(S, A), font(S, B)) << "\n";

  std::cout << "different letters (both font A):\n";
  std::cout << "  L vs 7: " << area(font(L, A), font(seven, A)) << "\n";
  std::cout << "  Z vs S: " << area(font(Z, A), font(S, A)) << "\n";

  std::cout << "extreme font:\n";
  std::cout << "  Z(A) vs Z(BIG): " << area(font(Z, A), font(Z, BIG)) << "\n";

  // Narrative images.
  write_svg(dir + "/success_same_L.svg", font(L, A), font(L, B));
  write_svg(dir + "/success_diff_L_7.svg", font(L, A), font(seven, A));
  write_svg(dir + "/failure_diff_Z_S.svg", font(Z, A), font(S, A));
  write_svg(dir + "/failure_same_Z_distorted.svg", font(Z, A), font(Z, BIG));
  return 0;
}
