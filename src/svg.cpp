// Copyright Ingo Proff 2016-2026.
// https://github.com/CrikeeIP/MinimalHomotopy
// Distributed under the MIT Software License (X11 license).

#include "minimal_homotopy/svg.hpp"

#include <CGAL/number_utils.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>

namespace mh {
namespace {

// World <-> SVG mapping. SVG y points down, so we flip; a fixed margin keeps
// strokes off the edge.
struct View {
  double minx, miny, maxx, maxy, scale, margin;
  double sx(double x) const { return margin + (x - minx) * scale; }
  double sy(double y) const { return margin + (maxy - y) * scale; }  // y-flip
  double width()  const { return (maxx - minx) * scale + 2 * margin; }
  double height() const { return (maxy - miny) * scale + 2 * margin; }
};

void grow(const Curve& c, double& minx, double& miny, double& maxx, double& maxy) {
  for (const auto& p : c) {
    const double x = CGAL::to_double(p.x()), y = CGAL::to_double(p.y());
    minx = std::min(minx, x); maxx = std::max(maxx, x);
    miny = std::min(miny, y); maxy = std::max(maxy, y);
  }
}

std::string points_attr(const Curve& c, const View& v) {
  std::ostringstream o;
  for (std::size_t i = 0; i < c.size(); ++i) {
    if (i) o << ' ';
    o << v.sx(CGAL::to_double(c[i].x())) << ',' << v.sy(CGAL::to_double(c[i].y()));
  }
  return o.str();
}

// Centroid of a ring (vertex average) in SVG coordinates — good enough to place
// a small multiplicity label.
std::pair<double, double> ring_centroid(const Curve& ring, const View& v) {
  double sx = 0, sy = 0;
  for (const auto& p : ring) {
    sx += v.sx(CGAL::to_double(p.x()));
    sy += v.sy(CGAL::to_double(p.y()));
  }
  const double n = static_cast<double>(ring.size());
  return {sx / n, sy / n};
}

} // namespace

std::string render_svg(const Curve& P, const Curve& Q,
                       const HomotopyDecomposition& decomp) {
  double minx = std::numeric_limits<double>::max(), miny = minx;
  double maxx = std::numeric_limits<double>::lowest(), maxy = maxx;
  grow(P, minx, miny, maxx, maxy);
  grow(Q, minx, miny, maxx, maxy);
  for (const auto& f : decomp.swept) grow(f.boundary, minx, miny, maxx, maxy);
  if (minx > maxx) { minx = miny = 0; maxx = maxy = 1; }  // empty guard

  const double span = std::max({maxx - minx, maxy - miny, 1e-9});
  View view{minx, miny, maxx, maxy, 600.0 / span, 40.0};

  std::ostringstream s;
  s << std::fixed;
  s.precision(2);
  s << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << view.width()
    << "\" height=\"" << view.height() << "\" viewBox=\"0 0 " << view.width()
    << ' ' << view.height() << "\" font-family=\"sans-serif\">\n";

  // Diagonal hatch used to fill the swept regions.
  s << "  <defs>\n"
       "    <pattern id=\"hatch\" patternUnits=\"userSpaceOnUse\" width=\"7\" "
       "height=\"7\" patternTransform=\"rotate(45)\">\n"
       "      <line x1=\"0\" y1=\"0\" x2=\"0\" y2=\"7\" stroke=\"#2a9d8f\" "
       "stroke-width=\"1.4\"/>\n"
       "    </pattern>\n"
       "  </defs>\n";

  s << "  <rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

  // Swept faces: hatched. Drawing a face |wn| times deepens the hatch where the
  // homotopy sweeps with higher multiplicity.
  for (const auto& f : decomp.swept) {
    const std::string pts = points_attr(f.boundary, view);
    const int reps = std::max(1, std::abs(f.winding));
    for (int r = 0; r < reps; ++r)
      s << "  <polygon points=\"" << pts
        << "\" fill=\"url(#hatch)\" fill-opacity=\"0.55\" stroke=\"#2a9d8f\" "
           "stroke-width=\"1\" stroke-opacity=\"0.6\"/>\n";
    if (std::abs(f.winding) > 1) {
      auto [cx, cy] = ring_centroid(f.boundary, view);
      s << "  <text x=\"" << cx << "\" y=\"" << cy
        << "\" font-size=\"15\" fill=\"#1d6f64\" text-anchor=\"middle\">×"
        << std::abs(f.winding) << "</text>\n";
    }
  }

  // The two curves, in distinct colours.
  s << "  <polyline points=\"" << points_attr(P, view)
    << "\" fill=\"none\" stroke=\"#1f77b4\" stroke-width=\"3\" "
       "stroke-linejoin=\"round\" stroke-linecap=\"round\"/>\n";
  s << "  <polyline points=\"" << points_attr(Q, view)
    << "\" fill=\"none\" stroke=\"#e4572e\" stroke-width=\"3\" "
       "stroke-linejoin=\"round\" stroke-linecap=\"round\"/>\n";

  // Shared endpoints (anchors of every homotopy).
  for (const Point_2* p : {&P.front(), &P.back()})
    s << "  <circle cx=\"" << view.sx(CGAL::to_double(p->x())) << "\" cy=\""
      << view.sy(CGAL::to_double(p->y())) << "\" r=\"4\" fill=\"#333\"/>\n";

  // Legend.
  const double lx = view.margin, ly = view.height() - 14;
  s << "  <text x=\"" << lx << "\" y=\"" << ly
    << "\" font-size=\"14\" fill=\"#1f77b4\">P</text>\n";
  s << "  <text x=\"" << lx + 18 << "\" y=\"" << ly
    << "\" font-size=\"14\" fill=\"#e4572e\">Q</text>\n";
  s << "  <text x=\"" << lx + 40 << "\" y=\"" << ly
    << "\" font-size=\"14\" fill=\"#1d6f64\">minimal homotopy area = " << decomp.area
    << "</text>\n";

  s << "</svg>\n";
  return s.str();
}

} // namespace mh
