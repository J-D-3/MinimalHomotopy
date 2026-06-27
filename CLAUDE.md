# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project computes

Given two simple polygonal curves **P** and **Q** in the plane that share their start and end
points, it computes the **minimal homotopy area**: the smallest total area swept while continuously
deforming P into Q. This is the curve-similarity measure of Erin Wolf Chambers & Yusu Wang,
*"Measuring Similarity Between Curves on 2-Manifolds via Homotopy Area"* (2013) — the paper PDF lives
in [`background/`](background/ChambersWang_article.pdf). We implement the **planar case** only.

The algorithm, in one paragraph: form the closed curve `C = P ∘ rev(Q)`. Its arrangement divides the
plane into faces, each with a *winding number* `wn`. If all winding numbers have the same sign the
loop is *consistent* and its optimal sweep cost is the *total winding number* `Tw = Σ |wn(f)|·area(f)`
(paper Lemma 4.3). If signs are mixed, the curve must be split at an *anchor point* (a shared
intersection every optimal homotopy keeps fixed), and a dynamic program finds the cheapest set of
anchors (paper §4.2).

## Build, test, run

This is a **CMake + vcpkg** project. vcpkg lives at `C:/Repositories/vcpkg` (wired in via
[`CMakePresets.json`](CMakePresets.json)); if it moves, update the preset's `toolchainFile`.

```sh
cmake --preset default          # configure; first run installs CGAL+Boost+GMP+MPFR via vcpkg (slow, one-time)
cmake --build build             # build the library, tests, and the demo CLI
ctest --test-dir build --output-on-failure   # run all tests
```

- **Run a single test case** (doctest filter):
  `./build/mh_tests --test-case="inconsistent winding forces an anchor split"`
  (on Windows the binary is under `build/Release/mh_tests.exe` for the VS generator).
- **List test cases:** `./build/mh_tests --list-test-cases`
- **Demo CLI:** `printf "3\n0 0\n1 0\n1 1\n3\n0 0\n0 1\n1 1\n" | ./build/mh_cli` → prints `1`.

If MSVC fails with "too many sections", the CGAL arrangement headers need `/bigobj` — it is already
set on the library and test targets in [`CMakeLists.txt`](CMakeLists.txt); add it to any new target
that includes the arrangement code.

## Architecture

Everything is under namespace `mh`. Headers in `include/minimal_homotopy/`, sources in `src/`.

| Unit | Role |
|------|------|
| [`kernel.hpp`](include/minimal_homotopy/kernel.hpp) | **The single place the CGAL kernel is named.** Defines `Kernel`, `Point_2`, `Segment_2`, and `Curve = std::vector<Point_2>`. |
| [`result.hpp`](include/minimal_homotopy/result.hpp) | `Result<T>` — tiny value-or-error type (the C++17 stand-in for the old `fplus::result`). |
| [`curve.{hpp,cpp}`](include/minimal_homotopy/curve.hpp) | Input validation: simple-polyline test and the three Chambers/Wang restrictions. |
| [`winding.{hpp,cpp}`](include/minimal_homotopy/winding.hpp) | **The CGAL core.** `analyze_closed_curve` builds the arrangement of one closed loop and returns `Σ|wn|·area` plus a consistency flag. |
| [`homotopy_area.{hpp,cpp}`](include/minimal_homotopy/homotopy_area.hpp) | Public API `calculate_min_homotopy_area(P, Q)`: finds P∩Q intersections and runs the anchor-point DP. |

Data flow: `calculate_min_homotopy_area` → finds & orders the intersections → for each candidate
sub-loop `C[i,j]` calls `analyze_closed_curve` → fills the DP table `T[i]` → returns `T[last]`.

### Key design decisions (read before changing the core)

- **Kernel = EPECK** (`Exact_predicates_exact_constructions_kernel`). Both geometric *decisions*
  (robust arrangement topology) and constructed *values* (intersection points, areas) are exact. The
  whole codebase is kernel-agnostic — for raw speed over exact areas, change the one typedef in
  `kernel.hpp` to EPICK (and restore a tolerance in the `homotopy_area.cpp` dedup; see below).
- **Winding numbers by flood-fill, not ray casting.** `winding.cpp` seeds the unbounded face with
  `wn = 0` and propagates: two faces sharing an edge differ by exactly 1, with the face to the *left*
  of C's direction being higher. Direction is recovered exactly via `Arr_curve_data_traits_2`, which
  tags each curve with its segment index so the data survives edge splitting. This replaces the
  old floating-point angle-summation (`winding_number`) and its epsilon hacks entirely.
- **Exact dedup.** Intersection points constructed from different segment pairs are bit-identical
  under EPECK, so `homotopy_area.cpp` deduplicates and orders them by exact point comparison
  (`PointLess`) and exact along-segment position (squared distance), with no tolerance. Under EPICK
  this must revert to a quantised coordinate key.

### Mapping to the old stump (now removed)

The original `include/MinimalHomtopyArea/` is gone (see git history). The CGAL rewrite replaces it
piece for piece: `graphs.{h,cpp}` Boost.Graph face traversal + `initialize_intersection_lut`/
`build_p2` → CGAL `Arrangement_2` in `winding.cpp`; `winding_number` (angle sum) → flood-fill;
`decompose_at_anchor_points` (recursion) → the forward DP in `homotopy_area.cpp`.

## Input restrictions (enforced by `validate_inputs`)

1. P and Q share start and end points.
2. Each of P, Q is a simple (non-self-intersecting) polyline with no zero-length edges.
3. P and Q do not overlap collinearly — intersections must be transversal.

Caveat: restriction 3 hard-rejects collinear *overlaps* (which would break winding propagation) but
does not fully detect tangential touch-but-not-cross contacts. Inputs are assumed transversal beyond
that. Pass `validate = false` to skip these checks when curves are known-admissible.
