cmake_minimum_required(VERSION 3.20)
project(MinimalHomotopy LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# CGAL is header-only since 5.0, but its config target (CGAL::CGAL) brings in the
# non-header-only dependencies (Boost, GMP, MPFR). vcpkg provides all of them.
find_package(CGAL CONFIG REQUIRED)

# ---- Core library -----------------------------------------------------------
add_library(minimal_homotopy
  src/curve.cpp
  src/winding.cpp
  src/homotopy_area.cpp
)
target_include_directories(minimal_homotopy
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(minimal_homotopy PUBLIC CGAL::CGAL)

# CGAL's arrangement templates produce very large object files on MSVC.
if(MSVC)
  target_compile_options(minimal_homotopy PRIVATE /bigobj)
endif()

# ---- Demo CLI ---------------------------------------------------------------
add_executable(mh_cli apps/main.cpp)
target_link_libraries(mh_cli PRIVATE minimal_homotopy)

# ---- Tests ------------------------------------------------------------------
enable_testing()
find_package(doctest CONFIG REQUIRED)
add_executable(mh_tests
  test/test_main.cpp
  test/test_winding.cpp
  test/test_homotopy_area.cpp
)
target_link_libraries(mh_tests PRIVATE minimal_homotopy doctest::doctest)
if(MSVC)
  target_compile_options(mh_tests PRIVATE /bigobj)
endif()
add_test(NAME mh_tests COMMAND mh_tests)
