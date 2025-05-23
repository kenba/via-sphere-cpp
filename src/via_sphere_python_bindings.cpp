//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2018-2025 Ken Barker.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
/// @file via_sphere_python_bindings.cpp
/// @brief Contains the via::sphere python interface
//////////////////////////////////////////////////////////////////////////////
// ensure numpy.h included before sphere.hpp
// clang-format off
#include <pybind11/numpy.h>
#include "via/sphere.hpp"
// clang-format on
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(via_sphere, m) {
  try {
    py::module::import("numpy");
  } catch (...) {
    return;
  }

  // Python bindings for sphere constants
  m.attr("MIN_VALUE") = via::great_circle::MIN_VALUE<double>;
  m.attr("MIN_SIN_ANGLE") = via::vector::MIN_SIN_ANGLE<double>;
  m.attr("MIN_SQ_NORM") = via::vector::MIN_SQ_NORM<double>;
  m.attr("MIN_SQ_DISTANCE") = via::vector::MIN_SQ_DISTANCE<double>;

  // Python bindings for great_circle functions
  m.def("e2gc_distance", &via::great_circle::e2gc_distance<double>,
        "Convert a Euclidean distance to a Great Circle distance (in "
        "`Radians`).");
  m.def(
      "gc2e_distance", &via::great_circle::gc2e_distance<double>,
      "Convert a Great Circle distance (in `Radians`) to a Euclidean distance");
  m.def("sq_euclidean_distance",
        &via::great_circle::sq_euclidean_distance<double>,
        "Calculate the square of the Euclidean distance (i.e. using "
        "Pythagoras) between two points from their Latitudes and their "
        "Longitude difference.");
  m.def("calculate_gc_distance",
        &via::great_circle::calculate_gc_distance<double>,
        "Calculate the Great Circle distance (angle from centre) between two "
        "points from their Latitudes and their Longitude difference.");
  m.def(
      "calculate_gc_azimuth", &via::great_circle::calculate_gc_azimuth<double>,
      "Calculate the azimuth (bearing) along the great circle of point b from  "
      "point a from their Latitudes and their Longitude difference.");
  m.def("calculate_sigma", &via::great_circle::calculate_sigma<double>,
        "Calculate the Great Circle distance (as an angle, sigma) between two "
        "points from their Latitudes and their Longitude difference.");
  m.def("calculate_latitude", &via::great_circle::calculate_latitude<double>,
        "Calculate the latitude at great circle distance, sigma.");
  m.def("calculate_delta_longitude",
        &via::great_circle::calculate_delta_longitude<double>,
        "Calculate the longitude difference at great circle distance, sigma.");

  // Python bindings for Vector3 functions
  m.def("perp_product", &via::vector::perp_product<double>,
        "2d vector perp product function: a x b.");
  m.def("dot2d", &via::vector::dot2d<double>,
        "2d vector dot product function: a . b.");
  m.def("to_point", &via::vector::to_point<double>,
        "Convert a latitude and longitude to a point on the unit sphere.");
  m.def("latitude", &via::vector::latitude<double>,
        "Calculate the latitude of a point on the unit sphere.");
  m.def("longitude", &via::vector::longitude<double>,
        "Calculate the longitude of a point on the unit sphere.");
  m.def("sq_distance", &via::vector::sq_distance<double>,
        "Calculate the square of the Euclidean distance between two points.");
  m.def("distance", &via::vector::distance<double>,
        "Calculate the shortest (Euclidean) distance between two points.");
  m.def("are_orthogonal", &via::vector::are_orthogonal<double>,
        "Determine whether two `Vector3`s are orthogonal (perpendicular).");
  m.def("delta_longitude", &via::vector::delta_longitude<double>,
        "Calculate the relative longitude of point a from point b.");
  m.def("is_west_of", &via::vector::is_west_of<double>,
        "Determine whether point a is West of point b.");
  m.def("calculate_pole", &via::vector::calculate_pole<double>,
        "Calculate the right hand pole vector of a Great Circle from an "
        "initial position and an azimuth.");
  m.def(
      "calculate_azimuth", &via::vector::calculate_azimuth<double>,
      "Calculate the azimuth at a point on the Great Circle defined by pole.");
  m.def("calculate_direction", &via::vector::calculate_direction<double>,
        "Calculate the direction vector along a Great Circle from an initial "
        "position and an azimuth.");
  m.def("direction", &via::vector::direction<double>,
        "Calculate the direction vector of a Great Circle arc.");
  m.def("position", &via::vector::position<double>,
        "Calculate the position of a point along a Great Circle arc.");
  m.def("rotate", &via::vector::rotate<double>,
        "Calculate the direction vector of a Great Circle rotated by angle.");
  m.def("rotate_position", &via::vector::rotate_position<double>,
        "Calculate the position of a point rotated by angle at radius.");
  m.def("sin_xtd", &via::vector::sin_xtd<double>,
        "The sine of the across track distance of a point relative to a Great "
        "Circle pole.");
  m.def(
      "cross_track_distance", &via::vector::cross_track_distance<double>,
      "The across track distance of a point relative to a Great Circle pole.");
  m.def("sq_cross_track_distance",
        &via::vector::sq_cross_track_distance<double>,
        "The square of the Euclidean cross track distance of a point relative "
        "to a Great Circle pole.");
  m.def("calculate_point_on_plane",
        &via::vector::calculate_point_on_plane<double>,
        "Calculate the closest point on a plane to the given point.");
  m.def("sin_atd", &via::vector::sin_atd<double>,
        "The sine of the along track distance of a point along a Great Circle "
        "arc.");
  m.def("calculate_great_circle_atd",
        &via::vector::calculate_great_circle_atd<double>,
        "Calculate the relative distance of two points on a Great Circle arc.");
  m.def("along_track_distance", &via::vector::along_track_distance<double>,
        "The Great Circle distance of a point along the arc relative to a, "
        "(+ve) ahead of a, (-ve) behind a.");
  m.def("sq_along_track_distance",
        &via::vector::sq_along_track_distance<double>,
        "Calculate the square of the Euclidean along track distance of a point "
        "from the start of an Arc.");
  m.def("calculate_atd_and_xtd", &via::vector::calculate_atd_and_xtd<double>,
        "Calculate the along track and across track distance of a point from "
        "the start of an Arc.");

  // Python bindings for vector intersection functions
  m.def("calculate_intersection",
        &via::vector::intersection::calculate_intersection<double>,
        "Calculate an intersection point between the poles of two Great "
        "Circles.");
  m.def("calculate_intersection_distances",
        &via::vector::intersection::calculate_intersection_distances<double>,
        "Calculate the great circle distances to an intersection point from "
        "the start points of a pair of great circle arcs, on different great "
        "circles.");

  // Python bindings for sphere functions
  m.def("is_valid_latitude", &via::is_valid_latitude<double>,
        "Test whether a latitude in degrees is a valid latitude.");
  m.def("is_valid_longitude", &via::is_valid_longitude<double>,
        "Test whether a longitude in degrees is a valid longitude.");

  // Python numpy binding for the LatLong class
  PYBIND11_NUMPY_DTYPE(via::LatLong<double>, lat_, lon_);

  // Python bindings for the LatLong class
  py::class_<via::LatLong<double>>(m, "LatLong")
      .def(py::init<via::Degrees<double>, via::Degrees<double>>())
      .def(py::init<via::vector::Vector3<double>>())

      .def("is_valid", &via::LatLong<double>::is_valid)
      .def("lat", &via::LatLong<double>::lat)
      .def("lon", &via::LatLong<double>::lon)
      .def("to_point", &via::LatLong<double>::to_point)
      .def("__repr__", &via::LatLong<double>::python_repr)

      .def(py::self == py::self);

  m.def("calculate_azimuth_and_distance",
        &via::calculate_azimuth_and_distance<double>,
        "Calculate the azimuth and distance along the great circle of point b "
        "from point a.");

  m.def(
      "haversine_distance", &via::haversine_distance<double>,
      "Calculate the distance along the great circle of point b from point a.");

  // Python bindings for the Arc class
  py::class_<via::Arc<double>>(m, "Arc")
      .def(py::init<via::vector::Vector3<double>, via::vector::Vector3<double>,
                    via::Radians<double>, via::Radians<double>>())
      .def(py::init<via::LatLong<double>, via::Angle<double>,
                    via::Radians<double>, via::Radians<double>>())
      .def(py::init<via::LatLong<double>, via::LatLong<double>,
                    via::Radians<double>>())
      .def(py::init<via::LatLong<double>, via::LatLong<double>>())

      .def("set_half_width", &via::Arc<double>::set_half_width)
      .def("is_valid", &via::Arc<double>::is_valid)

      .def("a", &via::Arc<double>::a)
      .def("pole", &via::Arc<double>::pole)
      .def("length", &via::Arc<double>::length)
      .def("half_width", &via::Arc<double>::half_width)
      .def("azimuth", &via::Arc<double>::azimuth)
      .def("direction", &via::Arc<double>::direction)
      .def("position", &via::Arc<double>::position)
      .def("b", &via::Arc<double>::b)
      .def("mid_point", &via::Arc<double>::mid_point)
      .def("perp_position", &via::Arc<double>::perp_position)
      .def("angle_position", &via::Arc<double>::angle_position)
      .def("end_arc", &via::Arc<double>::end_arc)
      .def("calculate_atd_and_xtd", &via::Arc<double>::calculate_atd_and_xtd)

      .def("__repr__", &via::Arc<double>::python_repr);

  m.def("calculate_intersection_distances",
        &via::calculate_intersection_distances<double>,
        "Calculate the great-circle distances along a pair of `Arc`s to their "
        "closest intersection point or their coincident arc distances if the "
        "`Arc`s are on coincident Great Circles.");
  m.def("calculate_intersection_point",
        &via::calculate_intersection_point<double>,
        "Calculate whether a pair of `Arc`s intersect and (if so) where.");
}
