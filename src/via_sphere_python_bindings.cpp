//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2024 Ken Barker. All Rights Reserved.
// (ken dot barker at via-technology dot co dot uk)
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

  m.def("is_valid_latitude", &via::is_valid_latitude<double>,
        "Test whether a latitude in degrees is a valid latitude.");
  m.def("is_valid_longitude", &via::is_valid_longitude<double>,
        "Test whether a longitude in degrees is a valid longitude.");

  // Python numpy binding for the LatLong class
  PYBIND11_NUMPY_DTYPE(via::LatLong<double>, lat_, lon_);

  // Python bindings for the LatLong class
  py::class_<via::LatLong<double>>(m, "LatLong")
      .def(py::init<via::Degrees<double>, via::Degrees<double>>())

      .def("is_valid", &via::LatLong<double>::is_valid)
      .def("lat", &via::LatLong<double>::lat)
      .def("lon", &via::LatLong<double>::lon)
      .def("__repr__", &via::LatLong<double>::python_repr)

      .def(py::self == py::self);

  m.def("to_point", &via::to_point<double>,
        "Convert a `LatLong` to a point on the unit sphere.");

  m.def("calculate_azimuth_and_distance",
        &via::calculate_azimuth_and_distance<double>,
        "Calculate the azimuth and distance along the great circle of point b from "
        "point a.");

  m.def("haversine_distance",
        &via::haversine_distance<double>,
        "Calculate the distance along the great circle of point b from point a.");

  // Python numpy binding for the Arc class
//   PYBIND11_NUMPY_DTYPE(via::Arc<double>, a_, pole_, length_, half_width_);

  // Python bindings for the Arc class
  py::class_<via::Arc<double>>(m, "Arc")
      .def(py::init<via::vector::Vector3<double>, via::vector::Vector3<double>,
                    via::Radians<double>, via::Radians<double>>())
      .def(py::init<via::LatLong<double>, via::Angle<double>,
                    via::Radians<double>, via::Radians<double>>())
      .def(py::init<via::LatLong<double>, via::LatLong<double>>(),
                    via::Radians<double>>())

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

  m.def("calculate_intersection_distances", &via::calculate_intersection_distances<double>,
        "Calculate the great-circle distances along a pair of `Arc`s to their "
        "closest intersection point or their coincident arc distances if the "
        "`Arc`s are on coincident Great Circles.");
  m.def("calculate_intersection_point", &via::calculate_intersection_point<double>,
        "Calculate whether a pair of `Arc`s intersect and (if so) where.");
}
