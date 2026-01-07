//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2018-2026 Ken Barker
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
//////////////////////////////////////////////////////////////////////////////
#include "via/sphere.hpp"
#include <boost/test/unit_test.hpp>

using namespace via;
using namespace via::vector;
using namespace via::vector::intersection;

namespace {
constexpr auto EPSILON{std::numeric_limits<double>::epsilon()};
constexpr auto CALCULATION_TOLERANCE{101 * EPSILON};
} // namespace

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(Test_iuntersection_double)

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_calculate_intersection) {
  const LatLong lat_lon_south(Degrees(-90.0), Degrees(0.0));
  const Vector3<double> south_pole{lat_lon_south.to_point()};

  const LatLong lat_lon_north(Degrees(90.0), Degrees(0.0));
  const Vector3<double> north_pole{lat_lon_north.to_point()};

  const LatLong lat_lon_idl(Degrees(0.0), Degrees(180.0));
  const Vector3<double> idl{lat_lon_idl.to_point()};

  const auto equator_intersection{
      calculate_intersection(south_pole, north_pole, MIN_SQ_NORM<double>)};
  BOOST_CHECK(!equator_intersection.has_value());

  const auto gc_intersection1{
      calculate_intersection(idl, north_pole, MIN_SQ_NORM<double>)};
  BOOST_CHECK(gc_intersection1.has_value());
  const auto gc_intersection2{
      calculate_intersection(idl, south_pole, MIN_SQ_NORM<double>)};
  BOOST_CHECK(gc_intersection2.has_value());
  BOOST_CHECK_EQUAL(gc_intersection1.value(), -gc_intersection2.value());
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_use_antipodal_point) {
  const Vector3<double> north_pole(0.0, 0.0, 1.0);
  const Vector3<double> south_pole(0.0, 0.0, -1.0);
  const Vector3<double> g_eq(1.0, 0.0, 0.0);

  BOOST_CHECK(!use_antipodal_point<double>(north_pole, north_pole));
  BOOST_CHECK(!use_antipodal_point<double>(north_pole, g_eq));
  BOOST_CHECK(use_antipodal_point<double>(north_pole, south_pole));
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(
    test_calculate_arc_reference_distances_and_angle_coincident_great_circles) {
  const Vector3<double> point_1(1.0, 0.0, 0.0);
  const Vector3<double> pole_1(0.0, 0.0, 1.0);

  // same mid points and great circles
  const auto result_0 = calculate_arc_reference_distances_and_angle(
      point_1, pole_1, point_1, pole_1);
  BOOST_CHECK_EQUAL(Radians(0.0), get<0>(result_0));
  BOOST_CHECK_EQUAL(Radians(0.0), get<1>(result_0));
  BOOST_CHECK_EQUAL(Degrees(0.0), get<2>(result_0).to_degrees());

  // opposite mid points and same great circles
  const Vector3<double> point_m1{-point_1};
  const auto result_1 = calculate_arc_reference_distances_and_angle(
      point_1, pole_1, point_m1, pole_1);
  BOOST_CHECK_CLOSE(-trig::PI_2<double>, get<0>(result_1).v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(trig::PI_2<double>, get<1>(result_1).v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(Degrees(0.0), get<2>(result_1).to_degrees());

  // opposite mid points and great circles
  const Vector3<double> pole_m1{-pole_1};
  const auto result_2 = calculate_arc_reference_distances_and_angle(
      point_1, pole_1, point_m1, pole_m1);
  BOOST_CHECK_CLOSE(-trig::PI_2<double>, get<0>(result_2).v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(-trig::PI_2<double>, get<1>(result_2).v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(Degrees(180.0), get<2>(result_2).to_degrees());
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(
    test_calculate_arc_reference_distances_and_angle_intersecting_great_circles) {
  const Vector3<double> point_1(1.0, 0.0, 0.0);
  const Vector3<double> pole_1(0.0, 0.0, 1.0);
  const Vector3<double> pole_2(0.0, 1.0, 0.0);

  // intersection, same mid points
  const auto result_0 = calculate_arc_reference_distances_and_angle(
      point_1, pole_1, point_1, pole_2);
  BOOST_CHECK_EQUAL(Radians(0.0), get<0>(result_0));
  BOOST_CHECK_EQUAL(Radians(0.0), get<1>(result_0));
  BOOST_CHECK_EQUAL(Degrees(90.0), get<2>(result_0).to_degrees());

  // intersection, same mid points, acute angle
  const Vector3<double> pole_3{(pole_1 + pole_2).normalized()};
  const auto result_1 = calculate_arc_reference_distances_and_angle(
      point_1, pole_1, point_1, pole_3);
  BOOST_CHECK_EQUAL(Radians(0.0), get<0>(result_1));
  BOOST_CHECK_EQUAL(Radians(0.0), get<1>(result_1));
  BOOST_CHECK_EQUAL(Degrees(45.0), get<2>(result_1).to_degrees());

  // intersection, same mid points, obtuse angle
  const Vector3<double> pole_m3{-pole_3};
  const auto result_2 = calculate_arc_reference_distances_and_angle(
      point_1, pole_1, point_1, pole_m3);
  BOOST_CHECK_EQUAL(Radians(0.0), get<0>(result_2));
  BOOST_CHECK_EQUAL(Radians(0.0), get<1>(result_2));
  BOOST_CHECK_EQUAL(Degrees(135.0), get<2>(result_2).to_degrees());

  // intersection, different mid points, acute angle
  const Vector3<double> point_2{position(point_1, direction(point_1, pole_3),
                                         Angle<double>().quarter_turn_cw())};
  const auto result_3 = calculate_arc_reference_distances_and_angle(
      point_1, pole_1, point_2, pole_3);
  BOOST_CHECK_EQUAL(Radians(0.0), get<0>(result_3));
  BOOST_CHECK_CLOSE(-trig::PI_2<double>, get<1>(result_3).v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(Degrees(45.0), get<2>(result_3).to_degrees());

  // intersection, different mid points, obtuse angle
  const auto result_4 = calculate_arc_reference_distances_and_angle(
      point_1, pole_1, point_2, pole_m3);
  BOOST_CHECK_EQUAL(Radians(0.0), get<0>(result_4));
  BOOST_CHECK_CLOSE(trig::PI_2<double>, get<1>(result_4).v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(Degrees(135.0), get<2>(result_4).to_degrees());
}
//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()
//////////////////////////////////////////////////////////////////////////////
