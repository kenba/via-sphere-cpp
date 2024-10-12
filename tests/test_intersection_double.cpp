//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2024 Ken Barker
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
#include <via/trig.hpp>

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
BOOST_AUTO_TEST_CASE(test_calculate_intersection_point) {
  const LatLong lat_lon_south(Degrees(-90.0), Degrees(0.0));
  const Vector3<double> south_pole{to_point(lat_lon_south)};

  const LatLong lat_lon_north(Degrees(90.0), Degrees(0.0));
  const Vector3<double> north_pole{to_point(lat_lon_north)};

  const LatLong lat_lon_idl(Degrees(0.0), Degrees(180.0));
  const Vector3<double> idl{to_point(lat_lon_idl)};

  const auto equator_intersection{
      calculate_intersection_point(south_pole, north_pole)};
  BOOST_CHECK(!equator_intersection.has_value());

  const auto gc_intersection1{calculate_intersection_point(idl, north_pole)};
  BOOST_CHECK(gc_intersection1.has_value());
  const auto gc_intersection2{calculate_intersection_point(idl, south_pole)};
  BOOST_CHECK(gc_intersection2.has_value());
  BOOST_CHECK_EQUAL(gc_intersection1.value(), -gc_intersection2.value());
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_calculate_intersection_distances) {
  const LatLong start1(Degrees(-1.0), Degrees(-1.0));
  const Vector3<double> a1{to_point(start1)};
  const Angle<double> azimuth1(Degrees(45.0));
  const Vector3<double> pole1{
      calculate_pole(Angle(start1.lat()), Angle(start1.lon()), azimuth1)};

  const LatLong start2(Degrees(1.0), Degrees(-1.0));
  const Vector3<double> a2{to_point(start2)};
  const Angle<double> azimuth2(Degrees(135.0));
  const Vector3<double> pole2{
      calculate_pole(Angle(start2.lat()), Angle(start2.lon()), azimuth2)};

  const auto c{calculate_intersection_point(pole1, pole2)};
  BOOST_CHECK(c.has_value());
  const auto [c1, c2]{
      calculate_intersection_distances(a1, pole1, a2, pole2, c.value())};
  BOOST_CHECK_CLOSE(-3.1169124762478333, c1.v(), CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(-3.1169124762478333, c2.v(), CALCULATION_TOLERANCE);

  // Calculate the centre of the arc start points
  const Vector3<double> sum{a1 + a2};
  const auto centre_point{normalise(sum)};
  BOOST_CHECK(centre_point.has_value());
  BOOST_CHECK(sq_distance(c.value(), centre_point.value()) > 2.0);

  // opposite intersection point
  const Vector3<double> d{-c.value()};
  BOOST_CHECK(sq_distance(d, centre_point.value()) <= 2.0);
  const auto [d1,
              d2]{calculate_intersection_distances(a1, pole1, a2, pole2, d)};
  BOOST_CHECK_CLOSE(0.024680177341956263, d1.v(), CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(0.024680177341956263, d2.v(), CALCULATION_TOLERANCE);

  // Same start points and intersection point
  const auto [e1,
              e2]{calculate_intersection_distances(a1, pole1, a1, pole1, a1)};
  BOOST_CHECK_EQUAL(Radians(0.0), e1);
  BOOST_CHECK_EQUAL(Radians(0.0), e2);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_is_within) {
  BOOST_CHECK(!is_within(-2.0 * EPSILON, 2.0));
  BOOST_CHECK(is_within(-EPSILON, 2.0));
  BOOST_CHECK(is_within(2.0 * (1.0 + EPSILON), 2.0));
  BOOST_CHECK(!is_within(2.0 * (1.0 + 3.0 * EPSILON), 2.0));
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_calculate_coincident_arc_distances) {
  const Radians zero{0.0};
  const Radians length1{0.25};
  const Radians length2{0.75};

  const auto result0{
      calculate_coincident_arc_distances(length2, true, length2, length1)};
  BOOST_CHECK_EQUAL(length2, get<0>(result0));
  BOOST_CHECK_EQUAL(zero, get<1>(result0));

  const auto result1{
      calculate_coincident_arc_distances(length2, true, length1, length2)};
  BOOST_CHECK_EQUAL(zero, get<0>(result1));
  BOOST_CHECK_EQUAL(length2, get<1>(result1));

  const auto result2{
      calculate_coincident_arc_distances(Radians(1.0), true, length1, length2)};
  BOOST_CHECK_EQUAL(length1, get<0>(result2));
  BOOST_CHECK_EQUAL(length2, get<1>(result2));

  const auto result3{
      calculate_coincident_arc_distances(Radians(1.5), true, length1, length2)};
  BOOST_CHECK_EQUAL(length2, get<0>(result3));
  BOOST_CHECK_EQUAL(length2, get<1>(result3));

  const auto result4{calculate_coincident_arc_distances(Radians(-1.5), true,
                                                        length1, length2)};
  BOOST_CHECK_EQUAL(Radians(-1.5), get<0>(result4));
  BOOST_CHECK_EQUAL(zero, get<1>(result4));

  const auto result5{calculate_coincident_arc_distances(Radians(-1.0), false,
                                                        length1, length2)};
  BOOST_CHECK_EQUAL(zero, get<0>(result5));
  BOOST_CHECK_EQUAL(Radians(1.0), get<1>(result5));

  const auto result6{calculate_coincident_arc_distances(Radians(1.0), false,
                                                        length1, length2)};
  BOOST_CHECK_EQUAL(Radians(1.0), get<0>(result6));
  BOOST_CHECK_EQUAL(zero, get<1>(result6));

  const auto result7{
      calculate_coincident_arc_distances(-length2, false, length1, length2)};
  BOOST_CHECK_EQUAL(zero, get<0>(result7));
  BOOST_CHECK_EQUAL(length2, get<1>(result7));

  const auto result8{
      calculate_coincident_arc_distances(length1, false, length1, length2)};
  BOOST_CHECK_EQUAL(length1, get<0>(result8));
  BOOST_CHECK_EQUAL(zero, get<1>(result8));
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

BOOST_AUTO_TEST_SUITE_END()
//////////////////////////////////////////////////////////////////////////////
