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

namespace {
constexpr auto EPSILON{std::numeric_limits<double>::epsilon()};
constexpr auto CALCULATION_TOLERANCE{101 * EPSILON};
} // namespace

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(Test_vector_double)

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_normalise) {
  const Vector3<double> zero(0.0, 0.0, 0.0);
  BOOST_CHECK(!normalise(zero));

  // Greenwich equator
  const Vector3<double> g_eq(1.0, 0.0, 0.0);
  BOOST_CHECK(normalise(g_eq));

  // A vector just too small to normalize
  const Vector3<double> too_small(
      16383 * std::numeric_limits<double>::epsilon(), 0.0, 0.0);
  BOOST_CHECK(!normalise(too_small));

  // A vector just large enough to normalize
  const Vector3<double> small(16384 * std::numeric_limits<double>::epsilon(),
                              0.0, 0.0);
  BOOST_CHECK(!is_unit(small));
  const auto result{normalise(small)};
  BOOST_CHECK(result.has_value());
  BOOST_CHECK(is_unit(result.value()));
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_point_lat_longs) {
  // Test South pole
  const LatLong lat_lon_south(Degrees(-90.0), Degrees(180.0));
  const Vector3<double> point_south{lat_lon_south.to_point()};

  BOOST_CHECK_EQUAL(Degrees(-90.0), latitude(point_south).to_degrees());
  BOOST_CHECK_EQUAL(Degrees(0.0), longitude(point_south).to_degrees());

  LatLong result{point_south};
  BOOST_CHECK_EQUAL(Degrees(-90.0), result.lat());
  // Note: longitude is now zero, since the poles do not have a Longitude
  BOOST_CHECK_EQUAL(Degrees(0.0), result.lon());

  // Test Greenwich equator
  const LatLong lat_lon_0_0(Degrees(0.0), Degrees(0.0));
  const Vector3<double> point_0{lat_lon_0_0.to_point()};
  BOOST_CHECK(is_unit(point_0));
  BOOST_CHECK_EQUAL(lat_lon_0_0, LatLong(point_0));

  // Test IDL equator
  const LatLong lat_lon_0_180(Degrees(0.0), Degrees(180.0));
  const Vector3<double> point_1{lat_lon_0_180.to_point()};
  BOOST_CHECK(is_unit(point_1));
  BOOST_CHECK_EQUAL(lat_lon_0_180, LatLong(point_1));
  BOOST_CHECK(!is_west_of(point_0, point_1));
  BOOST_CHECK_EQUAL(Radians(-trig::PI<double>),
                    delta_longitude(point_0, point_1).to_radians());

  const LatLong lat_lon_0_m180(Degrees(0.0), Degrees(-180.0));
  const Vector3<double> point_2{lat_lon_0_m180.to_point()};
  BOOST_CHECK(is_unit(point_2));
  // Converts back to +ve longitude
  BOOST_CHECK_EQUAL(lat_lon_0_180, LatLong(point_2));
  BOOST_CHECK(!is_west_of(point_0, point_2));
  BOOST_CHECK_EQUAL(Radians(-trig::PI<double>),
                    delta_longitude(point_0, point_2).to_radians());

  const LatLong lat_lon_0_r3(Degrees(0.0), Degrees(trig::rad2deg(3.0)));
  const Vector3<double> point_3{lat_lon_0_r3.to_point()};
  BOOST_CHECK(is_unit(point_3));

  result = LatLong(point_3);
  BOOST_CHECK_EQUAL(Degrees(0.0), result.lat());
  BOOST_CHECK_EQUAL(Degrees(trig::rad2deg(3.0)), result.lon());
  BOOST_CHECK(is_west_of(point_0, point_3));
  BOOST_CHECK_EQUAL(Radians(-3.0),
                    delta_longitude(point_0, point_3).to_radians());

  BOOST_CHECK(!is_west_of(point_1, point_3));
  BOOST_CHECK_CLOSE(trig::PI<double> - 3.0,
                    delta_longitude(point_1, point_3).to_radians().v(),
                    CALCULATION_TOLERANCE);

  const LatLong lat_lon_0_mr3(Degrees(0.0), Degrees(trig::rad2deg(-3.0)));
  const Vector3<double> point_4{lat_lon_0_mr3.to_point()};
  BOOST_CHECK(is_unit(point_4));
  BOOST_CHECK_EQUAL(Radians(3.0),
                    delta_longitude(point_0, point_4).to_radians());

  result = LatLong(point_4);
  BOOST_CHECK_EQUAL(Degrees(0.0), result.lat());
  BOOST_CHECK_EQUAL(Degrees(trig::rad2deg(-3.0)), result.lon());
  BOOST_CHECK(is_west_of(point_1, point_4));
  BOOST_CHECK_CLOSE(3.0 - trig::PI<double>,
                    delta_longitude(point_1, point_4).to_radians().v(),
                    CALCULATION_TOLERANCE);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_point_distance) {
  const LatLong lat_lon_south(Degrees(-90.0), Degrees(180.0));
  const Vector3<double> south_pole{lat_lon_south.to_point()};

  const LatLong lat_lon_north(Degrees(90.0), Degrees(180.0));
  const Vector3<double> north_pole{lat_lon_north.to_point()};

  BOOST_CHECK_EQUAL(0.0, sq_distance(south_pole, south_pole));
  BOOST_CHECK_EQUAL(0.0, sq_distance(north_pole, north_pole));
  BOOST_CHECK_EQUAL(4.0, sq_distance(south_pole, north_pole));

  BOOST_CHECK_EQUAL(0.0, distance(south_pole, south_pole));
  BOOST_CHECK_EQUAL(0.0, distance(north_pole, north_pole));
  BOOST_CHECK_EQUAL(2.0, distance(south_pole, north_pole));

  // Greenwich equator
  const Vector3<double> g_eq(1.0, 0.0, 0.0);

  // Test IDL equator
  const Vector3<double> idl_eq(-1.0, 0.0, 0.0);

  BOOST_CHECK_EQUAL(0.0, sq_distance(g_eq, g_eq));
  BOOST_CHECK_EQUAL(0.0, sq_distance(idl_eq, idl_eq));
  BOOST_CHECK_EQUAL(4.0, sq_distance(g_eq, idl_eq));

  BOOST_CHECK_EQUAL(0.0, distance(g_eq, g_eq));
  BOOST_CHECK_EQUAL(0.0, distance(idl_eq, idl_eq));
  BOOST_CHECK_EQUAL(2.0, distance(g_eq, idl_eq));
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_calculate_azimuth_at_poles) {
  // Greenwich equator
  const Vector3<double> g_eq(1.0, 0.0, 0.0);
  const Vector3<double> south_pole(0.0, 0.0, -1.0);
  auto result = calculate_azimuth(south_pole, g_eq);
  BOOST_CHECK_EQUAL(Angle(Degrees(0.0)), result);

  const Vector3<double> north_pole(0.0, 0.0, 1.0);
  result = calculate_azimuth(north_pole, g_eq);
  BOOST_CHECK_EQUAL(Angle(Degrees(180.0)), result);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_calculate_pole_azimuth_and_direction) {
  // Greenwich equator
  const Vector3<double> g_eq(1.0, 0.0, 0.0);

  // 90 degrees East on the equator
  const Vector3<double> e_eq(0.0, 1.0, 0.0);

  // 90 degrees West on the equator
  const Vector3<double> w_eq(0.0, -1.0, 0.0);

  const Angle<double> angle_90(Degrees(90.0));
  const auto pole_a{calculate_pole(Angle<double>(), Angle<double>(), angle_90)};
  BOOST_CHECK(are_orthogonal(g_eq, pole_a));

  const auto dir_a{
      calculate_direction(Angle<double>(), Angle<double>(), angle_90)};
  BOOST_CHECK(are_orthogonal(g_eq, dir_a));
  BOOST_CHECK(are_orthogonal(pole_a, dir_a));
  BOOST_CHECK_EQUAL(dir_a, direction(g_eq, pole_a));

  const Vector3<double> north_pole(0.0, 0.0, 1.0);
  BOOST_CHECK_EQUAL(north_pole, pole_a);

  auto result{calculate_azimuth(g_eq, pole_a)};
  BOOST_CHECK_EQUAL(angle_90, result);

  const auto pole_b{
      calculate_pole(Angle<double>(), Angle<double>(), -angle_90)};
  BOOST_CHECK(are_orthogonal(g_eq, pole_b));

  const auto dir_b{
      calculate_direction(Angle<double>(), Angle<double>(), -angle_90)};
  BOOST_CHECK(are_orthogonal(g_eq, dir_b));
  BOOST_CHECK(are_orthogonal(pole_a, dir_b));
  BOOST_CHECK_EQUAL(dir_b, direction(g_eq, pole_b));

  const Vector3<double> south_pole(0.0, 0.0, -1.0);
  BOOST_CHECK_EQUAL(south_pole, pole_b);

  BOOST_CHECK_EQUAL(south_pole, g_eq.cross(w_eq));

  result = calculate_azimuth(g_eq, pole_b);
  BOOST_CHECK_EQUAL(-angle_90, result);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_calculate_position) {
  // Greenwich equator
  const Vector3<double> g_eq(1.0, 0.0, 0.0);

  // 90 degrees East on the equator
  const Vector3<double> e_eq(0.0, 1.0, 0.0);

  const Vector3<double> pole_0{g_eq.cross(e_eq)};

  const Angle<double> angle_90(Degrees(90.0));

  const Vector3<double> pos_1{
      position(g_eq, direction(g_eq, pole_0), angle_90)};
  BOOST_CHECK_EQUAL(e_eq, pos_1);

  const Vector3<double> pos_2{
      rotate_position(g_eq, pole_0, Angle<double>(), angle_90)};
  BOOST_CHECK_EQUAL(e_eq, pos_2);

  const Vector3<double> pos_3{
      rotate_position(g_eq, pole_0, angle_90, angle_90)};
  BOOST_CHECK_EQUAL(pole_0, pos_3);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_calculate_cross_track_distance_and_square) {
  // Greenwich equator
  const Vector3<double> g_eq(1.0, 0.0, 0.0);

  // 90 degrees East on the equator
  const Vector3<double> e_eq(0.0, 1.0, 0.0);

  const Vector3<double> pole_0{g_eq.cross(e_eq)};

  const Degrees<double> longitude{1.0};

  for (auto latitude = -89; latitude < 90; ++latitude) {
    const LatLong<double> latlong(Degrees<double>(latitude), longitude);
    const Vector3<double> point{latlong.to_point()};

    const double expected{trig::deg2rad(static_cast<double>(latitude))};
    const auto xtd{cross_track_distance(pole_0, point)};

    // Accuracy reduces outside of this range
    const auto tolerance{(std::abs(latitude) <= 83
                              ? 2 * CALCULATION_TOLERANCE
                              : 32 * CALCULATION_TOLERANCE)};

    BOOST_CHECK_CLOSE(expected, xtd.v(), tolerance);

    const double expected_e{
        great_circle::gc2e_distance(Radians<double>(expected))};
    const double sq_expected{expected_e * expected_e};
    const auto xtd2{sq_cross_track_distance(pole_0, point)};
    BOOST_CHECK_SMALL(sq_expected - xtd2, CALCULATION_TOLERANCE);
  }
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_calculate_along_track_distance_and_square) {
  // Greenwich equator
  const Vector3<double> g_eq(1.0, 0.0, 0.0);

  // 90 degrees East on the equator
  const Vector3<double> e_eq(0.0, 1.0, 0.0);

  const Vector3<double> pole_0{g_eq.cross(e_eq)};

  // North of Equator
  const Degrees<double> latitude{1.0};

  for (auto longitude = -179; longitude < 180; ++longitude) {
    const LatLong<double> latlong(latitude, Degrees<double>(longitude));
    const Vector3<double> point{latlong.to_point()};

    const double expected{trig::deg2rad(static_cast<double>(longitude))};
    const auto atd{along_track_distance(g_eq, pole_0, point)};

    // Accuracy reduces outside of this range
    const auto tolerance{(std::abs(longitude) <= 153
                              ? 2 * CALCULATION_TOLERANCE
                              : 32 * CALCULATION_TOLERANCE)};

    BOOST_CHECK_CLOSE(expected, atd.v(), tolerance);

    const auto [atd1, xtd1] = calculate_atd_and_xtd(g_eq, pole_0, point);
    BOOST_CHECK_CLOSE(expected, atd1.v(), tolerance);
    BOOST_CHECK_CLOSE(trig::deg2rad(1.0), xtd1.v(), tolerance);

    const double expected_e{
        great_circle::gc2e_distance(Radians<double>(expected))};
    const double sq_expected{expected_e * expected_e};
    const double atd2 = sq_along_track_distance(g_eq, pole_0, point);
    BOOST_CHECK_SMALL(sq_expected - atd2, tolerance);
  }
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_special_cases) {
  // Greenwich equator
  const Vector3<double> g_eq(1.0, 0.0, 0.0);

  // 90 degrees East on the equator
  const Vector3<double> e_eq(0.0, 1.0, 0.0);

  const Vector3<double> pole_0{g_eq.cross(e_eq)};

  // points are at the poles, so atc and sq_atd are zero
  BOOST_CHECK_EQUAL(0.0, along_track_distance(g_eq, pole_0, pole_0).v());
  BOOST_CHECK_EQUAL(0.0, sq_along_track_distance(g_eq, pole_0, pole_0));

  const auto [atd, xtd] = calculate_atd_and_xtd(g_eq, pole_0, g_eq);
  BOOST_CHECK_EQUAL(0.0, atd.v());
  BOOST_CHECK_EQUAL(0.0, xtd.v());

  const auto [atd1, xtd1] = calculate_atd_and_xtd(g_eq, pole_0, pole_0);
  BOOST_CHECK_EQUAL(0.0, atd1.v());
  BOOST_CHECK_EQUAL(trig::PI_2<double>, xtd1.v());

  // Test for 100% code coverage
  const double offset{0.000001};
  const LatLong<double> near_north_pole(Degrees(90.0 - offset), Degrees(0.0));
  const Vector3<double> p{near_north_pole.to_point()};
  const auto [atd2, xtd2] = calculate_atd_and_xtd(g_eq, pole_0, p);
  BOOST_CHECK_EQUAL(0.0, atd2.v());
  BOOST_CHECK_CLOSE(trig::PI_2<double>, xtd2.v(), 100 * offset);
}
//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()
//////////////////////////////////////////////////////////////////////////////
