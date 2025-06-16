//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2018-2025 Ken Barker
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

namespace {
constexpr auto EPSILON{std::numeric_limits<double>::epsilon()};
constexpr auto CALCULATION_TOLERANCE{101 * EPSILON};
} // namespace

using Vector3d = vector::Vector3<double>;

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(Test_sphere_double)

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_is_valid_latitude) {
  // value < -90
  BOOST_CHECK(!is_valid_latitude(-90.0001));
  // value = -90
  BOOST_CHECK(is_valid_latitude(-90.0));
  // value = 90
  BOOST_CHECK(is_valid_latitude(90.0));
  // value > 90
  BOOST_CHECK(!is_valid_latitude(90.0001));
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_is_valid_longitude) {
  // value < -180
  BOOST_CHECK(!is_valid_longitude(-180.0001));
  // value = -180
  BOOST_CHECK(is_valid_longitude(-180.0));
  // value = 180
  BOOST_CHECK(is_valid_longitude(180.0));
  // value > 180
  BOOST_CHECK(!is_valid_longitude(180.0001));
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_latlong_class) {
  const LatLong<double> a(Degrees(0.0), Degrees(90.0));
  BOOST_CHECK(a.is_valid());
  BOOST_CHECK_EQUAL(Degrees(0.0), a.lat());
  BOOST_CHECK_EQUAL(Degrees(90.0), a.lon());

  //
  // const LatLong<double> invalid_lat(Degrees(91.0), Degrees(0.0));
  // const LatLong<double> invalid_lon(Degrees(0.0), Degrees(181.0));
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_vector3_double) {
  const LatLong<double> a(Degrees(0.0), Degrees(90.0));
  const auto point{a.to_point()};

  BOOST_CHECK_EQUAL(0.0, point(0));
  BOOST_CHECK_EQUAL(1.0, point(1));
  BOOST_CHECK_EQUAL(0.0, point(2));

  BOOST_CHECK_EQUAL(Degrees(0.0), vector::latitude(point).to_degrees());
  BOOST_CHECK_EQUAL(Degrees(90.0), vector::longitude(point).to_degrees());

  BOOST_CHECK_EQUAL(a, LatLong<double>(point));
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_great_circle_90n_0n_0e) {
  const LatLong<double> a(Degrees(90.0), Degrees(0.0));
  const LatLong<double> b(Degrees(0.0), Degrees(0.0));
  const auto [azimuth, distance]{calculate_azimuth_and_distance(a, b)};

  BOOST_CHECK_CLOSE(trig::PI_2<double>, distance.v(), CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(Degrees(180.0), azimuth.to_degrees());

  const auto dist{haversine_distance(a, b)};
  BOOST_CHECK_CLOSE(trig::PI_2<double>, dist.v(), CALCULATION_TOLERANCE);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_great_circle_90s_0n_50e) {
  const LatLong<double> a(Degrees(-90.0), Degrees(0.0));
  const LatLong<double> b(Degrees(0.0), Degrees(50.0));
  const auto [azimuth, distance]{calculate_azimuth_and_distance(a, b)};

  BOOST_CHECK_CLOSE(trig::PI_2<double>, distance.v(), CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(Degrees(0.0), azimuth.to_degrees());

  const auto dist{haversine_distance(a, b)};
  BOOST_CHECK_CLOSE(trig::PI_2<double>, dist.v(), CALCULATION_TOLERANCE);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_great_circle_0n_60e_0n_60w) {
  const LatLong<double> a(Degrees(0.0), Degrees(60.0));
  const LatLong<double> b(Degrees(0.0), Degrees(-60.0));
  const auto [azimuth, distance]{calculate_azimuth_and_distance(a, b)};

  BOOST_CHECK_CLOSE(2 * trig::PI_3<double>, distance.v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(Degrees(-90.0), azimuth.to_degrees());

  const auto dist{haversine_distance(a, b)};
  BOOST_CHECK_CLOSE(2 * trig::PI_3<double>, dist.v(), CALCULATION_TOLERANCE);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_arc) {
  // Greenwich equator
  const LatLong<double> g_eq(Degrees(0.0), Degrees(0.0));

  // 90 degrees East on the equator
  const LatLong<double> e_eq(Degrees(0.0), Degrees(90.0));

  const Arc<double> arc(g_eq, e_eq, Radians(0.01));
  BOOST_CHECK(arc.is_valid());
  BOOST_CHECK_EQUAL(Radians(0.01), arc.half_width());

  BOOST_CHECK_EQUAL(g_eq.to_point(), arc.a());
  BOOST_CHECK_EQUAL(Vector3d(0.0, 0.0, 1.0), arc.pole());
  BOOST_CHECK_CLOSE(trig::PI_2<double>, arc.length().v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(Degrees(90.0), arc.azimuth().to_degrees());
  const auto b{e_eq.to_point()};
  BOOST_CHECK_SMALL(vector::distance(b, arc.b()), CALCULATION_TOLERANCE);

  const auto mid_point{arc.mid_point()};
  BOOST_CHECK_EQUAL(0.0, mid_point(2));
  BOOST_CHECK_CLOSE(45.0, vector::longitude(mid_point).to_degrees().v(),
                    CALCULATION_TOLERANCE);

  const auto start_arc{arc.end_arc(false)};
  BOOST_CHECK_EQUAL(Radians(0.02), start_arc.length());

  const auto start_arc_a{start_arc.a()};
  BOOST_CHECK_EQUAL(start_arc_a, arc.perp_position(arc.a(), Radians(0.01)));

  const Angle<double> angle_90(Degrees(90.0));
  const vector::Vector3<double> pole_0(0.0, 0.0, 1.0);
  BOOST_CHECK_SMALL(vector::distance(pole_0, arc.angle_position(angle_90)),
                    CALCULATION_TOLERANCE);

  const auto end_arc{arc.end_arc(true)};
  BOOST_CHECK_EQUAL(Radians(0.02), end_arc.length());

  const auto end_arc_a{end_arc.a()};
  BOOST_CHECK_EQUAL(end_arc_a, arc.perp_position(arc.b(), Radians(0.01)));
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_north_and_south_poles) {
  const LatLong<double> north_pole(Degrees(90.0), Degrees(0.0));
  const LatLong<double> south_pole(Degrees(-90.0), Degrees(0.0));

  const auto [azimuth_s, distance_s]{
      calculate_azimuth_and_distance(south_pole, north_pole)};
  BOOST_CHECK_EQUAL(Degrees(0.0), azimuth_s.to_degrees());
  BOOST_CHECK_EQUAL(Radians(trig::PI<double>), distance_s);

  const auto [azimuth_n, distance_n]{
      calculate_azimuth_and_distance(north_pole, south_pole)};
  BOOST_CHECK_EQUAL(Degrees(180.0), azimuth_n.to_degrees());
  BOOST_CHECK_EQUAL(Radians(trig::PI<double>), distance_n);

  // 90 degrees East on the equator
  const LatLong<double> e_eq(Degrees(0.0), Degrees(90.0));

  const Arc<double> arc_e(north_pole, e_eq);
  BOOST_CHECK_SMALL(vector::latitude(arc_e.b()).to_degrees().v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(e_eq.lon().v(),
                    vector::longitude(arc_e.b()).to_degrees().v(),
                    CALCULATION_TOLERANCE);

  // 90 degrees West on the equator
  const LatLong<double> w_eq(Degrees(0.0), Degrees(-90.0));

  const Arc<double> arc_w(north_pole, w_eq);
  BOOST_CHECK_SMALL(vector::latitude(arc_w.b()).to_degrees().v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(w_eq.lon().v(),
                    vector::longitude(arc_w.b()).to_degrees().v(),
                    CALCULATION_TOLERANCE);

  const Arc<double> arc_s(south_pole, w_eq);
  BOOST_CHECK_SMALL(vector::latitude(arc_s.b()).to_degrees().v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(w_eq.lon().v(),
                    vector::longitude(arc_s.b()).to_degrees().v(),
                    CALCULATION_TOLERANCE);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_arc_atd_and_xtd) {
  // Greenwich equator
  const LatLong<double> g_eq(Degrees(0.0), Degrees(0.0));

  // 90 degrees East on the equator
  const LatLong<double> e_eq(Degrees(0.0), Degrees(90.0));

  const Arc<double> arc(g_eq, e_eq);
  BOOST_CHECK(arc.is_valid());

  const auto start_arc{arc.end_arc(false)};
  BOOST_CHECK_EQUAL(Radians(0.0), start_arc.length());

  const auto start_arc_a{start_arc.a()};
  BOOST_CHECK_EQUAL(arc.a(), start_arc_a);

  const Degrees<double> longitude(1.0);

  // Test across track distance
  // Accuracy drops off outside of this range
  for (auto lat{-83}; lat < 84; ++lat) {
    const Degrees<double> latitude(lat);
    const LatLong<double> latlong(latitude, longitude);
    const auto point{latlong.to_point()};

    const double expected{trig::deg2rad(latitude.v())};
    const auto [atd, xtd]{arc.calculate_atd_and_xtd(point)};
    BOOST_CHECK_CLOSE(trig::deg2rad(1.0), atd.v(), CALCULATION_TOLERANCE);
    BOOST_CHECK_CLOSE(expected, xtd.v(), 2 * CALCULATION_TOLERANCE);

    const auto d{arc.shortest_distance(point)};
    BOOST_CHECK_CLOSE(std::abs(expected), d.v(), 2 * CALCULATION_TOLERANCE);
  }

  auto point = g_eq.to_point();
  auto d = arc.shortest_distance(point);
  BOOST_CHECK_EQUAL(Radians(0.0), d);

  point = e_eq.to_point();
  d = arc.shortest_distance(point);
  BOOST_CHECK_EQUAL(Radians(0.0), d);

  const LatLong<double> latlong(Degrees(0.0), Degrees(-1.0));
  point = latlong.to_point();
  d = arc.shortest_distance(point);
  BOOST_CHECK_CLOSE(trig::deg2rad(1.0), d.v(), CALCULATION_TOLERANCE);

  point = -point;
  d = arc.shortest_distance(point);
  BOOST_CHECK_CLOSE(trig::deg2rad(89.0), d.v(), CALCULATION_TOLERANCE);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_arc_intersection_point) {
  // Karney's example:
  // Istanbul, Washington, Reyjavik and Accra
  // from:
  // <https://sourceforge.net/p/geographiclib/discussion/1026621/thread/21aaff9f/#fe0a>
  const LatLong<double> istanbul(Degrees(42.0), Degrees(29.0));
  const LatLong<double> washington(Degrees(39.0), Degrees(-77.0));
  const LatLong<double> reyjavik(Degrees(64.0), Degrees(-22.0));
  const LatLong<double> accra(Degrees(6.0), Degrees(0.0));

  const Arc<double> arc1(istanbul, washington);
  const Arc<double> arc2(reyjavik, accra);

  const auto intersection_point_1{calculate_intersection_point(arc1, arc2)};
  BOOST_CHECK(intersection_point_1.has_value());

  const LatLong<double> latlong_1{intersection_point_1.value()};
  // Geodesic intersection latitude is 54.7170296089477
  BOOST_CHECK_CLOSE(54.72, latlong_1.lat().v(), 100 * 0.05);
  // Geodesic intersection longitude is -14.56385574430775
  BOOST_CHECK_CLOSE(-14.56, latlong_1.lon().v(), 100 * 0.02);

  // Switch arcs
  const auto intersection_point_2{calculate_intersection_point(arc2, arc1)};
  BOOST_CHECK(intersection_point_2.has_value());

  const LatLong<double> latlong_2{intersection_point_2.value()};
  // Geodesic intersection latitude is 54.7170296089477
  BOOST_CHECK_CLOSE(54.72, latlong_2.lat().v(), 100 * 0.05);
  // Geodesic intersection longitude is -14.56385574430775
  BOOST_CHECK_CLOSE(-14.56, latlong_2.lon().v(), 100 * 0.02);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_arc_intersection_same_great_circles) {
  const LatLong<double> south_pole_1(Degrees(-88.0), Degrees(-180.0));
  const LatLong<double> south_pole_2(Degrees(-87.0), Degrees(0.0));

  const Arc<double> arc1(south_pole_1, south_pole_2);

  const auto intersection_lengths{calculate_intersection_distances(arc1, arc1)};
  BOOST_CHECK_EQUAL(Radians(0.0), std::get<0>(intersection_lengths));
  BOOST_CHECK_EQUAL(Radians(0.0), std::get<1>(intersection_lengths));

  const auto intersection_point_1{calculate_intersection_point(arc1, arc1)};
  BOOST_CHECK(intersection_point_1.has_value());

  const LatLong<double> south_pole_3(Degrees(-85.0), Degrees(0.0));
  const LatLong<double> south_pole_4(Degrees(-86.0), Degrees(0.0));
  const Arc<double> arc2(south_pole_3, south_pole_4);
  const auto intersection_point_2{calculate_intersection_point(arc1, arc2)};
  BOOST_CHECK(!intersection_point_2.has_value());
}
//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()
//////////////////////////////////////////////////////////////////////////////
