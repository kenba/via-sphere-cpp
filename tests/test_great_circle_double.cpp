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
#include "via/sphere/great_circle.hpp"
#include <boost/test/unit_test.hpp>

using namespace via;
using namespace via::great_circle;

namespace {
constexpr auto EPSILON{std::numeric_limits<double>::epsilon()};
constexpr auto CALCULATION_TOLERANCE{101 * EPSILON};
} // namespace

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(Test_great_circle_double)

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_distance_conversion_functions) {
  BOOST_CHECK_EQUAL(trig::PI<double>, e2gc_distance(2.0).v());
  BOOST_CHECK_CLOSE(trig::PI_2<double>, e2gc_distance(std::sqrt(2.0)).v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(std::sqrt(2.0), gc2e_distance(Radians(trig::PI_2<double>)),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(2.0, gc2e_distance(Radians(trig::PI<double>)));
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_distance_and_azimuth_functions) {
  const Angle<double> angle_0;
  const auto angle_45{Angle<double>::from_y_x(1.0, 1.0)};

  auto gc_distance{calculate_gc_distance(angle_0, angle_45, -angle_45)};
  auto haversine_distance{
      calculate_haversine_distance(angle_0, angle_45, -angle_45, -angle_45)};
  BOOST_CHECK_CLOSE(gc_distance.v(), haversine_distance.v(),
                    CALCULATION_TOLERANCE);

  auto gc_azimuth{calculate_gc_azimuth(angle_0, angle_45, -angle_45)};
  BOOST_CHECK_CLOSE(-35.264389682754654, gc_azimuth.to_degrees().v(),
                    CALCULATION_TOLERANCE);

  // Same point
  gc_distance = calculate_gc_distance(angle_45, angle_45, angle_0);
  BOOST_CHECK_EQUAL(Radians(0.0), gc_distance);
  haversine_distance =
      calculate_haversine_distance(angle_45, angle_45, angle_0, angle_0);
  BOOST_CHECK_EQUAL(Radians(0.0), haversine_distance);
  gc_azimuth = calculate_gc_azimuth(angle_45, angle_45, angle_0);
  BOOST_CHECK_EQUAL(Degrees(0.0), gc_azimuth.to_degrees());
}
//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_north_and_south_pole_azimuths) {
  const auto angle_90{Angle<double>::from_y_x(1.0, 0.0)};
  const Angle angle_m90{-angle_90};

  const Angle<double> angle_0;
  const auto angle_45{Angle<double>::from_y_x(1.0, 1.0)};
  const Angle angle_180{angle_0.opposite()};

  // From South Pole
  auto gc_azimuth{calculate_gc_azimuth(angle_m90, angle_45, angle_0)};
  BOOST_CHECK_EQUAL(angle_0, gc_azimuth);

  gc_azimuth = calculate_gc_azimuth(angle_90, angle_45, angle_0);
  BOOST_CHECK_EQUAL(angle_180, gc_azimuth);
}
//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()
//////////////////////////////////////////////////////////////////////////////
