#!/usr/bin/env python

# Copyright (c) 2018-2024 Ken Barker
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
#  @file test_Arc
#  @brief Contains unit tests for the via sphere Arc class.

import pytest
import numpy as np
from numpy.testing import assert_almost_equal, assert_array_equal
from via_angle import Angle, Degrees, Radians, deg2rad
from via_sphere import Arc, LatLong, calculate_azimuth_and_distance, \
    calculate_intersection_point, distance, latitude, longitude

def test_arc():
    # Greenwich equator
    g_eq = LatLong(Degrees(0.0), Degrees(0.0))

    # 90 degrees East on the equator
    e_eq = LatLong(Degrees(0.0), Degrees(90.0))

    arc = Arc(g_eq, e_eq, Radians(0.01))
    assert arc.is_valid()
    assert Radians(0.01) == arc.half_width()

    assert_array_equal(g_eq.to_point(), arc.a())
    pole_0 = LatLong(Degrees(90.0), Degrees(0.0)).to_point()
    assert_array_equal(pole_0, arc.pole())
    assert_almost_equal(np.pi / 2.0, arc.length().v())
    assert Degrees(90.0) == arc.azimuth().to_degrees()
    b = e_eq.to_point()
    assert_almost_equal(0.0, distance(b, arc.b()))

    mid_point = arc.mid_point()
    assert 0.0 == mid_point[2]
    assert_almost_equal(45.0, longitude(mid_point).to_degrees().v())

    start_arc = arc.end_arc(False)
    assert Radians(0.02) == start_arc.length()

    start_arc_a = start_arc.a()
    assert_array_equal(start_arc_a, arc.perp_position(arc.a(), Radians(0.01)))

    angle_90 = Angle(Degrees(90.0))
    assert_almost_equal(0.0, distance(pole_0, arc.angle_position(angle_90)))

    end_arc = arc.end_arc(True)
    assert Radians(0.02) == end_arc.length()

    end_arc_a = end_arc.a()
    assert_array_equal(end_arc_a, arc.perp_position(arc.b(), Radians(0.01)))

def test_north_and_south_poles():
    north_pole = LatLong(Degrees(90.0), Degrees(0.0))
    south_pole = LatLong(Degrees(-90.0), Degrees(0.0))

    azimuth_s, distance_s = calculate_azimuth_and_distance(south_pole, north_pole)
    assert Degrees(0.0) == azimuth_s.to_degrees()
    assert Radians(np.pi) == distance_s

    azimuth_n, distance_n = calculate_azimuth_and_distance(north_pole, south_pole)
    assert Degrees(180.0) == azimuth_n.to_degrees()
    assert Radians(np.pi) == distance_n

    # 90 degrees East on the equator
    e_eq = LatLong(Degrees(0.0), Degrees(90.0))

    arc_e = Arc(north_pole, e_eq)
    assert_almost_equal(0.0, latitude(arc_e.b()).to_degrees().v())
    assert_almost_equal(e_eq.lon().v(), longitude(arc_e.b()).to_degrees().v())

    # 90 degrees West on the equator
    w_eq = LatLong(Degrees(0.0), Degrees(-90.0))

    arc_w = Arc(north_pole, w_eq)
    assert_almost_equal(0.0, latitude(arc_w.b()).to_degrees().v())
    assert_almost_equal(w_eq.lon().v(), longitude(arc_w.b()).to_degrees().v())

    arc_s = Arc(south_pole, w_eq)
    assert_almost_equal(0.0, latitude(arc_s.b()).to_degrees().v())
    assert_almost_equal(w_eq.lon().v(), longitude(arc_s.b()).to_degrees().v())

def test_arc_atd_and_xtd():
    # Greenwich equator
    g_eq = LatLong(Degrees(0.0), Degrees(0.0))

    # 90 degrees East on the equator
    e_eq = LatLong(Degrees(0.0), Degrees(90.0))

    arc = Arc(g_eq, e_eq)
    assert arc.is_valid()

    start_arc = arc.end_arc(False)
    assert Radians(0.0) == start_arc.length()

    start_arc_a = start_arc.a()
    assert_array_equal(arc.a(), start_arc_a)

    longitude = Degrees(1.0)

    for i in range(-89, 90):
        latitude = Degrees(i)
        latlong = LatLong(latitude, longitude)
        point = latlong.to_point()

        expected = deg2rad(latitude.v())
        atd, xtd = arc.calculate_atd_and_xtd(point)
        assert_almost_equal(deg2rad(1.0), atd.v())
        assert_almost_equal(expected, xtd.v())

def test_arc_intersection_point():
    # Karney's example:
    # Istanbul, Washington, Reyjavik and Accra
    # from:
    # <https://sourceforge.net/p/geographiclib/discussion/1026621/thread/21aaff9f/#fe0a>
    istanbul = LatLong(Degrees(42.0), Degrees(29.0))
    washington = LatLong(Degrees(39.0), Degrees(-77.0))
    reyjavik = LatLong(Degrees(64.0), Degrees(-22.0))
    accra = LatLong(Degrees(6.0), Degrees(0.0))

    arc1 = Arc(istanbul, washington)
    arc2 = Arc(reyjavik, accra)

    intersection_point_1 = calculate_intersection_point(arc1, arc2)
    latlong_1 = LatLong(intersection_point_1)
    assert_almost_equal(54.72, latlong_1.lat().v(), 0.05)
    assert_almost_equal(-14.56, latlong_1.lon().v(), 0.02)

if __name__ == '__main__':
    pytest.main()
