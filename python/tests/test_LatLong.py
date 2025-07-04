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
#  @file test_LatLong
#  @brief Contains unit tests for the via sphere LatLong class.

import pytest
import numpy as np
from numpy.testing import assert_almost_equal
from via_angle import Degrees, PI_3
from via_sphere import calculate_azimuth_and_distance, haversine_distance, \
    is_valid_latitude, is_valid_longitude, latitude, longitude, LatLong

def test_is_valid_latitude():
    # value < -90
    assert not is_valid_latitude(-90.0001)
    # value = -90
    assert is_valid_latitude(-90.0)
    # value = 90
    assert is_valid_latitude(90.0)
    # value > 90
    assert not is_valid_latitude(90.0001)

def test_is_valid_longitude():
    # value < -180
    assert not is_valid_longitude(-180.0001)
    # value = -180
    assert is_valid_longitude(-180.0)
    # value = 180
    assert is_valid_longitude(180.0)
    # value > 180
    assert not is_valid_longitude(180.0001)

def test_latlong_class():
    a = LatLong(Degrees(0.0), Degrees(90.0))
    assert a.is_valid()
    assert Degrees(0.0) == a.lat()
    assert Degrees(90.0) == a.lon()

    assert not a.is_south_of(a)
    assert not a.is_west_of(a)

    b = LatLong(Degrees(-10.0), Degrees(-91.0))
    assert b.is_south_of(a)
    assert b.is_west_of(a)

def test_latlong_vector3_conversion():
    a = LatLong(Degrees(0.0), Degrees(90.0))
    point = a.to_point()

    assert 0.0 == point[0]
    assert 1.0 == point[1]
    assert 0.0 == point[2]

    assert Degrees(0.0) == latitude(point).to_degrees()
    assert Degrees(90.0) == longitude(point).to_degrees()
    assert a == LatLong(point)

def test_great_circle_90n_0n_0e():
    a = LatLong(Degrees(90.0), Degrees(0.0))
    b = LatLong(Degrees(0.0), Degrees(0.0))
    azimuth, distance = calculate_azimuth_and_distance(a, b)

    assert_almost_equal(np.pi / 2.0, distance.v())
    assert Degrees(180.0) == azimuth.to_degrees()

    distance = haversine_distance(a, b)
    assert_almost_equal(np.pi / 2.0, distance.v())

def test_great_circle_90s_0n_50e():
    a = LatLong(Degrees(-90.0), Degrees(0.0))
    b = LatLong(Degrees(0.0), Degrees(50.0))
    azimuth, distance = calculate_azimuth_and_distance(a, b)

    assert_almost_equal(np.pi / 2.0, distance.v())
    assert Degrees(0.0) == azimuth.to_degrees()

    distance = haversine_distance(a, b)
    assert_almost_equal(np.pi / 2.0, distance.v())

def test_great_circle_0n_60e_0n_60w():
    a = LatLong(Degrees(0.0), Degrees(60.0))
    b = LatLong(Degrees(0.0), Degrees(-60.0))
    azimuth, distance = calculate_azimuth_and_distance(a, b)

    assert_almost_equal(2.0 * PI_3, distance.v())
    assert Degrees(-90.0) == azimuth.to_degrees()

    distance = haversine_distance(a, b)
    assert_almost_equal(2.0 * PI_3, distance.v())

if __name__ == '__main__':
    pytest.main()
