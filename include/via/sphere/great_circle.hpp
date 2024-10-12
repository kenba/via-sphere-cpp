#pragma once

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
#include <via/angle.hpp>

namespace via {
namespace great_circle {

/// The minimum value for angles and distances.
template <typename T>
  requires std::floating_point<T>
constexpr T MIN_VALUE{2 * std::numeric_limits<T>::epsilon()};

/// Calculate the Great Circle distance between two points from their
/// Latitude and Longitude differences, see:
/// [Haversine formula](https://en.wikipedia.org/wiki/Haversine_formula).
/// This function is less accurate than `calculate_gc_distance`.
/// @param a_lat start point Latitude.
/// @param b_lat  finish point Latitude.
/// @param delta_long Longitude difference between start and finish points.
/// @param delta_lat Latitude difference between start and finish points.
///
/// @return the Great Circle distance between the points in Radians.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto
calculate_haversine_distance(const Angle<T> a_lat, const Angle<T> b_lat,
                             const Angle<T> delta_long,
                             const Angle<T> delta_lat) noexcept -> Radians<T> {
  const auto haversine_lat{trig::sq_sine_half(delta_lat.cos())};
  const auto haversine_lon{trig::sq_sine_half(delta_long.cos())};

  const auto a = std::clamp<T>(
      (haversine_lat + a_lat.cos().v() * b_lat.cos().v() * haversine_lon), 0,
      1);
  return (a < MIN_VALUE<T>) ? Radians<T>(0)
                            : Radians(2 * std::asin(std::sqrt(a)));
}

/// Convert a Euclidean distance to a Great Circle distance (in `Radians`).
/// e should satisfy: 0 <= e <= 2, if not it is clamped into range.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto e2gc_distance(const T e) -> Radians<T> {
  return (e < MIN_VALUE<T>)
             ? Radians<T>(0)
             : Radians<T>(2 *
                          std::asin(trig::UnitNegRange<T>::clamp(e / 2).v()));
}

/// Convert a Great Circle distance (in radians) to a Euclidean distance.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto gc2e_distance(const Radians<T> gc) -> T {
  return 2 * std::sin(gc.v() / 2);
}

/// Calculate the square of the Euclidean distance (i.e. using Pythagoras)
/// between two points from their Latitudes and their Longitude difference.
/// @param a_lat start point Latitude.
/// @param b_lat  finish point Latitude.
/// @param delta_long Longitude difference between start and finish points.
///
/// @return the square of the Euclidean distance between the points.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto sq_euclidean_distance(const Angle<T> a_lat, const Angle<T> b_lat,
                                     const Angle<T> delta_long) noexcept -> T {
  const auto delta_x{b_lat.cos().v() * delta_long.cos().v() - a_lat.cos().v()};
  const auto delta_y{b_lat.cos().v() * delta_long.sin().v()};
  const auto delta_z{b_lat.sin().v() - a_lat.sin().v()};

  const auto result{delta_x * delta_x + delta_y * delta_y + delta_z * delta_z};
  return std::clamp<T>(result, 0, 4);
}

/// Calculate the Great Circle distance (angle from centre) between two points
/// from their Latitudes and their Longitude difference.
/// This function is more accurate than `calculate_haversine_distance`.
/// @param a_lat start point Latitude.
/// @param b_lat  finish point Latitude.
/// @param delta_long Longitude difference between start and finish points.
///
/// @return the Great Circle distance between the points in Radians.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto
calculate_gc_distance(const Angle<T> a_lat, const Angle<T> b_lat,
                      const Angle<T> delta_long) noexcept -> Radians<T> {
  return e2gc_distance(
      std::sqrt(sq_euclidean_distance(a_lat, b_lat, delta_long)));
}

/// Calculate the azimuth (bearing) along the great circle of point b from
/// point a from their Latitudes and their Longitude difference.
/// @param a_lat start point Latitude.
/// @param b_lat  finish point Latitude.
/// @param delta_long Longitude difference between start and finish points.
///
/// @return the Great Circle azimuth relative to North of point b from point a
/// as an Angle.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto
calculate_gc_azimuth(const Angle<T> a_lat, const Angle<T> b_lat,
                     const Angle<T> delta_long) noexcept -> Angle<T> {
  // if start point is North or South pole
  if (a_lat.cos().v() < MIN_VALUE<T>) {
    // azimuth is zero or 180 degrees
    return std::signbit(a_lat.sin().v()) ? Angle<T>() : Angle<T>().opposite();
  } else {
    const auto sin_azimuth{b_lat.cos().v() * delta_long.sin().v()};
    const auto temp{(a_lat.sin().v() * b_lat.cos().v() * delta_long.sin().v() *
                     delta_long.sin().v()) /
                    (1 + delta_long.cos().abs().v())};
    const auto cos_azimuth{(std::signbit(delta_long.cos().v()))
                               ? b_lat.sin().v() * a_lat.cos().v() +
                                     a_lat.sin().v() * b_lat.cos().v() - temp
                               : b_lat.sin().v() * a_lat.cos().v() -
                                     a_lat.sin().v() * b_lat.cos().v() + temp};
    return Angle(sin_azimuth, cos_azimuth);
  }
}
} // namespace great_circle
} // namespace via