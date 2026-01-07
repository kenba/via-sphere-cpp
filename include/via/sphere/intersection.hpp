#pragma once

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
/// @file
/// @brief The via::vector::intersection namespace.
//////////////////////////////////////////////////////////////////////////////
#include "vector.hpp"
#include <cmath>

namespace via {
namespace vector {
namespace intersection {

/// Calculate an intersection point between the poles of two Great Circles.
/// See:
/// <http://www.movable-type.co.uk/scripts/latlong-vectors.html#intersection>
/// @param pole1, pole2 the poles.
/// @param min_sq_value minimum square of a vector length to normalize
///
/// @return an intersection point or None if the poles represent coincident
/// Great Circles.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto calculate_intersection(const Vector3<T> &pole1,
                                      const Vector3<T> &pole2,
                                      const T min_sq_value) noexcept
    -> std::optional<Vector3<T>> {
  return normalise(pole1.cross(pole2), min_sq_value);
}

/// Determine whether the antipodal point is closer to the centroid of the
/// `Arc`s.
///
/// @param point a great-circle intersection point.
/// @param centroid the centroid (geometric mean) of the `Arc`s mid points.
///
/// @return true if the antipodal intersection is closer to the `centroid`
/// of the `Arc`s otherwise returns false.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto use_antipodal_point(const Vector3<T> &point,
                                   const Vector3<T> &centroid) noexcept
    -> bool {
  return sq_distance<T>(centroid, -point) < sq_distance<T>(centroid, point);
}

/// Return the closer intersection point to the centroid of the `Arc`s.
///
/// @param point a great-circle intersection point.
/// @param centroid the centroid (geometric mean) of the `Arc`s mid points.
///
/// @return the antipodal point if it is closer to the `centroid`,
/// otherwise returns the point.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto closest_intersection_point(const Vector3<T> &point,
                                          const Vector3<T> &centroid) noexcept
    -> Vector3<T> {
  return use_antipodal_point(point, centroid) ? -point : point;
}

/// Determine the reference point of a pair of arcs.
/// I.e. the closest intersection point if they intersect or the
/// centroid normalized to lie on the unit sphere if they don't.
///
/// @param mid_point_0, mid_point_1 the mid points of the `Arc`s.
/// @param pole_0, `pole_1` the poles of the `Arc` great circles.
///
/// @return the closest intersection point or normalized centroid and the
/// sine of the angle between the arcs, zero if the arcs are coincident.
/// And the absolute relative angle at the intersection point or centroid.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto calculate_reference_point_and_angle(
    const Vector3<T> &mid_point_0, const Vector3<T> &pole_0,
    const Vector3<T> &mid_point_1, const Vector3<T> &pole_1) noexcept
    -> std::tuple<Vector3<T>, Angle<T>> {
  const Vector3<T> centroid{(mid_point_0 + mid_point_1)};
  const Vector3<T> point{pole_0.cross(pole_1)};
  const auto p{normalise(point, MIN_SQ_NORM<T>)};
  if (p.has_value()) {
    // the great circles intersect
    const auto x{closest_intersection_point(p.value(), centroid)};
    return {x, Angle<T>::from_y_x(point.norm(), pole_0.dot(pole_1))};
  } else {
    // the great circles are coincident
    const Vector3<T> c{normalise_centroid(centroid, mid_point_0, pole_0)};
    const Angle<T> angle{
        std::signbit(pole_0.dot(pole_1)) ? Angle<T>().opposite() : Angle<T>()};

    return {c, angle};
  }
}

/// Calculate signed great circle distances from two arc mid points to their
/// closest intersection point or normalized centroid if the arcs are on
/// coincident great circles.
///
/// @param mid_point_0, mid_point_1 the mid points of the arcs.
/// @param pole_0, pole_1 the poles of the arc great circles.
///
/// @return the signed great circle distances of the closest intersection
/// point or centroid  from the arc mid points in `Radians`,
/// and the relative angle between the arc great circles.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto calculate_arc_reference_distances_and_angle(
    const Vector3<T> &mid_point_0, const Vector3<T> &pole_0,
    const Vector3<T> &mid_point_1, const Vector3<T> &pole_1) noexcept
    -> std::tuple<Radians<T>, Radians<T>, Angle<T>> {
  const auto [point, angle]{calculate_reference_point_and_angle(
      mid_point_0, pole_0, mid_point_1, pole_1)};
  const Radians<T> distance_0{
      calculate_great_circle_atd(mid_point_0, pole_0, point)};
  const Radians<T> distance_1{
      calculate_great_circle_atd(mid_point_1, pole_1, point)};
  return {distance_0, distance_1, angle};
}

} // namespace intersection
} // namespace vector
} // namespace via
