#pragma once

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

/// Calculate the great circle distances to an intersection point from the
/// start points of a pair of great circle arcs, on different great circles.
/// @param a1, a2 the start points of the great circle arcs
/// @param pole1, pole2 the poles of the great circle arcs
/// @param c the intersection point
///
/// @return a pair of great circle distances along the arcs to the
/// intersection point in `Radians`.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto
calculate_intersection_distances(const Vector3<T> &a1, const Vector3<T> &pole1,
                                 const Vector3<T> &a2, const Vector3<T> &pole2,
                                 const Vector3<T> &c) noexcept
    -> std::tuple<Radians<T>, Radians<T>> {
  return {calculate_great_circle_atd(a1, pole1, c),
          calculate_great_circle_atd(a2, pole2, c)};
}

/// Whether an along track distance is within an `Arc` length including
/// tolerance.
/// @param distance the along track distance from the start of the `Arc`.
/// @param length the length of the `Arc`.
/// @param tolerance the distance tolerance.
///
/// @return true if the along track distance is within the length including
/// tolerance, false otherwise.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto is_alongside(const Radians<T> distance, const Radians<T> length,
                            const Radians<T> tolerance) noexcept -> bool {
  return (-tolerance.v() <= distance.v()) &&
         (distance.v() <= length.v() + tolerance.v());
}

/// Whether an intersection point is within an `Arc`.
/// @param distance the along track distance to the point from the start of the
/// `Arc`.
/// @param length the length of the `Arc`.
///
/// @return true if the intersection point is within the `Arc`, false otherwise.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto is_within(const T distance, const T length) noexcept -> bool {
  return (-std::numeric_limits<T>::epsilon() <= distance) &&
         (distance <=
          length + (std::numeric_limits<T>::epsilon() * (T(1) + length)));
}

/// Calculate the great-circle distances along a pair of `Arc`s on coincident
/// Great Circles to their closest (reference) points.
/// @param gc_d the great-circle distance between the arc start points.
/// @param reciprocal whether the arcs are in reciprocal directions.
/// @param arc1_length, arc2_length the `Arc` lengths in `Radians`.
///
/// @return the distances along the first `Arc` and second `Arc` to their
/// closest (reference) points in `Radians`.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto
calculate_coincident_arc_distances(const Radians<T> gc_d, const bool reciprocal,
                                   const Radians<T> arc1_length,
                                   const Radians<T> arc2_length) noexcept
    -> std::tuple<Radians<T>, Radians<T>> {
  if (reciprocal) {
    // if the arcs intersect
    if (is_alongside(gc_d,
                     Radians<T>(std::max(arc1_length.v(), arc2_length.v())),
                     Radians<T>(4 * std::numeric_limits<T>::epsilon()))) {
      // The start of the first `Arc` is within the second `Arc`
      if (gc_d.v() <= arc2_length.v())
        // The start of the first `Arc` is within the second `Arc`
        return {Radians<T>(0), gc_d.clamp(arc2_length)};
      else
        // The start of the second `Arc` is within the first `Arc`
        return {gc_d.clamp(arc1_length), Radians<T>(0)};
    } else {
      const auto abs_d{gc_d.abs()};

      // The distance between the `Arc` b ends
      const auto b_d{abs_d.v() - arc1_length.v() - arc2_length.v()};
      // The distance between the `Arc` b ends around the Great Circle
      const auto b_gc_d{(gc_d.v() > T(0)) ? b_d : trig::TAU<T> + gc_d.v()};

      if (b_gc_d < abs_d.v())
        // The end of the second `Arc` is beyond the end of first `Arc`
        return {Radians(b_gc_d) + arc1_length, arc2_length};
      else
        // The start of the second `Arc` is before the start of first `Arc`
        return {-abs_d, Radians<T>(0)};
    }
  } else {
    // The distance to the start of arc2 from the end of arc1
    const auto b1a2{(gc_d.v() > T(0))
                        ? gc_d.v() - arc1_length.v()
                        : trig::TAU<T> + gc_d.v() - arc1_length.v()};
    // The distance to the start of arc1 from the end of arc2
    const auto b2a1{(gc_d.v() > T(0))
                        ? trig::TAU<T> + gc_d.v() - arc2_length.v()
                        : -gc_d.v() - arc2_length.v()};

    if (b2a1 < b1a2)
      // The start of the first arc is within the second arc
      return {Radians<T>(0), Radians<T>(b2a1 + arc2_length.v())};
    else
      // The start of the second arc relative to the start of first arc.
      return {Radians<T>(b1a2 + arc1_length.v()), Radians<T>(0)};
  }
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

/// Calculate the great-circle distances along a pair of arcs to their
/// closest intersection point or their coincident arc distances if the
/// `Arc`s are on coincident Great Circles.
/// @param a1, a2 the `Arc` start points.
/// @param pole1, pole2 the `Arc` poles.
/// @param length1, length2 the `Arc` lengths.
/// @param centroid the centroid (geometric mean) of the `Arc`s mid points.
///
/// @return the distances along the first arc and second arc to the intersection
/// point or to their coincident arc distances if the arcs do not intersect.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto calculate_intersection_point_distances(
    const Vector3<T> &a1, const Vector3<T> &pole1, const Radians<T> length1,
    const Vector3<T> &a2, const Vector3<T> &pole2, const Radians<T> length2,
    const Vector3<T> &centroid) noexcept -> std::tuple<Radians<T>, Radians<T>> {
  // Calculate the square of the Euclidean distance between the start points.
  const auto sq_d{sq_distance<T>(a1, a2)};
  if (sq_d < MIN_SQ_DISTANCE<T>)
    return {Radians<T>(0), Radians<T>(0)};

  const auto c{calculate_intersection(pole1, pole2, MIN_SQ_NORM<T>)};
  if (c.has_value()) {
    const auto d{use_antipodal_point(c.value(), centroid) ? -c.value()
                                                          : c.value()};
    return calculate_intersection_distances<T>(a1, pole1, a2, pole2, d);
  } else {
    const auto gc_d{calculate_great_circle_atd(a1, pole1, a2)};
    const bool opposite{std::signbit(pole1.dot(pole2))};
    return calculate_coincident_arc_distances(gc_d, opposite, length1, length2);
  }
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
    // find the closest intersection point to both the arcs
    const Vector3<T> x{use_antipodal_point(p.value(), centroid) ? -p.value()
                                                                : p.value()};
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
