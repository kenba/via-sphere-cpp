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
#include "great_circle.hpp"
#include <Eigen/Geometry>
#include <optional>
#include <via/angle.hpp>

namespace via {
namespace vector {

/// The minimum value of the square of distance.
template <typename T>
  requires std::floating_point<T>
constexpr T MIN_SQ_DISTANCE{great_circle::MIN_VALUE<T> *
                            great_circle::MIN_VALUE<T>};

/// A `Vector3` is an [Eigen](https://eigen.tuxfamily.org/) `Vector3`.
template <typename T> using Vector3 = Eigen::Vector3<T>;

/// 2d vector perp product function: a x b.
/// @param a, b the two vectors.
///
/// @return the 2D perp product of the two vector x and y values.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto perp_product(const Vector3<T> &a, const Vector3<T> &b) noexcept
    -> T {
  return a(0) * b(1) - a(1) * b(0);
}

template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto dot2d(const Vector3<T> &a, const Vector3<T> &b) noexcept -> T {
  return a(0) * b(0) + a(1) * b(1);
}

/// Convert a latitude and longitude to a point on the unit sphere.
/// @pre |lat| <= 90.0 degrees.
/// @param lat the latitude.
/// @param lon the longitude.
///
/// @return a `Vector3` of the point on the unit sphere.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto to_point(const Angle<T> lat, const Angle<T> lon) noexcept
    -> Vector3<T> {
  // clang-format off
  return Vector3<T>(lat.cos().v() * lon.cos().v(),
                    lat.cos().v() * lon.sin().v(),
                    lat.sin().v());
  // clang-format on
}

/// Calculate the latitude of a point.
/// @param a the point.
///
/// @return the latitude of the point
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto latitude(const Vector3<T> &a) noexcept -> Angle<T> {
  const auto sin_a = trig::UnitNegRange<T>(a(2));
  return Angle(sin_a, trig::swap_sin_cos(sin_a));
}

/// Calculate the longitude of a point.
/// @param a the point.
///
/// @return the longitude of the point
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto longitude(const Vector3<T> &a) noexcept -> Angle<T> {
  return Angle(a(1), a(0));
}

/// Determine whether a `Vector3` is a unit vector.
/// @param a the vector.
///
/// @return true if `a` is a unit vector, false otherwise.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto is_unit(const Vector3<T> &a) noexcept -> bool {
  constexpr T MIN_POINT_SQ_LENGTH{1 - 12 * std::numeric_limits<T>::epsilon()};
  constexpr T MAX_POINT_SQ_LENGTH{1 + 12 * std::numeric_limits<T>::epsilon()};

  const auto sq_length{a.squaredNorm()};
  return MIN_POINT_SQ_LENGTH <= sq_length && sq_length <= MAX_POINT_SQ_LENGTH;
}

/// Normalize a vector to lie on the surface of the unit sphere.
/// Note: this function returns an `Option` so uses the British spelling of
/// `normalise` to differentiate it from the standard `normalize` function.
/// @param a the `Vector3`
///
/// @return the nomalized point or None if the vector is too small to normalize.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto normalise(const Vector3<T> &a) noexcept
    -> std::optional<Vector3<T>> {
  constexpr T MIN_LENGTH{16384 * std::numeric_limits<T>::epsilon()};
  constexpr T MIN_NORM{MIN_LENGTH * MIN_LENGTH};

  if (a.squaredNorm() >= MIN_NORM)
    return a.normalized();

  return std::nullopt;
}

/// Calculate the square of the Euclidean distance between two points.
/// Note: points do NOT need to be valid Points.
/// @post for unit vectors: result <= 4
/// @param a, b the points.
///
/// @return the square of the Euclidean distance between the points.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto sq_distance(const Vector3<T> &a, const Vector3<T> &b) noexcept
    -> T {
  return (b - a).squaredNorm();
}

/// Calculate the shortest (Euclidean) distance between two Points.
/// @post for unit vectors: result <= 2
/// @param a, b the points.
///
/// @return the shortest (Euclidean) distance between the points.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto distance(const Vector3<T> &a, const Vector3<T> &b) noexcept
    -> T {
  return (b - a).norm();
}

/// Determine whether two `Vector3`s are orthogonal (perpendicular).
/// @param a, b the `Vector3`s.
///
/// @return true if a and b are orthogonal, false otherwise.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto are_orthogonal(const Vector3<T> &a, const Vector3<T> &b) noexcept
    -> bool {
  constexpr T MAX_LENGTH{4 * std::numeric_limits<T>::epsilon()};

  const auto length{a.dot(b)};
  return -MAX_LENGTH <= length && length <= MAX_LENGTH;
}

/// Calculate the relative longitude of point a from point b.
/// @param a, b the points.
///
/// @return the relative longitude of point a from point b,
/// negative if a is West of b, positive otherwise.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto delta_longitude(const Vector3<T> &a,
                               const Vector3<T> &b) noexcept -> Angle<T> {
  return Angle<T>(perp_product(b, a), dot2d(b, a));
}

/// Determine whether point a is West of point b.
/// It calculates and compares the perp product of the two points.
/// @param a, b the points.
///
/// @return true if a is West of b, false otherwise.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto is_west_of(const Vector3<T> &a, const Vector3<T> &b) noexcept
    -> bool {
  return perp_product(b, a) <= -std::numeric_limits<T>::epsilon();
}

/// Calculate the right hand pole vector of a Great Circle from an initial
/// position and an azimuth.
/// See: <http://www.movable-type.co.uk/scripts/latlong-vectors.html#distance>
/// @param lat - start point Latitude.
/// @param lon - start point Longitude.
/// @param azi - start point azimuth.
///
/// @return the right hand pole vector of the great circle.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto calculate_pole(const Angle<T> lat, const Angle<T> lon,
                              const Angle<T> azi) noexcept -> Vector3<T> {
  // clang-format off
  const auto x{trig::UnitNegRange<T>::clamp(
      lon.sin().v() * azi.cos().v() - lat.sin().v() * lon.cos().v() * azi.sin().v()
  )};
  const auto y{trig::UnitNegRange<T>::clamp(
      T(0) - lon.cos().v() * azi.cos().v() - lat.sin().v() * lon.sin().v() * azi.sin().v()
  )};
  const auto z{trig::UnitNegRange<T>::clamp(
      lat.cos().v() * azi.sin().v()
  )};
  // clang-format on
  return Vector3<T>(x.v(), y.v(), z.v());
}

/// Calculate the azimuth at a point on the Great Circle defined by pole.
/// @param point - the point.
/// @param pole - the right hand pole of the Great Circle.
///
/// @return the azimuth at the point on the great circle.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto calculate_azimuth(const Vector3<T> &point,
                                 const Vector3<T> &pole) noexcept -> Angle<T> {
  constexpr auto MAX_LAT{T(1) - great_circle::MIN_VALUE<T>};

  const auto sin_lat{point(2)};
  // if the point is close to the North or South poles, azimuth is 180 or 0.
  if (MAX_LAT <= std::abs(sin_lat)) {
    return std::signbit(sin_lat) ? Angle<T>() : Angle<T>().opposite();
  }

  return Angle<T>(pole(2), perp_product(pole, point));
}

/// Calculate the direction vector along a Great Circle from an initial
/// position and an azimuth.
/// See: Panou and Korakitis equations: 30, 31, & 32a
/// <https://arxiv.org/abs/1811.03513>
/// @param lat - start point Latitude.
/// @param lon - start point Longitude.
/// @param azi - start point azimuth.
///
/// @return the direction vector at the point on the great circle.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto calculate_direction(const Angle<T> lat, const Angle<T> lon,
                                   const Angle<T> azi) noexcept -> Vector3<T> {
  // clang-format off
  const auto x{trig::UnitNegRange<T>::clamp(
      T(0) - lat.sin().v() * lon.cos().v() * azi.cos().v() - lon.sin().v() * azi.sin().v()
  )};
  const auto y{trig::UnitNegRange<T>::clamp(
      T(0) - lat.sin().v() * lon.sin().v() * azi.cos().v() + lon.cos().v() * azi.sin().v()
  )};
  const auto z{trig::UnitNegRange<T>::clamp(
      lat.cos().v() * azi.cos().v()
  )};
  // clang-format on
  return Vector3<T>(x.v(), y.v(), z.v());
}

/// Calculate the direction vector of a Great Circle arc.
/// @param a - the start point.
/// @param pole - the pole of a Great Circle.
///
/// @return the direction vector at the point on the great circle.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto direction(const Vector3<T> &a, const Vector3<T> &pole) noexcept
    -> Vector3<T> {
  return pole.cross(a);
}

/// Calculate the position of a point along a Great Circle arc.
/// @param a - the start point.
/// @param dir - the direction vector of a Great Circle at a.
/// @param distance - the a Great Circle as an Angle.
///
/// @return the position vector at the point on the great circle.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto position(const Vector3<T> &a, const Vector3<T> &dir,
                        const Angle<T> distance) noexcept -> Vector3<T> {
  return distance.cos().v() * a + distance.sin().v() * dir;
}

/// Calculate the direction vector of a Great Circle rotated by angle.
/// @param dir the direction vector of a Great Circle arc.
/// @param pole the pole of a Great Circle.
/// @param angle the angle to rotate the direction vector by.
///
/// @return the direction vector at the point on the great circle
/// rotated by angle.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto rotate(const Vector3<T> &dir, const Vector3<T> &pole,
                      const Angle<T> angle) noexcept -> Vector3<T> {
  return position(dir, pole, angle);
}

/// Calculate the position of a point rotated by angle at radius.
/// @param a the start point.
/// @param pole the pole of a Great Circle.
/// @param angle the angle to rotate the direction vector by.
/// @param radius the radius from the start point.
///
/// @return the position vector at angle and radius from the start point.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto rotate_position(const Vector3<T> &a, const Vector3<T> &pole,
                               const Angle<T> angle,
                               const Angle<T> radius) noexcept -> Vector3<T> {
  return position(a, rotate(direction(a, pole), pole, angle), radius);
}

/// The sine of the across track distance of a point relative to a Great Circle
/// pole. It is simply the dot product of the pole and the point: pole . point
/// @param pole the Great Circle pole.
/// @param point the point.
///
/// @return the sine of the across track distance of point relative to the pole.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto sin_xtd(const Vector3<T> &pole, const Vector3<T> &point) noexcept
    -> trig::UnitNegRange<T> {
  return trig::UnitNegRange<T>::clamp(pole.dot(point));
}

/// The across track distance of a point relative to a Great Circle pole.
/// @param pole the Great Circle pole.
/// @param point the point.
///
/// @return the across track distance of point relative to pole, in `Radians`.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto cross_track_distance(const Vector3<T> &pole,
                                    const Vector3<T> &point) noexcept
    -> Radians<T> {
  const auto sin_d{sin_xtd(pole, point)};
  return Radians<T>((std::abs(sin_d.v()) < std::numeric_limits<T>::epsilon())
                        ? 0
                        : std::asin(sin_d.v()));
}

/// The square of the Euclidean cross track distance of a point relative to a
/// Great Circle pole.
/// @param pole the Great Circle pole.
/// @param point the point.
///
/// @return the square of the euclidean distance of point relative to pole.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto sq_cross_track_distance(const Vector3<T> &pole,
                                       const Vector3<T> &point) noexcept -> T {
  const auto sin_d{sin_xtd(pole, point)};
  return (std::abs(sin_d.v()) < std::numeric_limits<T>::epsilon())
             ? T(0)
             : 2 * (T(1) - trig::swap_sin_cos(sin_d).v());
}

/// Calculate the closest point on a plane to the given point.
/// See: [Closest Point on
/// Plane](https://gdbooks.gitbooks.io/3dcollisions/content/Chapter1/closest_point_on_plane.html)
/// @param pole the Great Circle pole (aka normal) of the plane.
/// @param point the point.
///
/// @return the closest point on a plane to the given point.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto calculate_point_on_plane(const Vector3<T> &pole,
                                        const Vector3<T> &point) noexcept
    -> Vector3<T> {
  const auto t{sin_xtd(pole, point)};
  return point - t.v() * pole;
}

/// The sine of the along track distance of a point along a Great Circle arc.
/// It is the triple product of the pole, a and the point:
/// (pole X a) . point = pole . (a X point)
/// @param a the start point of the Great Circle arc.
/// @param pole the pole of the Great Circle arc.
/// @param point the point.
///
/// @return the sine of the along track distance of point relative to the start
/// of a great circle arc.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto sin_atd(const Vector3<T> &a, const Vector3<T> &pole,
                       const Vector3<T> &point) noexcept
    -> trig::UnitNegRange<T> {
  return trig::UnitNegRange<T>::clamp(pole.cross(a).dot(point));
}

/// Calculate the relative distance of two points on a Great Circle arc.
/// @pre both points must be on the Great Circle defined by `pole`.
/// @param a the start point of the Great Circle arc.
/// @param pole the pole of the Great Circle arc.
/// @param point the point.
///
/// @return the Great Circle along track distance in `Radians`.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto calculate_great_circle_atd(const Vector3<T> &a,
                                          const Vector3<T> &pole,
                                          const Vector3<T> &point) noexcept
    -> Radians<T> {
  const auto sq_atd{sq_distance(a, point)};
  return Radians<T>(
      (sq_atd < MIN_SQ_DISTANCE<T>)
          ? 0
          : std::copysign(great_circle::e2gc_distance(std::sqrt(sq_atd)).v(),
                          sin_atd(a, pole, point).v()));
}

/// The Great Circle distance of a point along the arc relative to a,
/// (+ve) ahead of a, (-ve) behind a.
/// @param a the start point of the Great Circle arc.
/// @param pole the pole of the Great Circle arc.
/// @param point the point.
///
/// @return the along track distance of point relative to the start of a great
/// circle arc.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto along_track_distance(const Vector3<T> &a, const Vector3<T> &pole,
                                    const Vector3<T> &point) noexcept
    -> Radians<T> {
  const auto plane_point{normalise(calculate_point_on_plane(pole, point))};
  return plane_point.has_value()
             ? calculate_great_circle_atd(a, pole, plane_point.value())
             : Radians<T>(0);
}

/// Calculate the square of the Euclidean along track distance of a point
/// from the start of an Arc.
/// It is calculated using the closest point on the plane to the point.
/// @param a the start point of the Great Circle arc.
/// @param pole the pole of the Great Circle arc.
/// @param point the point.
///
/// @return the square of the Euclidean along track distance
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto sq_along_track_distance(const Vector3<T> &a,
                                       const Vector3<T> &pole,
                                       const Vector3<T> &point) noexcept -> T {
  const auto plane_point{normalise<T>(calculate_point_on_plane(pole, point))};
  return plane_point.has_value() ? sq_distance(a, plane_point.value()) : T(0);
}

template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto calculate_atd_and_xtd(const Vector3<T> &a,
                                     const Vector3<T> &pole,
                                     const Vector3<T> &p) noexcept
    -> std::tuple<Radians<T>, Radians<T>> {
  auto atd{Radians<T>(0)};
  auto xtd{Radians<T>(0)};

  const auto sq_d{sq_distance(a, p)};
  if (sq_d >= MIN_SQ_DISTANCE<T>) {
    const auto sine_xtd{sin_xtd(pole, p).v()};
    if (std::abs(sine_xtd) >= std::numeric_limits<T>::epsilon())
      xtd = Radians<T>(std::asin(sine_xtd));

    const auto plane_point{normalise<T>(p - sine_xtd * pole)};
    atd = plane_point.has_value()
              ? calculate_great_circle_atd(a, pole, plane_point.value())
              : Radians<T>(0);
  }

  return {atd, xtd};
}
} // namespace vector
} // namespace via
