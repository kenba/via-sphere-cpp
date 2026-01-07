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
/// @brief The via-sphere-cpp library header file.
//////////////////////////////////////////////////////////////////////////////
/// @mainpage via-sphere-cpp
///
/// The library uses a combination of spherical trigonometry and vector geometry
/// to perform [great-circle
/// navigation](https://en.wikipedia.org/wiki/Great-circle_navigation) on the
/// surface of a unit sphere.
///
//////////////////////////////////////////////////////////////////////////////
#include "via/sphere/intersection.hpp"
#include <cmath>
#include <via/angle.hpp>

namespace via {

/// Test whether a latitude in degrees is a valid latitude.
/// I.e. whether it lies in the range: -90.0 <= degrees <= 90.0
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto is_valid_latitude(const T degrees) -> T {
  return (T(-90) <= degrees) && (degrees <= T(90));
}

/// Test whether a longitude in degrees is a valid longitude.
/// I.e. whether it lies in the range: -180.0 <= degrees <= 180.0
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto is_valid_longitude(const T degrees) -> T {
  return (T(-180) <= degrees) && (degrees <= T(180));
}

/// LatLong represents a position by its latitude and longitude components.
/// @invariant latitude lies in the range: -90.0 <= lat_ <= 90.0
/// and longitude lies in the range: -180.0 <= lon_ <= 180.0.
template <typename T>
  requires std::floating_point<T>
class LatLong final {
#ifdef PYBIND11_NUMPY_DTYPE
public:
#endif
  Degrees<T> lat_ = {0};
  Degrees<T> lon_ = {0};

#ifndef PYBIND11_NUMPY_DTYPE
public:
#endif

  /// Contructor from latitude and longitude in Degrees.
  /// @pre the latitude and longitude must be valid.
  /// @param lat, lon the latitude and longitude in Degrees.
  constexpr LatLong(const Degrees<T> lat, const Degrees<T> lon)
      : lat_{lat}, lon_{lon} {
    Ensures(is_valid());
  }

  /// Contructor from a unit vector.
  /// @pre a must be a unit vector.
  /// @param a the unit vector.
  explicit constexpr LatLong(const vector::Vector3<T> &a)
      : LatLong(vector::latitude(a).to_degrees(),
                vector::longitude(a).to_degrees()) {
    Expects(vector::is_unit(a));
  }

  /// Function to determine whether the LatLong is valid.
  [[nodiscard("Pure Function")]]
  constexpr auto is_valid() const noexcept -> bool {
    return is_valid_latitude(lat_.v()) && is_valid_longitude(lon_.v());
  }

  /// Accessor for the latitude.
  [[nodiscard("Pure Function")]]
  constexpr auto lat() const noexcept -> Degrees<T> {
    return lat_;
  }

  /// Accessor for the longitude.
  [[nodiscard("Pure Function")]]
  constexpr auto lon() const noexcept -> Degrees<T> {
    return lon_;
  }

  /// Determine whether the `LatLong` is South of a.
  /// It compares the latitude of the two points.
  /// @param a the other point.
  ///
  /// @return true if a is South of b, false otherwise.
  [[nodiscard("Pure Function")]]
  constexpr auto is_south_of(const LatLong<T> &a) const noexcept -> bool {
    return lat_.v() < a.lat_.v();
  }

  /// Determine whether the `LatLong` is West of a.
  /// It compares the longitude difference of the two points.
  /// @param a the other point.
  ///
  /// @return true if a is West of b, false otherwise.
  [[nodiscard("Pure Function")]]
  constexpr auto is_west_of(const LatLong<T> &a) const noexcept -> bool {
    return (a.lon_ - lon_).v() < T();
  }

  /// Convert a `LatLong` to a point on the unit sphere
  [[nodiscard("Pure Function")]]
  constexpr auto to_point() const -> vector::Vector3<T> {
    return vector::to_point(Angle<T>(lat_), Angle<T>(lon_));
  }

  /// A Python representation of an LatLong type.
  /// I.e.: LatLong([lat, lon])
  /// @return a string in Python repr format.
  constexpr auto python_repr() const noexcept -> std::string {
    return "LatLong([ " + std::to_string(lat_.v()) + ", " +
           std::to_string(lon_.v()) + " ])";
  }
};

/// LatLong equality operator
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto operator==(const LatLong<T> &lhs, const LatLong<T> &rhs) noexcept
    -> bool {
  return lhs.lat() == rhs.lat() && lhs.lon() == rhs.lon();
}

/// LatLong ostream << operator
template <typename T>
  requires std::floating_point<T>
constexpr auto operator<<(std::ostream &os, const LatLong<T> &a)
    -> std::ostream & {
  return os << '(' << a.lat() << ',' << a.lon() << ')';
}

/// Calculate the azimuth and distance along the great circle of point b from
/// point a.
/// @param a, b the start and end positions
///
/// @return the great-circle azimuth relative to North and distance of point b
/// from point a.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto calculate_azimuth_and_distance(const LatLong<T> a,
                                              const LatLong<T> b) noexcept
    -> std::tuple<Angle<T>, Radians<T>> {
  const Angle a_lat{a.lat()};
  const Angle b_lat{b.lat()};
  const Angle delta_long(b.lon(), a.lon());
  return {great_circle::calculate_gc_azimuth(a_lat, b_lat, delta_long),
          great_circle::calculate_gc_distance(a_lat, b_lat, delta_long)};
}

/// Calculate the distance along the great circle of point b from point a.
///
/// See: [Haversine formula](https://en.wikipedia.org/wiki/Haversine_formula).
/// This function is less accurate than `calculate_azimuth_and_distance`.
/// @param a, b the start and end positions
///
/// @return the great-circle distance of point b from point a in `Radians`.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto haversine_distance(const LatLong<T> a,
                                  const LatLong<T> b) noexcept -> Radians<T> {
  const Angle a_lat{a.lat()};
  const Angle b_lat{b.lat()};
  const Angle delta_long(b.lon(), a.lon());
  const Angle delta_lat(b.lat(), a.lat());
  return great_circle::calculate_haversine_distance(a_lat, b_lat, delta_long,
                                                    delta_lat);
}

/// An `Arc` of a Great Circle on a unit sphere.
/// @invariant `a_` and `pole_` are orthogonal unit vectors,
/// length_ and half_width_ are not negative.
template <typename T>
  requires std::floating_point<T>
class Arc final {
#ifdef PYBIND11_NUMPY_DTYPE
public:
#endif
  /// The start point of the `Arc`.
  vector::Vector3<T> a_ = {};
  /// The right hand pole of the Great Circle of the `Arc`.
  vector::Vector3<T> pole_ = {};
  /// The length of the `Arc`.
  Radians<T> length_ = {Radians<T>(0)};
  /// The half width of the `Arc`.
  Radians<T> half_width_ = {Radians<T>(0)};

#ifndef PYBIND11_NUMPY_DTYPE
public:
#endif
  /// Arc constructor
  /// @param a the start point of the `Arc`.
  /// @param pole the right hand pole of the Great Circle of the `Arc`.
  /// @param length the length of the `Arc`.
  /// @param half_width the half width of the `Arc`.
  constexpr Arc(const vector::Vector3<T> a, const vector::Vector3<T> pole,
                const Radians<T> length,
                const Radians<T> half_width = Radians<T>(0))
      : a_{a}, pole_{pole}, length_{length}, half_width_{half_width} {}

  /// Arc constructor
  /// @param a start position
  /// @param azimuth the azimuth at a.
  /// @param length the length of the `Arc`.
  /// @param half_width the half width of the `Arc`.
  constexpr Arc(const LatLong<T> a, const Angle<T> azimuth,
                const Radians<T> length,
                const Radians<T> half_width = Radians<T>(0))
      : a_{a.to_point()}, pole_{vector::calculate_pole<T>(
                              Angle(a.lat()), Angle(a.lon()), azimuth)},
        length_{length}, half_width_{half_width} {}

  /// Arc constructor from start and end positions
  /// @param a, b start and end positions
  /// @param half_width the half width of the `Arc`.
  constexpr Arc(const LatLong<T> a, const LatLong<T> b,
                const Radians<T> half_width = Radians<T>(0)) {
    const auto [azimuth, length]{calculate_azimuth_and_distance<T>(a, b)};
    const auto a_lat{Angle(a.lat())};
    // if a is at the North or South pole
    if (a_lat.cos().v() < great_circle::MIN_VALUE<T>)
      *this = Arc<T>(LatLong(a.lat(), b.lon()), azimuth, length, half_width);
    else
      *this = Arc<T>(a, azimuth, length, half_width);
  }

  /// Set the `half_width` of an `Arc`
  /// @param  half_width the half width of the `Arc`.
  constexpr auto set_half_width(const Radians<T> half_width) noexcept -> void {
    half_width_ = half_width;
  }

  /// Test whether an `Arc` is valid.
  /// @return true if both a and pole are orthogonal unit vectors
  /// and both length and `half_width` are not negative.
  [[nodiscard("Pure Function")]]
  constexpr auto is_valid() const noexcept -> bool {
    return vector::is_unit(a_) && vector::is_unit(pole_) &&
           vector::are_orthogonal(a_, pole_) && !std::signbit(length_.v()) &&
           !std::signbit(half_width_.v());
  }

  /// The start point of the `Arc`.
  [[nodiscard("Pure Function")]]
  constexpr auto a() const noexcept -> vector::Vector3<T> {
    return a_;
  }

  /// The right hand pole of the Great Circle at the start point of the `Arc`.
  [[nodiscard("Pure Function")]]
  constexpr auto pole() const noexcept -> vector::Vector3<T> {
    return pole_;
  }

  /// The length of the `Arc`.
  [[nodiscard("Pure Function")]]
  constexpr auto length() const noexcept -> Radians<T> {
    return length_;
  }

  /// The half width of the `Arc`.
  [[nodiscard("Pure Function")]]
  constexpr auto half_width() const noexcept -> Radians<T> {
    return half_width_;
  }

  /// The azimuth at the start point.
  [[nodiscard("Pure Function")]]
  constexpr auto azimuth() const noexcept -> Angle<T> {
    return vector::calculate_azimuth(a_, pole_);
  }

  /// The direction vector of the `Arc` at the start point.
  [[nodiscard("Pure Function")]]
  constexpr auto direction() const noexcept -> vector::Vector3<T> {
    return vector::direction(a_, pole_);
  }

  /// A position vector at distance along the `Arc`.
  [[nodiscard("Pure Function")]]
  constexpr auto position(const Radians<T> distance) const noexcept
      -> vector::Vector3<T> {
    return vector::position(a_, direction(), Angle<T>(distance));
  }

  /// The end point of the `Arc`.
  [[nodiscard("Pure Function")]]
  constexpr auto b() const noexcept -> vector::Vector3<T> {
    return position(length_);
  }

  /// The mid point of the `Arc`.
  [[nodiscard("Pure Function")]]
  constexpr auto mid_point() const noexcept -> vector::Vector3<T> {
    return position(length_.half());
  }

  /// The position of a perpendicular point at distance from the `Arc`.
  /// @param point a point on the `Arc`'s great circle.
  /// @param distance the perpendicular distance from the `Arc`'s great
  /// circle.
  ///
  /// @return the point at perpendicular distance from point.
  [[nodiscard("Pure Function")]]
  constexpr auto perp_position(const vector::Vector3<T> &point,
                               const Radians<T> distance) const noexcept
      -> vector::Vector3<T> {
    return vector::position(point, pole_, Angle<T>(distance));
  }

  /// The position of a point at angle from the `Arc` start, at `Arc` length.
  /// @param angle the angle from the `Arc` start.
  ///
  /// @return the point at angle from the `Arc` start, at `Arc` length.
  [[nodiscard("Pure Function")]]
  constexpr auto angle_position(const Angle<T> angle) const noexcept
      -> vector::Vector3<T> {
    return vector::rotate_position(a_, pole_, angle, Angle<T>(length_));
  }

  /// The `Arc` at the end of an `Arc`, just the point if `half_width` is
  /// zero.
  /// @param `at_b` if true the `Arc` at b, else the `Arc` at a.
  ///
  /// @return the end `Arc` at a or b.
  [[nodiscard("Pure Function")]]
  constexpr auto end_arc(const bool at_b) const noexcept -> Arc<T> {
    const vector::Vector3<T> p{at_b ? b() : a_};
    const vector::Vector3<T> pole{vector::direction(p, pole_)};
    if (half_width_.v() < great_circle::MIN_VALUE<T>) {
      return Arc(p, pole, Radians<T>(0.0));
    } else {
      const vector::Vector3<T> a{perp_position(p, half_width_)};
      return Arc(a, pole, half_width_ + half_width_);
    }
  }

  /// Calculate great-circle along and across track distances of point from
  /// the `Arc`.
  /// @param point the point.
  ///
  /// @return the along and across track distances of the point in Radians.
  [[nodiscard("Pure Function")]]
  constexpr auto
  calculate_atd_and_xtd(const vector::Vector3<T> &point) const noexcept
      -> std::tuple<Radians<T>, Radians<T>> {
    return vector::calculate_atd_and_xtd(a_, pole_, point);
  }

  /// Calculate the shortest great-circle distance of a point from the `Arc`.
  /// @param point the point.
  ///
  /// @return the shortest distance of a point from the `Arc` in Radians.
  [[nodiscard("Pure Function")]]
  constexpr auto
  shortest_distance(const vector::Vector3<T> &point) const noexcept
      -> Radians<T> {
    const auto [atd, xtd]{calculate_atd_and_xtd(point)};
    if (-great_circle::MIN_VALUE<T> <= atd.v() &&
        atd.v() <= length_.v() + great_circle::MIN_VALUE<T>) {
      return xtd.abs();
    } else {
      // adjust atd to measure the distance from the centre of the Arc
      const auto atd_centre{atd - length_.half()};
      const auto p{std::signbit(atd_centre.v()) ? a_ : b()};
      return great_circle::e2gc_distance(vector::distance(p, point));
    }
  }

  /// A Python representation of an Arc.
  /// @return a string in Python repr format.
  [[nodiscard("Pure Function")]]
  constexpr auto python_repr() const noexcept -> std::string {
    constexpr auto ARC_START("Arc([[ ");
    constexpr auto DELIM(" ");
    constexpr auto MID_POINT("],[");
    constexpr auto END_POINT("]],");
    constexpr auto FINISH(")");

    return ARC_START + std::to_string(a_(0)) + DELIM + std::to_string(a_(1)) +
           DELIM + std::to_string(a_(2)) + MID_POINT +
           std::to_string(pole_(0)) + DELIM + std::to_string(pole_(1)) + DELIM +
           std::to_string(pole_(2)) + END_POINT + std::to_string(length_.v()) +
           DELIM + std::to_string(half_width_.v()) + FINISH;
  }
};

/// Arc equality operator
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto operator==(const Arc<T> &lhs, const Arc<T> &rhs) noexcept
    -> bool {
  return lhs.a() == rhs.a() && lhs.pole() == rhs.pole() &&
         lhs.length() == rhs.length() && lhs.half_width() == rhs.half_width();
}

/// Calculate the great-circle distances along a pair of `Arc`s to their
/// closest intersection point or their coincident arc distances if the
/// `Arc`s are on coincident Great Circles.
/// @param arc_0, arc_1 the `Arc`s.
///
/// @return the distances along the first `Arc` and second `Arc` to the
/// intersection point or to their coincident arc distances if the `Arc`s do
/// not intersect.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto calculate_intersection_distances(const Arc<T> &arc_0,
                                                const Arc<T> &arc_1) noexcept
    -> std::tuple<Radians<T>, Radians<T>> {
  const auto [distance_0, distance_1, angle]{
      vector::intersection::calculate_arc_reference_distances_and_angle(
          arc_0.mid_point(), arc_0.pole(), arc_1.mid_point(), arc_1.pole())};

  return {distance_0 + arc_1.length().half(),
          distance_1 + arc_1.length().half()};
}

/// Calculate whether a pair of `Arc`s intersect and (if so) where.
/// @param arc_0, arc_1 the `Arc`s.
///
/// @return the intersection point or `std::nullopt` if they don't intersect.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto calculate_intersection_point(const Arc<T> &arc_0,
                                            const Arc<T> &arc_1) noexcept
    -> std::optional<vector::Vector3<T>> {
  const auto [point,
              angle]{vector::intersection::calculate_reference_point_and_angle(
      arc_0.mid_point(), arc_0.pole(), arc_1.mid_point(), arc_1.pole())};

  const Radians<T> distance_0{vector::calculate_great_circle_atd(
      arc_0.mid_point(), arc_0.pole(), point)};
  const Radians<T> distance_1{vector::calculate_great_circle_atd(
      arc_1.mid_point(), arc_1.pole(), point)};

  const bool arcs_are_coincident{angle.sin().v() == T()};
  const bool arcs_intersect_or_overlap{
      arcs_are_coincident
          ? distance_0.abs().v() + distance_1.abs().v() <=
                arc_0.length().half().v() + arc_1.length().half().v() +
                    great_circle::MIN_VALUE<T>
          : (distance_0.abs().v() <=
             arc_0.length().half().v() + great_circle::MIN_VALUE<T>) &&
                (distance_1.abs().v() <=
                 arc_1.length().half().v() + great_circle::MIN_VALUE<T>)};
  if (arcs_intersect_or_overlap) {
    return point;
  }

  return std::nullopt;
}

} // namespace via
