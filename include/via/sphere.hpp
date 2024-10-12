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
#include "via/sphere/intersection.hpp"

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

/// The Angle represents an angle by it's sine and cosine components.
/// @invariant sin() * sin() + cos() * cos() = 1
template <typename T>
  requires std::floating_point<T>
class LatLong {
#ifdef PYBIND11_NUMPY_DTYPE
public:
#endif
  Degrees<T> lat_ = {0};
  Degrees<T> lon_ = {0};

#ifndef PYBIND11_NUMPY_DTYPE
public:
#endif

  constexpr LatLong(const Degrees<T> lat, const Degrees<T> lon)
      : lat_{lat}, lon_{lon} {
#ifndef PYBIND11_VERSION_MAJOR
    Ensures(is_valid());
#endif
  }

  /// Function to determine whether the LatLong is valid.
  [[nodiscard("Pure Function")]]
  constexpr auto is_valid() const noexcept -> bool {
    return is_valid_latitude(lat_.v()) && is_valid_longitude(lon_.v());
  }

  [[nodiscard("Pure Function")]]
  constexpr auto lat() const noexcept -> Degrees<T> {
    return lat_;
  }

  [[nodiscard("Pure Function")]]
  constexpr auto lon() const noexcept -> Degrees<T> {
    return lon_;
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

/// Convert a `LatLong` to a point on the unit sphere
/// @param a - the LatLong.
///
/// returns a `Vector3` of the point on the unit sphere.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto to_point(const LatLong<T> &a) -> vector::Vector3<T> {
  return vector::to_point(Angle<T>(a.lat()), Angle<T>(a.lon()));
}

/// Convert a point on the unit sphere to a `LatLong`
/// @param a - the `Vector3` of the point .
///
/// returns the LatLong of the point.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto to_lat_long(const vector::Vector3<T> &a) -> LatLong<T> {
  const auto lat{vector::latitude(a)};
  const auto lon{vector::longitude(a)};
  return LatLong<T>(lat.to_degrees(), lon.to_degrees());
}

/// Calculate the azimuth and distance along the great circle of point b from
/// point a.
/// @aparam a, b the start and end positions
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

/// Calculate the distance along the great circle of point b from point a,
/// see: [Haversine formula](https://en.wikipedia.org/wiki/Haversine_formula).
/// This function is less accurate than `calculate_azimuth_and_distance`.
/// @aparam a, b the start and end positions
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

template <typename T>
  requires std::floating_point<T>
class Arc {
#ifdef PYBIND11_NUMPY_DTYPE
public:
#endif
  vector::Vector3<T> a_ = {};
  vector::Vector3<T> pole_ = {};
  Radians<T> length_ = {Radians<T>(0)};
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
      : a_{to_point<T>(a)}, pole_{vector::calculate_pole<T>(
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

  [[nodiscard("Pure Function")]]
  constexpr auto is_valid() const noexcept -> bool {
    return vector::is_unit(a_) && vector::is_unit(pole_) &&
           vector::are_orthogonal(a_, pole_) && T(0) <= length_.v() &&
           T(0) <= half_width_.v();
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
    return position(Radians<T>(length_.v() / 2));
  }

  /// The position of a perpendicular point at distance from the `Arc`.
  /// @param point a point on the `Arc`'s great circle.
  /// @param distance the perpendicular distance from the `Arc`'s great circle.
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

  /// The `Arc` at the end of an `Arc`, just the point if `half_width` is zero.
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
};

/// Calculate the great-circle distances along a pair of `Arc`s to their
/// closest intersection point or their coincident arc distances if the
/// `Arc`s are on coincident Great Circles.
/// @param arc1, arc2 the `Arc`s.
///
/// @return the distances along the first `Arc` and second `Arc` to the
/// intersection point or to their coincident arc distances if the `Arc`s do not
/// intersect.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto calculate_intersection_distances(const Arc<T> &arc1,
                                                const Arc<T> &arc2) noexcept
    -> std::tuple<Radians<T>, Radians<T>> {
  const auto centroid{(arc1.mid_point() + arc2.mid_point()) / 2};
  return vector::intersection::calculate_intersection_point_distances<T>(
      arc1.a(), arc1.pole(), arc1.length(), arc2.a(), arc2.pole(),
      arc2.length(), centroid);
}

/// Calculate whether a pair of `Arc`s intersect and (if so) where.
/// @param arc1, arc2 the `Arc`s.
///
/// @return the intersection point or `std::nullopt` if they don't intersect.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto calculate_intersection_point(const Arc<T> &arc1,
                                            const Arc<T> &arc2) noexcept
    -> std::optional<vector::Vector3<T>> {
  const auto [distance1,
              distance2]{calculate_intersection_distances(arc1, arc2)};
  if (vector::intersection::is_within(distance1.v(), arc1.length().v()) &&
      vector::intersection::is_within(distance2.v(), arc2.length().v())) {
    return arc1.position(distance1);
  }

  return std::nullopt;
}

} // namespace via