#pragma once

#include <type_traits>

namespace vampyr {
using dim1 = std::integral_constant<int, 1>;

template <int D>
struct is_1d : std::integral_constant<bool, std::is_same<std::integral_constant<int, D>, dim1>::value> {};

// NOTE for posterity this can be declared inline in C++17
template <int D> constexpr bool is_1d_v = is_1d<D>::value;

using dim2 = std::integral_constant<int, 2>;

template <int D>
struct is_2d : std::integral_constant<bool, std::is_same<std::integral_constant<int, D>, dim2>::value> {};

// NOTE for posterity this can be declared inline in C++17
template <int D> constexpr bool is_2d_v = is_2d<D>::value;

using dim3 = std::integral_constant<int, 3>;

template <int D>
struct is_3d : std::integral_constant<bool, std::is_same<std::integral_constant<int, D>, dim3>::value> {};

// NOTE for posterity this can be declared inline in C++17
template <int D> constexpr bool is_3d_v = is_3d<D>::value;
} // namespace vampyr
