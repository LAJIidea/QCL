//
// Created by BY210033 on 2022/8/25.
//

#ifndef QCL_UTIL_H
#define QCL_UTIL_H

#include <algorithm>
#include <array>
#include <bitset>
#include <complex>
#include <cmath>
#include <deque>
#include <functional>
#include <iostream>
#include <list>
#include <map>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include <unordered_map>
#include <vector>

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include "random_engine.h"

#ifdef QCL_Release
#define QRAM_NOEXCEPT noexcept
#else
#define QRAM_NOEXCEPT
#endif

namespace qram_simulator {

    using complex_t = std::complex<double>;
    using memory_entry_t = size_t;
    using memory_t = std::vector<memory_entry_t>;
    using bus_t = memory_entry_t;

    constexpr double pi = 3.141592653589793238462643383279502884L;
    constexpr double sqrt2 = 1.41421356237309504880168872420969807856967;
    constexpr double sqrt2inv = 1.0 / sqrt2;
    constexpr double epsilon = 1.e-7;

    template<typename Ty>
    auto abs_sqr(const std::complex<Ty>& c) -> Ty {
        return c.real() * c.real() + c.imag() * c.imag();
    }

    inline constexpr size_t pow2(size_t n) { return (1ull) << (n); }

    inline constexpr bool get_digit(size_t n, size_t digit) {
        return (n >> digit) & 1;
    }

    inline constexpr bool get_digit_reverse(size_t n, size_t digit, size_t max_digit) {
        return (n >> (max_digit - digit - 1)) % 2;
    }

    inline constexpr size_t log2(size_t n) {
        size_t ret = 0;
        while (n > 1) {
            ret++;
            n /= 2;
        }
        return ret;
    }

    inline bool operator<(const complex_t& lhs, const complex_t& rhs)
    {
        return (abs_sqr(lhs) < abs_sqr(rhs));
    }

    template<typename Ty>
    Ty amp_sum(const std::vector<std::complex<Ty>>& amps)
    {
        Ty val = 0;
        for (auto a : amps) {
            val += abs_sqr(a);
        }
        return val;
    }

    constexpr inline bool digit1(size_t i, size_t digit)
    {
        return (i >> digit) & 1;
    }

    constexpr inline bool digit0(size_t i, size_t digit)
    {
        return !digit1(i, digit);
    }

    constexpr size_t flip_digit(size_t i, size_t digit)
    {
        auto m = (1 << digit);
        return digit1(i, digit) ? (i -= m) : (i += m);
    }

}

#endif //QCL_UTIL_H
