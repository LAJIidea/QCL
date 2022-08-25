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

    constexpr bool ignorable(const double v) {
        if (v > -epsilon && v < epsilon) return true;
        else return false;
    }

    template<typename Ty>
    constexpr bool ignorable(const std::complex<Ty>& v) {
        Ty value = abs_sqr(v);
        if (ignorable(value)) return true;
        else return false;
    }

    constexpr std::pair<size_t, size_t> get_layer_range(size_t layer_id)
    {
        size_t lower = pow2(layer_id + 1) -2;
        size_t upper = pow2(layer_id + 2) -3;
        return {lower, upper};
    }

    constexpr std::pair<size_t, size_t> get_layer_node_range(size_t layer_id)
    {
        size_t lower = pow2(layer_id) - 1;
        size_t upper = pow2(layer_id + 1) -2;
        return {lower, upper};
    }

    template<typename T>
    struct _remove_cv_ref {
        using type = std::remove_cv_t<std::remove_reference_t<T>>;
    };

    template<typename T>
    using _remove_cverf_t = typename _remove_cv_ref<T>::type;

    template<typename Ty>
    void* to_void_ptr(Ty ptr) {
        using T_ptr_t = _remove_cv_ref<Ty>;
        using T = _remove_cverf_t<std::remove_pointer_t<T_ptr_t>>;
        using clear_pointer_type = T *;
        return reinterpret_cast<void*>(const_cast<clear_pointer_type>(ptr));
    }

    template<typename KeyTy, typename ValTy>
    void map2vec(std::vector<std::pair<void*, void*>>& vec, const std::map<KeyTy, ValTy>& map1) {
        vec.clear();
        vec.reserve(map1.size());
        for (const auto& item : map1) {
            void* key_ptr = to_void_ptr(&(item.first));
            void* val_ptr = to_void_ptr(&(item.second));
            // vec.push_back({key_ptr, val_ptr});
            vec.emplace_back(key_ptr, val_ptr);
        }
    }

    template<typename EngineType, typename MemoryContainer>
    void random_memory(MemoryContainer& memory, size_t memory_size, EngineType &engine) {
        size_t size = memory.size();
        std::uniform_int_distribution<memory_entry_t> ud(0, pow2(memory_size) - 1);

        for (auto iter = std::begin(memory); iter != std::end(memory); ++iter) {
            *iter = ud(engine);
        }
    }

    inline void random_memory(std::vector<size_t>& memory, size_t memory_size) {
        random_memory(memory, memory_size,random_engine::get_engine());
    }

    template<typename FwdIt, typename Pred, typename Func>
    FwdIt unique_and_merge(FwdIt first, FwdIt last, Pred pred, Func fn)
    {
        if (first == last) return last;

        FwdIt result = first;
        while (++first != last)
        {
            if (!pred(*result, *first))
                *(++result) = *first;
            else
                fn(*result, *first);
        }
        return ++result;
    }

    template< class Key, class T, class Compare, class Alloc, class Pred >
    void erase_if(std::map<Key, T, Compare, Alloc>& c, Pred pred) {
        for (auto i = c.begin(), last = c.end(); i != last; ) {
            if (pred(*i)) {
                i = c.erase(i);
            }
            else {
                ++i;
            }
        }
    }

    template<typename Rng>
    void choice_from(std::set<size_t>& samples, int size, size_t n_samples, Rng& g)
    {
        samples.clear();
        std::uniform_int_distribution<size_t> ud(0, 1ull << size);
        while (n_samples > 0) {
            if (samples.insert(ud(g)).second) { n_samples--; };
        }
    }

    inline std::vector<double> lin_space(double min, double max, size_t points) {
        double delta = (max - min) / (points - 1);
        std::vector<double> ret;
        ret.reserve(points);
        for (size_t i = 0; i < points; ++i) {
            ret.push_back(min + delta * i);
        }
        return ret;
    }

    inline std::pair<double, double> mean_std(const std::vector<double> &m) {
        auto sq = [](double m, double y) {
            return m + y * y;
        };

        double sum = std::accumulate(m.begin(), m.end(), 0.0);
        double sumsq = std::accumulate(m.begin(), m.end(), 0.0, sq);
        double mean = sum / m.size();
        double meansq = sumsq / m.size();
        return { mean, sqrt(meansq - mean * mean) };

    }

    inline const char* bool2char(bool x)
    {
        return x ? "1" : "0";
    }

    inline constexpr int bool2int(bool x)
    {
        return x ? 1 : 0;
    }

    inline const char* bool2char_pm_basis(bool x)
    {
        return x ? "+" : "-";
    }


    inline std::string complex2str(const std::complex<double>& x)
    {
        if (x.imag() > 0)
            return fmt::format("{}+{}j", x.real(), x.imag());
        else
            return fmt::format("{}-{}j", x.real(), -x.imag());
    }

    inline int bit_count(size_t n)
    {
        int count = 0;
        while (n) {
            count++;
            n &= (n - 1);
        }
        return count;
    }

    /* Check the nan */
    inline void check_nan(const complex_t &m)
    {
#ifndef QRAM_Release
        if (std::isnan(m.real()) || std::isnan(m.imag()))
        {
            throw std::runtime_error("Nan!");
        }
#endif
    }

    /* Check the nan */
    inline void check_nan(const double &m)
    {
#ifndef QRAM_Release
        if (std::isnan(m))
        {
            throw std::runtime_error("Nan!");
        }
#endif
    }

    inline void throw_not_implemented()
    {
#ifndef QRAM_Release
        throw std::runtime_error("Not implemented.");
#endif
    }

    inline void throw_bad_switch_case()
    {
#ifndef QRAM_Release
        throw std::runtime_error("Impossible switch branch.");
#endif
    }

    inline void throw_invalid_input()
    {
#ifndef QRAM_Release
        throw std::runtime_error("Invalid input.");
#endif
    }

    inline void throw_bad_result()
    {
#ifndef QRAM_Release
        throw std::runtime_error("Bad result.");
#endif

    }

    template<typename Ty = void>
    struct array_length
    {};

    template<typename DTy, int sz>
    struct array_length<std::array<DTy, sz>>
    {
        static constexpr int value = sz;
    };

    template<typename DTy, int sz>
    struct array_length<DTy[sz]>
    {
        static constexpr int value = sz;
    };

} // namespace qram_simulator


template <> struct fmt::formatter<std::complex<double>> {
    // Presentation format: 'f' - fixed, 'e' - exponential.
    char presentation = 'f';

    // Parses format specifications of the form ['f' | 'e'].
    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
        // [ctx.begin(), ctx.end()) is a character range that contains a part of
        // the format string starting from the format specifications to be parsed,
        // e.g. in
        //
        //   fmt::format("{:f} - point of interest", point{1, 2});
        //
        // the range will contain "f} - point of interest". The formatter should
        // parse specifiers until '}' or the end of the range. In this example
        // the formatter should parse the 'f' specifier and return an iterator
        // pointing to '}'.

        // Parse the presentation format and store it in the formatter:
        auto it = ctx.begin(), end = ctx.end();
        if (it != end && (*it == 'f' || *it == 'e')) presentation = *it++;

        // Check if reached the end of the range:
        if (it != end && *it != '}') throw format_error("invalid format");

        // Return an iterator past the end of the parsed range:
        return it;
    }

    // Formats the point p using the parsed format specification (presentation)
    // stored in this formatter.
    template <typename FormatContext>
    auto format(const std::complex<double>& p, FormatContext& ctx) -> decltype(ctx.out()) {
        // ctx.out() is an output iterator to write to.
        if (p.imag() >= 0)
            return presentation == 'f'
                   ? format_to(ctx.out(), "{:f}+{:f}i", p.real(), p.imag())
                   : format_to(ctx.out(), "{:e}+{:e}i", p.real(), p.imag());
        else
            return presentation == 'f'
                   ? format_to(ctx.out(), "{:f}{:f}i", p.real(), p.imag())
                   : format_to(ctx.out(), "{:e}{:e}i", p.real(), p.imag());
    }
};

namespace std {

    template<size_t sz>
    inline bool operator<(const std::bitset<sz>& lhs, const std::bitset<sz>& rhs)
    {
        return lhs.to_ullong() < rhs.to_ullong();
    }

    template<typename Ty>
    std::vector<Ty> operator+(const std::vector<Ty>& lhs, const std::vector<Ty>& rhs)
    {
        std::vector<Ty> ret;
        if (lhs.size() == rhs.size())
        {
            ret.resize(lhs.size());
            for (size_t i = 0; i < lhs.size(); ++i)
            {
                ret[i] = lhs[i] + rhs[i];
            }
        }
        return lhs;
    }

    template<typename Ty>
    std::vector<Ty>& operator+=(std::vector<Ty>& lhs, const std::vector<Ty>& rhs)
    {
        if (lhs.size() == rhs.size())
            for (size_t i = 0; i < lhs.size(); ++i)
            {
                lhs[i] += rhs[i];
            }
        return lhs;
    }
}



#endif //QCL_UTIL_H
