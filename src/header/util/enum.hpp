/**
* @author Tim Luchterhand
* @date 04.10.24
* @brief enum utilities
*/

#ifndef TEMPO_ENUM_HPP
#define TEMPO_ENUM_HPP

#include <ostream>
#include <iterator>
#include <array>
#include <ranges>
#include <type_traits>
#include <string_view>

/**
 * Converts enum to underlying type
 * @tparam E enum type
 * @param e enum to convert
 * @return value of underlying type
 * @note replace with std impl in c++23: https://en.cppreference.com/w/cpp/utility/to_underlying
 */
template<typename E>
constexpr auto to_underlying(E e) noexcept {
    return static_cast<std::underlying_type_t<E>>(e);
}

namespace enum_detail {
    template<const std::string_view &String, char Delim>
    class split {
        static consteval auto numParts() {
            unsigned ret = 1;
            for (char c: String) {
                ret += (c == Delim);
            }

            return ret;
        }

        static consteval auto extractString(auto begin, auto end) {
            auto str = std::ranges::subrange(begin, end) |
                       std::views::drop_while([](char c) { return c == ' '; });
            return std::string_view(str.begin(), str.end());
        }

        static consteval auto doSplit() {
            std::array<std::string_view, numParts()> ret;
            auto start = String.begin();
            std::size_t arrayIdx = 0;
            for (auto curr = String.begin(); curr != String.end(); ++curr) {
                if (*curr == Delim) {
                    ret[arrayIdx++] = extractString(start, curr);
                    start = curr + 1;
                }
            }

            ret[arrayIdx] = extractString(start, String.end());
            return ret;
        }

    public:
        static constexpr auto value = doSplit();
    };

    template<typename T>
    struct enum_conversion {
        static constexpr bool enabled = false;
    };

    template<typename T>
    concept penum = std::is_enum_v<T> and requires(T val, T &ref, std::string_view sv) {
        str_to_penum(sv, ref);
        { penum_to_string(val) } -> std::same_as<std::string>;
    };
}

#define PENUM(NAME, ...) enum class NAME { __VA_ARGS__ };                                                       \
inline constexpr std::string_view __##NAME##_enum_str_vals__ = #__VA_ARGS__;                                    \
inline constexpr std::array __##NAME##_converter__{enum_detail::split<__##NAME##_enum_str_vals__, ','>::value };\
inline std::ostream &operator<<(std::ostream &os, NAME e) {                                                     \
    os << __##NAME##_converter__.at(to_underlying(e));                                                          \
    return os;                                                                                                  \
}                                                                                                               \
constexpr auto penum_to_string(NAME e) {                                                                        \
    return std::string(__##NAME##_converter__.at(to_underlying(e)));                                            \
}                                                                                                               \
constexpr void str_to_penum(std::string_view str, NAME &out) {                                                  \
    using namespace std::ranges;                                                                                \
    using namespace std::literals;                                                                              \
    auto res = find(__##NAME##_converter__, str);                                                               \
    if (res != __##NAME##_converter__.end()) {                                                                  \
        out =  static_cast<NAME>(distance(__##NAME##_converter__.begin(), res));                                \
        return;                                                                                                 \
    }                                                                                                           \
                                                                                                                \
    throw std::runtime_error("cannot convert '"s + std::string(str) + "' to enum "s + #NAME);                   \
}                                                                                                               \
constexpr auto str_to_##NAME(std::string_view str) {                                                            \
NAME ret;                                                                                                       \
    str_to_penum(str, ret);                                                                                     \
    return ret;                                                                                                 \
}

#endif //TEMPO_ENUM_HPP
