/**
* @author Tim Luchterhand
* @date 03.09.24
* @brief
*/

#ifndef TEMPO_INTFINITY_HPP
#define TEMPO_INTFINITY_HPP

#include <concepts>
#include <limits>
#include <compare>
#include <ostream>
#include <cmath>

namespace tempo {
    namespace detail {
        template<typename T, bool B, T TVal, T FVal>
        struct conditional {
        };

        template<typename T, T TVal, T FVal>
        struct conditional<T, true, TVal, FVal> {
            static constexpr T value = TVal;
        };

        template<typename T, T TVal, T FVal>
        struct conditional<T, false, TVal, FVal> {
            static constexpr T value = FVal;
        };

        template<typename T, bool B, T TVal, T FVal>
        inline constexpr auto conditional_v = conditional<T, B, TVal, FVal>::value;

        template<typename T>
        constexpr T sgn(T x) noexcept {
            if (x < 0) {
                return -1;
            }

            if (x > 0) {
                return 1;
            }

            return 0;
        }
    }

    template<std::integral T, bool UnsignedUnderflow = false>
    class intfinity {
        T value;
        static constexpr bool TwoCompl = not std::is_signed_v<T> or T(-1) == compl T(0);
        static_assert(TwoCompl, "This implementation only works on machines that use the 2's complement for signed"
                                "integers. I suggest using a machine from past WW2 :D");
        static constexpr T NotANumber = detail::conditional_v<T, std::is_signed_v<T>,
                std::numeric_limits<T>::min(), std::numeric_limits<T>::max()>;
        static constexpr T Infinity =
                std::numeric_limits<T>::max() - detail::conditional_v<T, std::is_signed_v<T>, 0, 1>;
    public:
        static constexpr intfinity Inf() { return Infinity; }

        static constexpr intfinity Nan() { return NotANumber; }

        constexpr intfinity() noexcept = default;

        constexpr intfinity(T value) noexcept: value(value) {}

        template<std::floating_point F>
        constexpr intfinity(F value) noexcept: value(value) {
            constexpr auto FInf = static_cast<F>(Infinity);
            if (std::isnan(value)) {
                this->value = NotANumber;
                return;
            }

            if constexpr (std::is_signed_v<T>) {
                if (std::abs(value) > FInf) {
                    this->value = static_cast<T>(detail::sgn(value)) * FInf;
                }
            } else {
                if (value < 0) {
                    if constexpr (UnsignedUnderflow) {
                        this->value = Infinity;
                    } else {
                        this->value = 0;
                    }
                } else if (value > FInf) {
                    this->value = Infinity;
                }
            }
        }

        [[nodiscard]] constexpr bool isInf() const noexcept {
            if constexpr (std::is_signed_v<T>) {
                return std::abs(value) == Infinity;
            } else {
                return value == Infinity;
            }
        }

        [[nodiscard]] constexpr bool isNan() const noexcept {
            return value == NotANumber;
        }

        constexpr bool operator==(intfinity other) const noexcept {
            if (isNan() or other.isNan()) { return false; }
            return value == other.value;
        }

        template<std::floating_point F>
        constexpr bool operator==(F other) const noexcept {
            return static_cast<F>(*this) == other;
        }

        constexpr std::partial_ordering operator<=>(intfinity other) const noexcept {
            if (isNan() or other.isNan()) { return std::partial_ordering::unordered; }
            return value <=> other.value;
        }

        template<std::floating_point F>
        constexpr auto operator<=>(F other) const noexcept {
            return static_cast<F>(*this) <=> other;
        }

        constexpr intfinity &operator+=(intfinity other) noexcept {
            if (isNan() or other.isNan()) {
                value = NotANumber;
                return *this;
            }

            if constexpr (std::is_signed_v<T>) {
                if (isInf() and other.isInf() and detail::sgn(value) != detail::sgn(other.value)) {
                    value = NotANumber;
                    return *this;
                }
            }

            if (isInf()) {
                return *this;
            }

            if constexpr (std::is_signed_v<T>) {
                if (other.value > 0 and Infinity - other.value < value) {
                    value = Infinity;
                    return *this;
                } else if (other.value < 0 and -Infinity - other.value > value) {
                    value = -Infinity;
                    return *this;
                }
            } else {
                if (Infinity - other.value < value) {
                    value = Infinity;
                    return *this;
                }
            }

            value += other.value;
            return *this;
        }

        template<std::floating_point F>
        constexpr intfinity &operator+=(F other) noexcept {
            *this = static_cast<F>(*this) + other;
            return *this;
        }

        constexpr intfinity &operator-=(intfinity other) noexcept {
            if (isNan() or other.isNan()) {
                value = NotANumber;
                return *this;
            }

            if constexpr (std::is_signed_v<T>) {
                if (isInf() and other.isInf() and detail::sgn(value) == detail::sgn(other.value)) {
                    value = NotANumber;
                    return *this;
                }
            } else {
                if (isInf() and other.isInf()) {
                    value = NotANumber;
                    return *this;
                }
            }

            if (isInf()) {
                return *this;
            }

            if constexpr (std::is_signed_v<T>) {
                if (other.value < 0 and Infinity + other.value < value) {
                    value = Infinity;
                    return *this;
                } else if (other.value > 0 and -Infinity + other.value > value) {
                    value = -Infinity;
                    return *this;
                }
            } else {
                if (value < other.value) {
                    if constexpr (UnsignedUnderflow) {
                        value = Infinity;
                    } else {
                        value = 0;
                    }

                    return *this;
                }
            }

            value -= other.value;
            return *this;
        }

        template<std::floating_point F>
        constexpr intfinity &operator-=(F other) {
            *this = static_cast<F>(*this) - other;
            return *this;
        }

        constexpr friend intfinity operator+(intfinity lhs, intfinity rhs) noexcept {
            lhs += rhs;
            return lhs;
        }

        template<std::floating_point F>
        friend constexpr F operator+(intfinity lhs, F rhs) noexcept {
            return static_cast<F>(lhs) + rhs;
        }

        template<std::floating_point F>
        friend constexpr F operator+(F lhs, intfinity rhs) noexcept {
            return lhs + static_cast<F>(rhs);
        }

        constexpr friend intfinity operator-(intfinity lhs, intfinity rhs) noexcept {
            lhs -= rhs;
            return lhs;
        }

        template<std::floating_point F>
        friend constexpr F operator-(intfinity lhs, F rhs) noexcept {
            return static_cast<F>(lhs) - rhs;
        }

        template<std::floating_point F>
        friend constexpr F operator-(F lhs, intfinity rhs) noexcept {
            return lhs - static_cast<F>(rhs);
        }

        constexpr intfinity &operator++() noexcept {
            return operator+=(1);
        }

        constexpr intfinity operator++(int) noexcept {
            auto old = *this;
            operator++();
            return old;
        }

        constexpr intfinity &operator--() noexcept {
            return operator-=(1);
        }

        constexpr intfinity operator--(int) noexcept {
            auto old = *this;
            operator--();
            return old;
        }

        constexpr intfinity &operator*=(intfinity other) noexcept {
            if (isNan() or other.isNan()) {
                value = NotANumber;
                return *this;
            }

            if ((isInf() and other.value == 0) or (value == 0 and other.isInf())) {
                value = NotANumber;
                return *this;
            }

            T lhs;
            T rhs;
            if constexpr (std::is_signed_v<T>) {
                lhs = std::abs(value);
                rhs = std::abs(other.value);
            } else {
                lhs = value;
                rhs = other.value;
            }

            if (rhs != 0 and Infinity / rhs < lhs) {
                value = Infinity;
                if constexpr (std::is_signed_v<T>) {
                    value *= detail::sgn(value) * detail::sgn(other.value);
                }

                return *this;
            }

            value *= other.value;
            return *this;
        }

        template<std::floating_point F>
        constexpr intfinity &operator*=(F other) noexcept {
            *this = static_cast<F>(*this) * other;
            return *this;
        }

        constexpr intfinity &operator/=(intfinity other) noexcept {
            if (isNan() or other.isNan()) {
                value = NotANumber;
                return *this;
            }

            if ((isInf() and other.isInf()) or (value == 0 and other.value == 0)) {
                value = NotANumber;
                return *this;
            }

            if (isInf()) {
                return *this;
            }

            if (other.value == 0) {
                if constexpr (std::is_signed_v<T>) {
                    value = detail::sgn(value) * Infinity;
                    return *this;
                } else {
                    value = Infinity;
                    return *this;
                }
            }

            value /= other.value;
            return *this;
        }

        template<std::floating_point F>
        constexpr intfinity &operator/=(F other) noexcept {
            *this = static_cast<F>(*this) / other;
            return *this;
        }

        constexpr friend intfinity operator*(intfinity lhs, intfinity rhs) noexcept {
            lhs *= rhs;
            return lhs;
        }

        template<std::floating_point F>
        friend constexpr F operator*(intfinity lhs, F rhs) noexcept {
            return static_cast<F>(lhs) * rhs;
        }

        template<std::floating_point F>
        friend constexpr F operator*(F lhs, intfinity rhs) noexcept {
            return lhs * static_cast<F>(rhs);
        }


        constexpr friend intfinity operator/(intfinity lhs, intfinity rhs) noexcept {
            lhs /= rhs;
            return lhs;
        }

        template<std::floating_point F>
        friend constexpr F operator/(intfinity lhs, F rhs) noexcept {
            return static_cast<F>(lhs) / rhs;
        }

        template<std::floating_point F>
        friend constexpr F operator/(F lhs, intfinity rhs) noexcept {
            return lhs / static_cast<F>(rhs);
        }

        constexpr intfinity operator+() const noexcept {
            return value;
        }

        template<std::signed_integral = T>
        constexpr intfinity operator-() const noexcept {
            if (not isNan()) {
                return -value;
            }

            return value;
        }

        template<std::floating_point F>
        explicit constexpr operator F() const noexcept {
            if (isNan()) {
                return std::numeric_limits<F>::quiet_NaN();
            }

            if (isInf()) {
                if constexpr (std::is_signed_v<T>) {
                    return static_cast<F>(detail::sgn(value)) * std::numeric_limits<F>::infinity();
                } else {
                    return std::numeric_limits<F>::infinity();
                }
            }

            return static_cast<F>(value);
        }

        friend std::ostream &operator<<(std::ostream &os, intfinity val) {
            if (val.isNan()) {
                os << "NaN";
            } else if (val.isInf()) {
                os << (val.get() < 0 ? "-" : "") << "Infinity";
            } else {
                os << val.get();
            }

            return os;
        }

        constexpr T get() const noexcept {
            return value;
        }
    };
}
namespace std {
    template<typename T, bool B>
    class numeric_limits<tempo::intfinity<T, B>> {
    public:
        static constexpr bool is_specialist = true;
        static constexpr bool is_signed = std::numeric_limits<T>::is_signed;
        static constexpr bool is_integer = true;
        static constexpr bool is_exact = true;
        static constexpr bool has_infinity = true;
        static constexpr bool has_quiet_NaN = true;
        static constexpr bool has_signaling_NaN = false;
        static constexpr auto has_denorm = std::denorm_absent;
        static constexpr bool has_denorm_loss = false;
        static constexpr auto round_style = std::round_toward_zero;
        static constexpr bool is_iec559 = false;
        static constexpr bool is_bounded = true;
        static constexpr bool is_modulo = false;
        static constexpr auto digits = std::numeric_limits<T>::digits;
        static constexpr auto digits10 = std::numeric_limits<T>::digits10;
        static constexpr auto max_digits10 = std::numeric_limits<T>::max_digits10;
        static constexpr auto radix = std::numeric_limits<T>::radix;
        static constexpr auto min_exponent = std::numeric_limits<T>::min_exponent;
        static constexpr auto min_exponent10 = std::numeric_limits<T>::min_exponent10;
        static constexpr auto max_exponent = std::numeric_limits<T>::max_exponent;
        static constexpr auto max_exponent10 = std::numeric_limits<T>::max_exponent10;
        static constexpr auto traps = false;
        static constexpr auto tinyness_before = false;

        static constexpr tempo::intfinity<T, B> min() noexcept {
            if constexpr (std::is_signed_v<T>) {
                return std::numeric_limits<T>::min() + 2;
            } else {
                return 0;
            }
        }

        static constexpr tempo::intfinity<T, B> max() noexcept {
            return std::numeric_limits<T>::max() - tempo::detail::conditional_v<T, std::is_signed_v<T>, 1, 2>;
        }

        static constexpr tempo::intfinity<T, B> lowest() noexcept {
            return min();
        }

        static constexpr tempo::intfinity<T, B> epsilon() noexcept {
            return 0;
        }

        static constexpr tempo::intfinity<T, B> round_error() noexcept {
            return 0;
        }

        static constexpr tempo::intfinity<T, B> denorm_min() noexcept {
            return 0;
        }

        static constexpr auto infinity() noexcept {
            return tempo::intfinity<T, B>::Inf();
        }

        static constexpr auto quiet_NaN() noexcept {
            return tempo::intfinity<T, B>::Nan();
        }
    };

    template<typename T, bool B>
    constexpr bool isinf(tempo::intfinity<T, B> val) noexcept {
        return val.isInf();
    }

    template<typename T, bool B>
    constexpr bool isfinite(tempo::intfinity<T, B> val) noexcept {
        return not val.isInf() and not val.isNan();
    }

    template<typename T, bool B>
    constexpr bool isnan(tempo::intfinity<T, B> val) noexcept {
        return val.isNan();
    }
}

#endif //TEMPO_INTFINITY_HPP
