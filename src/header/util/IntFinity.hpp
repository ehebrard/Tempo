/**
* @author Tim Luchterhand
* @date 03.09.24
* @brief Integer-like class that supports infinity and nan values
*/

#ifndef TEMPO_INTFINITY_HPP
#define TEMPO_INTFINITY_HPP

#include <concepts>
#include <limits>
#include <compare>
#include <ostream>
#include <cmath>
#include <utility>

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

    enum class UnsignedUnderflow {
        ToNan, ToZero, ToInfinity, WrapAround
    };

    /**
     * @brief Integer-like class that supports infinity and nan
     * @details @copybrief
     * Over- / underflows are well defined and produce infinity values. Behaves like ieee745 concerning nan values.
     * @note arithmetic operations on built-in integer types are usually faster. Use this class when extreme performance
     * is not necessary or when arithmetic operations are not the bottleneck.
     * @tparam T underlying integer type
     * @tparam UnderflowMode behavior in case of unsigned underflow (default: produce nan)
     */
    template<std::integral T, UnsignedUnderflow UnderflowMode = UnsignedUnderflow::ToNan>
    class intfinity {
        T value;
        static constexpr bool TwoCompl = not std::is_signed_v<T> or T(-1) == compl T(0);
        static_assert(TwoCompl, "This implementation only works on machines that use the 2's complement for signed"
                                "integers. I suggest using a machine from past WW2 :D");
        static constexpr T NotANumber = detail::conditional_v<T, std::is_signed_v<T>,
                std::numeric_limits<T>::min(), std::numeric_limits<T>::max()>;
        static constexpr T Infinity =
                std::numeric_limits<T>::max() - detail::conditional_v<T, std::is_signed_v<T>, 0, 1>;

        constexpr intfinity(T value, int) noexcept: value(value) {}

        static constexpr T underflow() noexcept {
            static_assert(not std::is_signed_v<T>);
            using enum UnsignedUnderflow;
            if constexpr (UnderflowMode == ToNan) {
                return NotANumber;
            } else if (UnderflowMode == ToInfinity) {
                return Infinity;
            } else if constexpr (UnderflowMode == WrapAround) {
                return Infinity - 1;
            } else {
                return 0;
            }
        }

        template<std::integral I>
        static constexpr T convert(I i) noexcept {
            if constexpr (std::is_signed_v<T> == std::is_signed_v<I> and sizeof(T) >= sizeof(I)) {
                return static_cast<T>(i);
            } else {
                if (std::cmp_greater(i, Infinity)) {
                    return Infinity;
                }

                if constexpr (std::is_signed_v<T>) {
                    if (std::cmp_less(i, -Infinity)) {
                        return -Infinity;
                    }
                } else {
                    if (std::cmp_less(i, 0)) {
                        return underflow();
                    }
                }

                return static_cast<T>(i);
            }
        }

    public:
        /**
         * Positive infinity
         * @return
         */
        static constexpr intfinity Inf() { return {Infinity, 0}; }

        /**
         * Not a number
         * @return
         */
        static constexpr intfinity Nan() { return {NotANumber, 0}; }

        /**
         * Default Ctor. Initializes to UNDEFINED value
         */
        constexpr intfinity() noexcept = default;

        /**
         * Ctor. Converts to infinity if number does not fit into underlying type
         * @param value init value
         */
        constexpr intfinity(T value) noexcept: value(value) {
            if (value == NotANumber) {
                if constexpr (std::is_signed_v<T>) {
                    this->value = detail::sgn(value) * Infinity;
                } else {
                    this->value = Infinity;
                }
            }
        }

        /**
         * Conversion Ctor. Converts to infinity if number does not fit into underlying type
         * @tparam I source integer type
         * @param value init value
         */
        template<std::integral I>
        constexpr intfinity(I value) noexcept: value(convert(value)) {}

        /**
         * Explicit conversion Ctor overload. Converts to infinity if number does not fit into underlying type.
         * @note behaves like cast()
         * @tparam I source integer type
         * @param value init value
         */
        template<std::integral I>
        explicit constexpr intfinity(intfinity<I> value) noexcept: value(convert(value.get())) {}

        /**
         * Float conversion Ctor. Converts to infinity if number does not fit into underlying type. Keeps nan values
         * @tparam F floating point type
         * @param value init value
         */
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
                    this->value = underflow();
                } else if (value > FInf) {
                    this->value = Infinity;
                }
            }
        }

        /**
         * Checks whether value is infinite
         * @return true on positive or negative infinity, false otherwise
         */
        [[nodiscard]] constexpr bool isInf() const noexcept {
            if constexpr (std::is_signed_v<T>) {
                return std::abs(value) == Infinity;
            } else {
                return value == Infinity;
            }
        }

        /**
         * Checks whether value is not a number
         * @return true if value is NaN, false otherwise
         */
        [[nodiscard]] constexpr bool isNan() const noexcept {
            return value == NotANumber;
        }

        /**
         * Equality comparison.
         * @note two nan values are always unequal
         * @param other
         * @return true if values are equal, false otherwise
         */
        constexpr bool operator==(intfinity other) const noexcept {
            if (isNan() or other.isNan()) { return false; }
            return value == other.value;
        }

        /**
         * Equality comparison with float. Performs conversion to float type and then does comparison
         * @tparam F floating point type
         * @param other
         * @return true if values are equal, false otherwise
         */
        template<std::floating_point F>
        constexpr bool operator==(F other) const noexcept {
            return static_cast<F>(*this) == other;
        }

        /**
         * Relational operators
         * @param other
         * @return
         */
        constexpr std::partial_ordering operator<=>(intfinity other) const noexcept {
            if (isNan() or other.isNan()) { return std::partial_ordering::unordered; }
            return value <=> other.value;
        }

        /**
         * Relational operators with float. Performs conversion to float and then does comparison.
         * @tparam F floating point type
         * @param other
         * @return
         */
        template<std::floating_point F>
        constexpr auto operator<=>(F other) const noexcept {
            return static_cast<F>(*this) <=> other;
        }

        /**
         * Compound addition. Performs bounds checking and supports NaN
         * @param other
         * @return
         */
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

        /**
         * Compound addition with float. Performs conversion to float, then performs addition, then converts back to
         * int. Performs bounds checking and supports NaN
         * @tparam F floating point type
         * @param other
         * @return
         */
        template<std::floating_point F>
        constexpr intfinity &operator+=(F other) noexcept {
            *this = static_cast<F>(*this) + other;
            return *this;
        }

        /**
         * Compound subtraction. Performs bounds checking and supports NaN
         * @param other
         * @return
         */
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
                    value = underflow();
                    return *this;
                }
            }

            value -= other.value;
            return *this;
        }

        /**
         * Compound subtraction with float. Performs conversion to float, then performs subtraction, then converts back
         * to int. Performs bounds checking and supports NaN
         * @tparam F floating point type
         * @param other
         * @return
         */
        template<std::floating_point F>
        constexpr intfinity &operator-=(F other) {
            *this = static_cast<F>(*this) - other;
            return *this;
        }

        /**
         * Binary addition. See compound addition
         * @param lhs
         * @param rhs
         * @return
         */
        constexpr friend intfinity operator+(intfinity lhs, intfinity rhs) noexcept {
            lhs += rhs;
            return lhs;
        }

        /**
         * Binary addition with float. See compound addition
         * @param lhs
         * @param rhs
         * @return
         */
        template<std::floating_point F>
        friend constexpr F operator+(intfinity lhs, F rhs) noexcept {
            return static_cast<F>(lhs) + rhs;
        }

        /**
         * Binary addition with float. See compound addition
         * @param lhs
         * @param rhs
         * @return
         */
        template<std::floating_point F>
        friend constexpr F operator+(F lhs, intfinity rhs) noexcept {
            return lhs + static_cast<F>(rhs);
        }

        /**
         * Binary subtraction. See compound addition
         * @param lhs
         * @param rhs
         * @return
         */
        constexpr friend intfinity operator-(intfinity lhs, intfinity rhs) noexcept {
            lhs -= rhs;
            return lhs;
        }

        /**
         * Binary subtraction with float. See compound subtraction
         * @param lhs
         * @param rhs
         * @return
         */
        template<std::floating_point F>
        friend constexpr F operator-(intfinity lhs, F rhs) noexcept {
            return static_cast<F>(lhs) - rhs;
        }

        /**
         * Binary subtraction with float. See compound subtraction
         * @param lhs
         * @param rhs
         * @return
         */
        template<std::floating_point F>
        friend constexpr F operator-(F lhs, intfinity rhs) noexcept {
            return lhs - static_cast<F>(rhs);
        }


        /**
         * prefix increment
         * @return
         */
        constexpr intfinity &operator++() noexcept {
            return operator+=(1);
        }

        /**
         * postfix increment
         * @return
         */
        constexpr intfinity operator++(int) noexcept {
            auto old = *this;
            operator++();
            return old;
        }

        /**
         * prefix decrement
         * @return
         */
        constexpr intfinity &operator--() noexcept {
            return operator-=(1);
        }

        /**
         * postfix decrement
         * @return
         */
        constexpr intfinity operator--(int) noexcept {
            auto old = *this;
            operator--();
            return old;
        }


        /**
         * Compound multiplication. Performs bounds checking and supports NaN
         * @param other
         * @return
         */
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

        /**
         * Compound multiplication with float. Performs conversion to float, then performs multiplication, then converts
         * back to int. Performs bounds checking and supports NaN
         * @tparam F floating point type
         * @param other
         * @return
         */
        template<std::floating_point F>
        constexpr intfinity &operator*=(F other) noexcept {
            *this = static_cast<F>(*this) * other;
            return *this;
        }

        /**
         * Compound division. Performs bounds checking and supports NaN
         * @param other
         * @return
         */
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

        /**
         * Compound division with float. Performs conversion to float, then performs division, then converts
         * back to int. Performs bounds checking and supports NaN
         * @tparam F floating point type
         * @param other
         * @return
         */
        template<std::floating_point F>
        constexpr intfinity &operator/=(F other) noexcept {
            *this = static_cast<F>(*this) / other;
            return *this;
        }

        /**
         * Binary multiplication. See compound multiplication
         * @param lhs
         * @param rhs
         * @return
         */
        constexpr friend intfinity operator*(intfinity lhs, intfinity rhs) noexcept {
            lhs *= rhs;
            return lhs;
        }

        /**
         * Binary multiplication with float. See compound multiplication
         * @param lhs
         * @param rhs
         * @return
         */
        template<std::floating_point F>
        friend constexpr F operator*(intfinity lhs, F rhs) noexcept {
            return static_cast<F>(lhs) * rhs;
        }

        /**
         * Binary multiplication with float. See compound multiplication
         * @param lhs
         * @param rhs
         * @return
         */
        template<std::floating_point F>
        friend constexpr F operator*(F lhs, intfinity rhs) noexcept {
            return lhs * static_cast<F>(rhs);
        }

        /**
         * Binary division. See compound division
         * @param lhs
         * @param rhs
         * @return
         */
        constexpr friend intfinity operator/(intfinity lhs, intfinity rhs) noexcept {
            lhs /= rhs;
            return lhs;
        }

        /**
         * Binary division with float. See compound division
         * @param lhs
         * @param rhs
         * @return
         */
        template<std::floating_point F>
        friend constexpr F operator/(intfinity lhs, F rhs) noexcept {
            return static_cast<F>(lhs) / rhs;
        }

        /**
         * Binary division with float. See compound division
         * @param lhs
         * @param rhs
         * @return
         */
        template<std::floating_point F>
        friend constexpr F operator/(F lhs, intfinity rhs) noexcept {
            return lhs / static_cast<F>(rhs);
        }

        /**
         * Unary +. Returns *this
         * @return *this
         */
        constexpr intfinity operator+() const noexcept {
            return *this;
        }

        /**
         * Unary -. Inverts sign on signed types
         * @return
         */
        template<std::signed_integral = T>
        constexpr intfinity operator-() const noexcept {
            if (not isNan()) {
                return -value;
            }

            return *this;
        }

        /**
         * Explicit float conversion. If the value is too large for the destiantion type, the behavior is undefined
         * @tparam F floating point type
         * @return
         */
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

        friend std::istream &operator>>(std::istream &is, intfinity &val) {
            T data;
            is >> data;
            val = data;
            return is;
        }

        /**
         * Gets the raw integer value
         * @return
         */
        constexpr T get() const noexcept {
            return value;
        }

        /**
         * Cast to different integer type. Yields infinity if too large for target type. Supports nan
         * @tparam I target type
         * @tparam UMode target underflow mode
         * @return cast result
         */
        template<std::integral I, UnsignedUnderflow UMode = UnsignedUnderflow::ToNan>
        constexpr auto cast() const noexcept {
            return intfinity<I, UMode>(value);
        }
    };
}
namespace std {
    template<typename T, tempo::UnsignedUnderflow U>
    class numeric_limits<tempo::intfinity<T, U>> {
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

        static constexpr tempo::intfinity<T, U> min() noexcept {
            if constexpr (std::is_signed_v<T>) {
                return std::numeric_limits<T>::min() + 2;
            } else {
                return 0;
            }
        }

        static constexpr tempo::intfinity<T, U> max() noexcept {
            return std::numeric_limits<T>::max() - tempo::detail::conditional_v<T, std::is_signed_v<T>, 1, 2>;
        }

        static constexpr tempo::intfinity<T, U> lowest() noexcept {
            return min();
        }

        static constexpr tempo::intfinity<T, U> epsilon() noexcept {
            return 0;
        }

        static constexpr tempo::intfinity<T, U> round_error() noexcept {
            return 0;
        }

        static constexpr tempo::intfinity<T, U> denorm_min() noexcept {
            return 0;
        }

        static constexpr auto infinity() noexcept {
            return tempo::intfinity<T, U>::Inf();
        }

        static constexpr auto quiet_NaN() noexcept {
            return tempo::intfinity<T, U>::Nan();
        }
    };

    template<typename T, tempo::UnsignedUnderflow U>
    constexpr bool isinf(tempo::intfinity<T, U> val) noexcept {
        return val.isInf();
    }

    template<typename T, tempo::UnsignedUnderflow U>
    constexpr bool isfinite(tempo::intfinity<T, U> val) noexcept {
        return not val.isInf() and not val.isNan();
    }

    template<typename T, tempo::UnsignedUnderflow U>
    constexpr bool isnan(tempo::intfinity<T, U> val) noexcept {
        return val.isNan();
    }
}

#endif //TEMPO_INTFINITY_HPP
