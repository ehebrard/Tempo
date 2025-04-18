/**
* @author Tim Luchterhand
* @date 04.09.24
* @brief
*/

#include <gtest/gtest.h>

#include "testing.hpp"
#include "util/IntFinity.hpp"

using namespace tempo;

template<typename T>
void testRelationsNan() {
    auto nan = intfinity<T>::Nan();
    auto inf = intfinity<T>::Inf();
    EXPECT_FALSE(nan == 3);
    EXPECT_TRUE(nan != 3);
    EXPECT_FALSE(nan <= 3);
    EXPECT_FALSE(nan >= 3);
    EXPECT_FALSE(nan < 3);
    EXPECT_FALSE(nan > 3);

    EXPECT_FALSE(nan == inf);
    EXPECT_TRUE(nan != inf);
    EXPECT_FALSE(nan <= inf);
    EXPECT_FALSE(nan >= inf);
    EXPECT_FALSE(nan < inf);
    EXPECT_FALSE(nan > inf);

    EXPECT_FALSE(nan == nan);
    EXPECT_TRUE(nan != nan);
    EXPECT_FALSE(nan <= nan);
    EXPECT_FALSE(nan >= nan);
    EXPECT_FALSE(nan < nan);
    EXPECT_FALSE(nan > nan);
}

template<typename T>
void testRelationsInf() {
    auto inf = intfinity<T>::Inf();
    EXPECT_TRUE(inf > 1000);
    EXPECT_TRUE(inf >= 1000);
    EXPECT_TRUE(inf != 1000);
    EXPECT_FALSE(inf == 1000);
    EXPECT_FALSE(inf <= 1000);
    EXPECT_FALSE(inf < 1000);

    EXPECT_FALSE(inf > inf);
    EXPECT_FALSE(inf != inf);
    EXPECT_FALSE(inf < inf);
    EXPECT_TRUE(inf == inf);
    EXPECT_TRUE(inf >= inf);
    EXPECT_TRUE(inf <= inf);
}

template<typename T>
void testRelationsNanFloat() {
    auto nan = intfinity<T>::Nan();
    auto nanf = std::numeric_limits<float>::quiet_NaN();
    intfinity<T> number = 5;
    EXPECT_FALSE(nan == 3.3);
    EXPECT_TRUE(nan != 3.f);
    EXPECT_FALSE(nan <= 3.f);
    EXPECT_FALSE(nan >= 3.5);
    EXPECT_FALSE(nan < 3.f);
    EXPECT_FALSE(nan > 3.f);
    EXPECT_TRUE(number != nanf);
    EXPECT_FALSE(number == nanf);
    EXPECT_FALSE(number <= nanf);
    EXPECT_FALSE(number >= nanf);
    EXPECT_FALSE(number < nanf);
    EXPECT_FALSE(number > nanf);
}

template<typename T>
void testRelationsInfFloat() {
    auto inf = intfinity<T>::Inf();
    auto inff = -std::numeric_limits<double>::infinity();
    intfinity<T> number = 4;
    EXPECT_TRUE(inf > 1000.f);
    EXPECT_TRUE(inf >= 1000.f);
    EXPECT_TRUE(inf != 1000.4);
    EXPECT_FALSE(inf == 1000.);
    EXPECT_FALSE(inf <= 1000.f);
    EXPECT_FALSE(inf < 1000.);

    EXPECT_FALSE(number == inff);
    EXPECT_FALSE(number <= inff);
    EXPECT_FALSE(number < inff);
    EXPECT_TRUE(number != inff);
    EXPECT_TRUE(number >= inff);
    EXPECT_TRUE(number > inff);
}

template<typename T, bool Sign>
void testArithmeticsInf() {
    auto inf = intfinity<T>::Inf();
    intfinity<T> zero = 0;
    intfinity<T> three = 3;
    if constexpr (not Sign) {
        inf = -intfinity<T>::Inf();
        three = -3;
    }

    EXPECT_EQ(inf += 1000000, inf);
    EXPECT_EQ(inf -= 1000000, inf);
    EXPECT_EQ(inf *= 1000000, inf);
    EXPECT_EQ(inf /= 1000000, inf);
    EXPECT_EQ(inf + 1000000, inf);
    EXPECT_EQ(inf - 1000000, inf);
    EXPECT_EQ(inf * 1000000, inf);
    EXPECT_EQ(inf / 1000000, inf);

    EXPECT_TRUE((inf - inf).isNan());
    EXPECT_TRUE((inf / inf).isNan());
    EXPECT_TRUE((inf * 0).isNan());
    EXPECT_TRUE((three / 0).isInf());
    EXPECT_TRUE((zero / 0).isNan());
    EXPECT_EQ(+inf, inf);
    if constexpr(std::is_signed_v<T>) {
        EXPECT_TRUE((-inf).isInf());
        EXPECT_EQ(-inf, -inf.get());
    }
}

template<typename T>
void testArithmeticsNan() {
    auto nan = intfinity<T>::Nan();
    auto inf = intfinity<T>::Inf();
    EXPECT_TRUE((nan += 1000000).isNan());
    EXPECT_TRUE((nan -= 1000000).isNan());
    EXPECT_TRUE((nan *= 1000000).isNan());
    EXPECT_TRUE((nan /= 1000000).isNan());
    EXPECT_TRUE((nan + 1000000).isNan());
    EXPECT_TRUE((nan - 1000000).isNan());
    EXPECT_TRUE((nan * 1000000).isNan());
    EXPECT_TRUE((nan / 1000000).isNan());

    EXPECT_TRUE((nan += inf).isNan());
    EXPECT_TRUE((nan -= inf).isNan());
    EXPECT_TRUE((nan *= inf).isNan());
    EXPECT_TRUE((nan /= inf).isNan());
    EXPECT_TRUE((nan + inf).isNan());
    EXPECT_TRUE((nan - inf).isNan());
    EXPECT_TRUE((nan * inf).isNan());
    EXPECT_TRUE((nan / inf).isNan());

    EXPECT_TRUE((nan += nan).isNan());
    EXPECT_TRUE((nan -= nan).isNan());
    EXPECT_TRUE((nan *= nan).isNan());
    EXPECT_TRUE((nan /= nan).isNan());
    EXPECT_TRUE((nan + nan).isNan());
    EXPECT_TRUE((nan - nan).isNan());
    EXPECT_TRUE((nan * nan).isNan());
    EXPECT_TRUE((nan / nan).isNan());
    EXPECT_TRUE((+nan).isNan());
    if constexpr(std::is_signed_v<T>) {
        EXPECT_TRUE((-nan).isNan());
    }
}

template<typename T>
void testPlusOverflow() {
    auto number = std::numeric_limits<intfinity<T>>::max();
    number += 0;
    EXPECT_FALSE(number.isInf());
    EXPECT_FALSE(number.isNan());
    auto smaller = number - 100;
    EXPECT_TRUE((number + 1).isInf());
    EXPECT_FALSE((smaller + 80).isInf());
    EXPECT_TRUE((smaller + 80000).isInf());
    EXPECT_TRUE((intfinity<T>{} + intfinity<T>::Inf()).isInf());
}

template<typename T>
void testMultNormal() {
    intfinity<T> number = 5;
    EXPECT_EQ(number * 1, number);
    EXPECT_EQ(number * 3, 15);
    EXPECT_EQ(number *= 2, 10);
    EXPECT_EQ(number *= 0, 0);
    EXPECT_EQ(number * number, 0);
}

template<typename T>
void testAddFloat() {
    intfinity<T> number = 5;
    EXPECT_EQ(number += 0.3, 5);
    EXPECT_EQ(number + 0.3f, 5.3f);
    EXPECT_EQ(0.3f + number, 5.3f);
    EXPECT_TRUE((std::same_as<float, decltype(number + 0.1f)>));
}

template<typename T>
void testAddFloatSpecial() {
    constexpr auto Inff = std::numeric_limits<double>::infinity();
    auto large = std::numeric_limits<intfinity<T>>::max();
    large += 0.7;
    EXPECT_FALSE(large.isInf());
    large += 1.03;
    EXPECT_EQ(large, intfinity<T>::Inf());
    EXPECT_EQ(large + 1000.f, intfinity<T>::Inf());
    EXPECT_TRUE(std::isinf(large + 1000.f));
    if constexpr (std::is_signed_v<T>) {
        EXPECT_TRUE(std::isnan(-intfinity<T>::Inf() + std::numeric_limits<double>::infinity()));
    }

    intfinity<T> number = 5;
    EXPECT_EQ(number + Inff, Inff);
    EXPECT_TRUE(std::isnan(number + std::numeric_limits<float>::signaling_NaN()));
}

template<typename T>
void testSubFloat() {
    intfinity<T> number = 5;
    EXPECT_EQ(number -= 0.3, 4);
    EXPECT_EQ(number - 0.3f, 3.7f);
    EXPECT_EQ(0.3f - number, -3.7f);
    EXPECT_TRUE((std::same_as<float, decltype(number - 0.1f)>));
}

template<typename T>
void testMultFloat() {
    intfinity<T> number = 5;
    EXPECT_EQ(number *= 0.3, 1);
    EXPECT_EQ(number * 0.3f, 0.3f);
    EXPECT_EQ(0.3f * number, 0.3f);
    EXPECT_TRUE((std::same_as<float, decltype(number * 0.1f)>));
}

template<typename T>
void testDivFloat() {
    intfinity<T> number = 5;
    EXPECT_EQ(number /= 2.5, 2);
    EXPECT_EQ(number / 0.5f, 4.f);
    EXPECT_EQ(0.5f / number, 0.25f);
    EXPECT_TRUE((std::same_as<float, decltype(number / 0.1f)>));
}

TEST(util, intfinity_default) {
    intfinity<int> zero{};
    EXPECT_EQ(zero.get(), 0);
    EXPECT_EQ(zero, 0);
    EXPECT_FALSE(zero.isInf());
    EXPECT_FALSE(zero.isNan());
}

TEST(util, intfinity_ctor) {
    intfinity<int> number = 5;
    EXPECT_EQ(number.get(), 5);
    EXPECT_EQ(number, 5);
    EXPECT_FALSE(number.isInf());
    EXPECT_FALSE(number.isNan());
}

TEST(util, intfinity_ctor_overflow) {
    intfinity<int> inf = std::numeric_limits<int>::min();
    EXPECT_EQ(inf, -intfinity<int>::Inf());
}

TEST(util, intfinity_ctor_conversion) {
    intfinity<int> inf = std::numeric_limits<long>::max();
    EXPECT_TRUE(inf.isInf());
    intfinity<unsigned> nan = -4;
    EXPECT_TRUE(nan.isNan());
    intfinity<unsigned, UnsignedUnderflow::ToZero> zero = -3;
    EXPECT_EQ(zero, 0);
    intfinity<unsigned, UnsignedUnderflow::ToInfinity> uinf = -1;
    EXPECT_TRUE(uinf.isInf());
    intfinity<short> sInf = std::numeric_limits<std::int32_t>::min();
    EXPECT_EQ(sInf, -intfinity<short>::Inf());
}

TEST(util, intfinity_unsigned_underflow) {
    using enum UnsignedUnderflow;
    intfinity<unsigned, ToZero> z = -4;
    intfinity<unsigned, ToNan> n = -4;
    intfinity<unsigned, ToInfinity> i = -4;
    intfinity<unsigned, WrapAround> w = -4;
    EXPECT_EQ(z, 0);
    EXPECT_TRUE(n.isNan());
    EXPECT_TRUE(i.isInf());
    EXPECT_EQ(w, std::numeric_limits<decltype(w)>::max());
}

TEST(util, intfinity_float_conversion) {
    intfinity<int> number = 1.5;
    intfinity<unsigned, tempo::UnsignedUnderflow::ToZero> uNumber = 1.5;
    intfinity<unsigned, tempo::UnsignedUnderflow::ToInfinity> uONumber = 1.5;
    EXPECT_EQ(number, 1);
    EXPECT_EQ(uNumber, 1);
    number = -14.2f;
    uNumber = -14.2f;
    uONumber = -14.2f;
    EXPECT_EQ(number, -14);
    EXPECT_EQ(uNumber, 0);
    EXPECT_TRUE(uONumber.isInf());
}

TEST(util, intfinity_float_conversion_special) {
    intfinity<int> number = std::numeric_limits<double>::infinity();
    intfinity<unsigned, tempo::UnsignedUnderflow::ToZero> uNumber = std::numeric_limits<float>::infinity();
    intfinity<unsigned, tempo::UnsignedUnderflow::ToInfinity> uONumber = -std::numeric_limits<double>::infinity();
    EXPECT_TRUE(number.isInf());
    EXPECT_EQ(number, intfinity<int>::Inf());
    EXPECT_TRUE(uNumber.isInf());
    EXPECT_TRUE(uONumber.isInf());
    number = -1 / 0.0;
    EXPECT_EQ(number, -intfinity<int>::Inf());
    uNumber = -1 / 0.0;
    EXPECT_EQ(uNumber, 0);
    number = std::numeric_limits<float>::quiet_NaN();
    uNumber = std::numeric_limits<double>::quiet_NaN();
    uONumber = std::numeric_limits<double>::quiet_NaN();
    EXPECT_TRUE(number.isNan());
    EXPECT_TRUE(uNumber.isNan());
    EXPECT_TRUE(uONumber.isNan());
    number = std::numeric_limits<float>::signaling_NaN();
    uNumber = std::numeric_limits<double>::signaling_NaN();
    uONumber = std::numeric_limits<double>::signaling_NaN();
    EXPECT_TRUE(number.isNan());
    EXPECT_TRUE(uNumber.isNan());
    EXPECT_TRUE(uONumber.isNan());
    number = static_cast<double>(intfinity<int>::Inf().get());
    uNumber = static_cast<double>(intfinity<unsigned>::Inf().get());
    uONumber = static_cast<double>(intfinity<unsigned, tempo::UnsignedUnderflow::ToInfinity>::Inf().get());
    EXPECT_EQ(number, intfinity<int>::Inf());
    EXPECT_EQ(uNumber, (intfinity<unsigned, tempo::UnsignedUnderflow::ToZero>::Inf()));
    EXPECT_EQ(uONumber, (intfinity<unsigned, tempo::UnsignedUnderflow::ToInfinity>::Inf()));
}

TEST(util, float_conversion_underflow_mode) {
    intfinity<unsigned> uNan = -1.4;
    intfinity<unsigned, tempo::UnsignedUnderflow::ToZero> uZero = -4.3;
    intfinity<unsigned, tempo::UnsignedUnderflow::ToInfinity> uInf = -1.3;
    EXPECT_TRUE(uNan.isNan());
    EXPECT_EQ(uZero, 0);
    EXPECT_TRUE(uInf.isInf());
}

TEST(util, intfinity_to_float) {
    constexpr auto FInf = std::numeric_limits<float>::infinity();
    constexpr auto DInf = std::numeric_limits<double>::infinity();
    intfinity number = -4;
    EXPECT_EQ(static_cast<double>(number), -4.0);
    number = 4;
    EXPECT_EQ(static_cast<float>(number), 4.0f);
    number = intfinity<int>::Inf();
    EXPECT_EQ(static_cast<float>(number), FInf);
    EXPECT_EQ(static_cast<double>(-number), -DInf);
    number = intfinity<int>::Nan();
    EXPECT_TRUE(std::isnan(static_cast<double>(number)));
    EXPECT_TRUE(std::isnan(static_cast<float>(-number)));
}

TEST(util, intfinity_comparisons_normal) {
    intfinity five = 5;
    intfinity minusTwo = -2;
    EXPECT_EQ(five, five);
    EXPECT_FALSE(five < five);
    EXPECT_FALSE(five > five);
    EXPECT_FALSE(five != five);
    EXPECT_TRUE(five <= five);
    EXPECT_TRUE(five >= five);

    EXPECT_TRUE(five != minusTwo);
    EXPECT_FALSE(five == minusTwo);
    EXPECT_TRUE(five > minusTwo);
    EXPECT_TRUE(five >= minusTwo);
    EXPECT_FALSE(five <= minusTwo);
    EXPECT_FALSE(five < minusTwo);
}

TEST(util, intfinity_comparisons_nan) {
    testRelationsNan<int>();
    testRelationsNan<unsigned>();
}

TEST(util, intfinity_comparisons_inf) {
    testRelationsInf<int>();
    testRelationsInf<unsigned>();
}

TEST(util, intfinity_comparisons_float) {
    intfinity number = 5;
    EXPECT_TRUE(number < 5.3f);
    EXPECT_TRUE(number <= 5.3);
    EXPECT_TRUE(number == 5.0);
    EXPECT_TRUE(number >= 5.f);
    EXPECT_TRUE(number > 4.45f);
    EXPECT_TRUE(number != 5.00001);
    testRelationsNanFloat<int>();
    testRelationsNanFloat<unsigned>();
    testRelationsInfFloat<int>();
    testRelationsInfFloat<unsigned>();
}

TEST(util, intfinity_arithmetics_unsigned_inf) {
    testArithmeticsInf<unsigned, true>();
}

TEST(util, intfinity_arithmetics_signed_inf) {
    testArithmeticsInf<int, true>();
    testArithmeticsInf<int, false>();
}

TEST(util, intfinity_arithmetics_nan) {
    testArithmeticsNan<unsigned>();
    testArithmeticsNan<int>();
}

TEST(util, intfinity_arithmetics_unsigned_plus_normal) {
    intfinity number = 4u;
    number += 3;
    EXPECT_EQ(number, 7);
    EXPECT_EQ(number + 3, 10);
    EXPECT_EQ(++number, 8);
    EXPECT_EQ(number++, 8);
    EXPECT_EQ(number, 9);
    EXPECT_EQ(+number, number);
}

TEST(util, intfinity_arithmetics_unsigned_plus_overflow) {
    testPlusOverflow<unsigned>();
}

TEST(util, intfinity_arithmetics_signed_plus_normal) {
    intfinity number = -4;
    number += 3;
    EXPECT_EQ(number, -1);
    EXPECT_EQ(number + 3, 2);
    EXPECT_EQ(++number, 0);
    EXPECT_EQ(number++, 0);
    EXPECT_EQ(number, 1);
    EXPECT_EQ(+number, number);
}

TEST(util, intfinity_arithmetics_signed_plus_overflow) {
    testPlusOverflow<int>();
    auto nInf = -intfinity<int>::Inf();
    EXPECT_TRUE((nInf + intfinity<int>::Inf()).isNan());
    auto low = std::numeric_limits<intfinity<int>>::min();
    EXPECT_EQ(low + (-3), -intfinity<int>::Inf());
    EXPECT_EQ(low += (-3), -intfinity<int>::Inf());
}


TEST(util, intfinity_arithmetics_plus_float_normal) {
    testAddFloat<int>();
    testAddFloat<unsigned>();
}

TEST(util, intfinity_arithmetics_plus_float_special) {
    testAddFloatSpecial<int>();
    testAddFloatSpecial<unsigned>();
}

TEST(util, intfinity_arithmetics_unsigned_minus_normal) {
    intfinity number = 15u;
    number -= 5;
    EXPECT_EQ(number, 10);
    EXPECT_EQ(number - 10, 0);
    EXPECT_EQ(--number, 9);
    EXPECT_EQ(number--, 9);
    EXPECT_EQ(number, 8);
}

TEST(util, intfinity_arithmetics_unsigned_minus_underflow) {
    intfinity<unsigned, tempo::UnsignedUnderflow::ToZero> number = 5u;
    EXPECT_EQ(number -8, 0);
    EXPECT_EQ(number -= 7, 0);
    EXPECT_EQ(number + 5, 5);

    intfinity<unsigned, tempo::UnsignedUnderflow::ToInfinity> uFlowNumber(5);
    EXPECT_TRUE((uFlowNumber - 8).isInf());
    EXPECT_TRUE((uFlowNumber -= 8).isInf());
}

TEST(util, intfinity_arithmetics_signed_minus_normal) {
    intfinity number = 5;
    number -= 5;
    EXPECT_EQ(number, 0);
    EXPECT_EQ(number - 10, -10);
    EXPECT_EQ(--number, -1);
    EXPECT_EQ(number--, -1);
    EXPECT_EQ(number, -2);
    EXPECT_EQ(-number, 2);
}

TEST(util, intfinity_arithmetics_signed_minus_underflow) {
    auto low = std::numeric_limits<intfinity<int>>::min();
    EXPECT_TRUE((low - 1).isInf());
    EXPECT_EQ(low - 1, -intfinity<int>::Inf());
    EXPECT_TRUE((low -= 1000).isInf());
    EXPECT_EQ(-low, intfinity<int>::Inf());
    auto high = std::numeric_limits<intfinity<int>>::max();
    EXPECT_EQ(high - (-5), intfinity<int>::Inf());
}

TEST(util, intfinity_arithmetics_minus_float_normal) {
    testSubFloat<int>();
    testSubFloat<unsigned>();
}

TEST(util, intfinity_arithmetics_minus_float_special) {
    auto number = std::numeric_limits<intfinity<int>>::min();
    EXPECT_EQ(number -= 0.7, number);
    EXPECT_EQ(number -= 1.4, -intfinity<int>::Inf());
    EXPECT_TRUE(std::isinf(number - 10));
    intfinity<unsigned, tempo::UnsignedUnderflow::ToZero> uNumber = 4u;
    EXPECT_DOUBLE_EQ(uNumber - 5.6, -1.6);
    EXPECT_EQ(uNumber -= 5.6, 0);
    intfinity<unsigned, tempo::UnsignedUnderflow::ToInfinity> uONumber = 4;
    EXPECT_EQ(uONumber -= 4.6, (intfinity<unsigned, tempo::UnsignedUnderflow::ToInfinity>::Inf()));
    EXPECT_TRUE(std::isnan(intfinity<int>::Inf() - std::numeric_limits<float>::infinity()));
    EXPECT_TRUE(std::isnan(std::numeric_limits<float>::signaling_NaN() - uONumber));
}

TEST(util, intfinity_arithmetics_unsigned_mult_normal) {
    testMultNormal<unsigned>();
}

TEST(util, intfinity_arithmetics_unsigned_mult_overflow) {
    auto number = std::numeric_limits<intfinity<unsigned>>::max() / 2;
    EXPECT_TRUE((number * 3).isInf());
}

TEST(util, intfinity_arithmetics_signed_mult_normal) {
    testMultNormal<int>();
    intfinity number = 5;
    EXPECT_EQ(number * -3, -15);
    EXPECT_EQ(number *= -2, -10);
    EXPECT_EQ(number * -4, 40);
}

TEST(util, intfinity_arithmetics_signed_mult_overflow) {
    auto number = std::numeric_limits<intfinity<int>>::max() / 2;
    EXPECT_TRUE((number * 3).isInf());
    EXPECT_TRUE((number * (-3)).isInf());
    EXPECT_EQ(number * (-3), -intfinity<int>::Inf());
}

TEST(util, intfinity_arithmetics_signed_mult_float_normal) {
    testMultFloat<int>();
    testMultFloat<unsigned>();
}

TEST(util, intfinity_arithmetics_signed_mult_float_special) {
    intfinity number = 1000;
    EXPECT_EQ(number *= -1e10, -intfinity<int>::Inf());
    EXPECT_EQ(0.3 * number * -0.4, std::numeric_limits<double>::infinity());
    intfinity<unsigned, tempo::UnsignedUnderflow::ToZero> uNumber = 1000u;
    EXPECT_EQ(uNumber *= -1e10, 0);
    EXPECT_TRUE(std::isnan(uNumber * 0.3 * number));
}


TEST(util, intfinity_arithmetics_unsigned_div) {
    intfinity number = 18u;
    EXPECT_EQ(number /= 4, 4);
    EXPECT_EQ(number / 2, 2);
    EXPECT_EQ(number / 0, intfinity<unsigned>::Inf());
}

TEST(util, intfinity_arithmetics_signed_div) {
    intfinity number = 18;
    EXPECT_EQ(number /= 4, 4);
    EXPECT_EQ(number / -2, -2);
    EXPECT_EQ(-number / 0, -intfinity<int>::Inf());
    EXPECT_EQ(number / 0, intfinity<int>::Inf());
}

TEST(util, intfinity_arithmetics_signed_div_float_normal) {
    testDivFloat<int>();
    testDivFloat<unsigned>();
}

TEST(util, intfinity_arithmetics_signed_div_float_special) {
    intfinity number = 5;
    EXPECT_EQ(number / 0.0, std::numeric_limits<double>::infinity());
    number /= 1000;
    EXPECT_EQ(0.4f / number, std::numeric_limits<float>::infinity());
    EXPECT_TRUE((number /= 0.0).isNan());
    EXPECT_TRUE(std::isnan(number / 2.0));
}

template<typename T, UnsignedUnderflow B>
void testNumericLimits() {
    using L = std::numeric_limits<intfinity<T, B>>;
    EXPECT_TRUE(L::is_specialist);
    EXPECT_EQ(L::is_signed, std::is_signed_v<T>);
    EXPECT_TRUE(L::is_integer);
    EXPECT_TRUE(L::is_exact);
    EXPECT_TRUE(L::has_infinity);
    EXPECT_TRUE(L::has_quiet_NaN);
    EXPECT_FALSE(L::has_signaling_NaN);
    EXPECT_EQ(L::has_denorm, std::denorm_absent);
    EXPECT_FALSE(L::has_denorm_loss);
    EXPECT_EQ(L::round_style, std::round_toward_zero);
    EXPECT_FALSE(L::is_iec559);
    EXPECT_TRUE(L::is_bounded);
    EXPECT_FALSE(L::is_modulo);
    EXPECT_FALSE(L::traps);
    EXPECT_FALSE(L::tinyness_before);
    EXPECT_EQ(L::infinity(), (intfinity<T, B>::Inf()));
    EXPECT_TRUE(L::quiet_NaN().isNan());
    EXPECT_EQ(L::epsilon(), 0);
    EXPECT_EQ(L::denorm_min(), 0);
    EXPECT_EQ(L::round_error(), 0);
    auto max = std::numeric_limits<intfinity<T, B>>::max();
    EXPECT_FALSE(max.isInf());
    EXPECT_EQ(max + 1, (intfinity<T, B>::Inf()));
    EXPECT_EQ(L::lowest(), L::min());
    auto min = std::numeric_limits<intfinity<T, B>>::min();
    EXPECT_FALSE(min.isInf());
    if constexpr (std::is_signed_v<T>) {
        EXPECT_EQ(min - 1, (-intfinity<T, B>::Inf()));
    } else {
        EXPECT_EQ(min, 0);
    }
}

TEST(util, intfinity_numeric_limits) {
    using enum tempo::UnsignedUnderflow;
    testNumericLimits<int, ToNan>();
    testNumericLimits<unsigned, ToNan>();
    testNumericLimits<unsigned, ToInfinity>();
    testNumericLimits<unsigned, ToZero>();
}

TEST(util, intfinity_std) {
    auto nan = intfinity<int>::Nan();
    auto inf = intfinity<int>::Inf();
    intfinity number = 4;
    EXPECT_TRUE(std::isnan(nan));
    EXPECT_FALSE(std::isinf(nan));
    EXPECT_FALSE(std::isfinite(nan));
    EXPECT_FALSE(std::isnan(inf));
    EXPECT_TRUE(std::isinf(inf));
    EXPECT_FALSE(std::isfinite(inf));
    EXPECT_FALSE(std::isnan(number));
    EXPECT_FALSE(std::isinf(number));
    EXPECT_TRUE(std::isfinite(number));
}

TEST(util, intfinity_cast) {
    intfinity number = 17;
    auto sNumber = number.cast<short>();
    EXPECT_EQ(number.get(), sNumber.get());
    number = std::numeric_limits<short>::min();
    sNumber = number.cast<short>();
    EXPECT_EQ(sNumber, -intfinity<short>::Inf());
    number = std::numeric_limits<int>::max() - 2;
    auto lNumber = number.cast<long>();
    EXPECT_EQ(number, lNumber.cast<int>());
    EXPECT_EQ(number.get(), lNumber.get());
    number = -17;
    auto uNumber = number.cast<unsigned>();
    EXPECT_TRUE(uNumber.isNan());
    auto uzNumber = number.cast<unsigned, UnsignedUnderflow::ToZero>();
    EXPECT_EQ(uzNumber, 0);
}

template<typename T, typename U, UnsignedUnderflow UM>
void testCast() {
    using namespace tempo::testing;
    intfinity src = random_int<T>(std::numeric_limits<T>::min(), std::numeric_limits<T>::max());
    auto dest1 = src.template cast<U, UM>();
    auto dest2 = static_cast<intfinity<U, UM>>(src);
    EXPECT_TRUE(dest1 == dest2 or (dest1.isNan() and dest2.isNan()));
}

TEST(util, intfinity_cast2) {
    using enum UnsignedUnderflow;
    for (int i = 0; i < 10000; ++i) {
        testCast<int, unsigned, ToZero>();
        testCast<unsigned long, short, ToNan>();
        testCast<int, unsigned, ToNan>();
        testCast<long, unsigned, ToInfinity>();
        testCast<long long, long, ToNan>();
        testCast<short, unsigned, WrapAround>();
    }
}