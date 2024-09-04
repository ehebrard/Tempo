/**
* @author Tim Luchterhand
* @date 04.09.24
* @brief
*/

#include <gtest/gtest.h>

#include "util/IntFinity.hpp"

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

TEST(util, intfinity_default) {
    intfinity<int> zero;
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
    intfinity<unsigned> inf = std::numeric_limits<double>::quiet_NaN();
    int infi = std::numeric_limits<double>::infinity();
    std::cout << inf;
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
    auto nan = intfinity<int>::Nan();
    auto inf = intfinity<int>::Inf();
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

TEST(util, intfinity_comparisons_inf) {
    auto inf = intfinity<int>::Inf();
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
    intfinity number = 5u;
    EXPECT_EQ(number -8, 0);
    EXPECT_EQ(number -= 7, 0);
    EXPECT_EQ(number + 5, 5);

    intfinity<unsigned, true> uFlowNumber(5);
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
