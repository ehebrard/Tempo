/**
* @author Tim Luchterhand
* @date 09.01.25
* @file Lookup.cpp
* @brief
*/

#include <gtest/gtest.h>
#include <Iterators.hpp>

#include "util/Lookup.hpp"
#include "Literal.hpp"

TEST(util, Lookup_default_ctor) {
    using namespace tempo;
    Lookup<int, int> lookup;
    EXPECT_EQ(lookup.size(), 0);
    EXPECT_EQ(lookup.keyOffset(), 0);
    EXPECT_THROW(lookup.at(0), std::out_of_range);
}

TEST(util, Lookup_ctor) {
    using namespace tempo;
    auto keys = {3, 1, 7};
    std::vector values{9, 2, 1};
    Lookup lookup(keys, 0, values);
    EXPECT_EQ(lookup.keyOffset(), 1);
    EXPECT_EQ(lookup.maxKey(), 7);
    EXPECT_EQ(lookup.size(), 7);
    EXPECT_EQ(lookup.data(), (std::vector{2, 0, 9, 0, 0, 0, 1}));
}

TEST(util, Lookup_ctor1) {
    using namespace tempo;
    auto keys = {3, 1, 7};
    Lookup lookup(keys, 17);
    EXPECT_EQ(lookup.keyOffset(), 1);
    EXPECT_EQ(lookup.maxKey(), 7);
    EXPECT_EQ(lookup.size(), 7);
    for (auto v : lookup.data()) {
        EXPECT_EQ(v, 17);
    }
}

template<typename L, typename Keys, typename Values, typename Key>
void testAccess(L &lookup, const Keys &keys, const Values &values, Key k1, Key k2) {
    for (auto [k, v] : iterators::const_zip(keys, values)) {
        EXPECT_EQ(lookup.at(k), v);
        EXPECT_TRUE(lookup.contains(k));
        EXPECT_EQ(lookup[k], v);
    }

    lookup.at(k1) = 19;
    EXPECT_EQ(lookup.at(k1), 19);
    ASSERT_TRUE(lookup.contains(k2));
    lookup[k2] = -3;
    EXPECT_EQ(lookup[k2], -3);
}

TEST(util, Lookup_access) {
    using namespace tempo;
    auto keys = {3, 1, 7};
    std::vector values{9, 2, 1};
    Lookup lookup(keys, 0, values);
    testAccess(lookup, keys, values, 4, 2);
}

TEST(util, Lookup_projection) {
    using namespace tempo;
    auto literals = {
        makeBooleanLiteral<int>(true, 4), makeBooleanLiteral<int>(false, 4), makeBooleanLiteral<int>(true, 8)
    };

    auto values = {9, 2, 1};
    Lookup lookup(literals, 0, values, {}, IdProjection{});
    testAccess(lookup, literals, values, makeBooleanLiteral<int>(true, 5), makeBooleanLiteral<int>(false, 7));
}