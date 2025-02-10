/**
* @author Tim Luchterhand
* @date 05.02.25
* @file SolutionPool.cpp
* @brief
*/

#include <gtest/gtest.h>
#include <unordered_set>

#include "heuristics/LNS/SolutionPool.hpp"


TEST(util, Pool_ctor) {
    using namespace tempo::lns;
    Pool<int> pool;
    EXPECT_EQ(pool.size(), 0);
    EXPECT_TRUE(pool.empty());
    EXPECT_EQ(pool.capacity(), std::numeric_limits<size_t>::max());
    pool = Pool<int>(10);
    EXPECT_EQ(pool.size(), 0);
    EXPECT_TRUE(pool.empty());
    EXPECT_EQ(pool.capacity(), 10);
}

TEST(util, Pool_add_pop_unlimited) {
    using namespace tempo::lns;
    Pool<int> pool;
    for (int i = 0; i < 100; i++) {
        pool.emplace_back(i);
        EXPECT_EQ(pool.size(), i + 1);
        EXPECT_EQ(pool.peekLast(), i);
        EXPECT_FALSE(pool.empty());
    }

    int last = pool.peekLast();
    while (!pool.empty()) {
        EXPECT_EQ(pool.popLast(), last);
        EXPECT_EQ(pool.size(), last);
        --last;
    }

    pool.emplace_back(15);
    pool.emplace_back(17);
    pool.popLast();
    pool.emplace_back(23);
    EXPECT_EQ(pool.size(), 2);
    EXPECT_EQ(pool.peekLast(), 23);
}

TEST(util, Pool_add_pop_limited) {
    using namespace tempo::lns;
    Pool<int> pool(5);
    for (int i = 0; i < 7; i++) {
        pool.emplace_back(i);
        EXPECT_EQ(pool.size(), std::min(i + 1, 5));
        EXPECT_EQ(pool.peekLast(), i);
        EXPECT_FALSE(pool.empty());
    }

    int last = pool.peekLast();
    auto size = pool.size();
    while (!pool.empty()) {
        EXPECT_EQ(pool.popLast(), last--);
        EXPECT_EQ(pool.size(), --size);
    }

    pool = Pool<int>(2);
    pool.emplace_back(3);
    pool.emplace_back(8);
    pool.emplace_back(17);
    pool.popLast();
    pool.emplace_back(23);
    EXPECT_EQ(pool.popLast(), 23);
    EXPECT_EQ(pool.popLast(), 8);
}

TEST(util, Pool_clear) {
    using namespace tempo::lns;
    Pool<int> pool;
    for (int i = 0; i < 100; i++) {
        pool.emplace_back(i);
    }

    EXPECT_EQ(pool.size(), 100);
    pool.clear();
    EXPECT_EQ(pool.size(), 0);
    EXPECT_TRUE(pool.empty());
    pool.emplace_back(15);
    pool.emplace_back(17);
    EXPECT_EQ(pool.size(), 2);
    EXPECT_EQ(pool.popLast(), 17);
    EXPECT_EQ(pool.popLast(), 15);

    pool = Pool<int>(2);
    pool.emplace_back(8);
    pool.emplace_back(17);
    pool.emplace_back(23);
    EXPECT_EQ(pool.size(), 2);
    pool.clear();
    EXPECT_EQ(pool.size(), 0);
    EXPECT_TRUE(pool.empty());
    pool.emplace_back(15);
    pool.emplace_back(17);
    EXPECT_EQ(pool.size(), 2);
    EXPECT_EQ(pool.popLast(), 17);
    EXPECT_EQ(pool.popLast(), 15);
}

TEST(util, Pool_random_unlimited) {
    using namespace tempo::lns;
    Pool<int> pool;
    std::unordered_set<int> values;
    for (int i = 0; i < 100; i++) {
        pool.emplace_back(i);
        values.emplace(i);
    }

    while (not pool.empty()) {
        EXPECT_TRUE(values.contains(pool.peekRandom()));
        int res = pool.popRandom();
        ASSERT_TRUE(values.contains(res));
        values.erase(res);
    }

    EXPECT_EQ(pool.size(), 0);
    EXPECT_TRUE(pool.empty());
    EXPECT_TRUE(values.empty());
}

TEST(util, Pool_random_limited) {
    using namespace tempo::lns;
    Pool<int> pool(20);
    auto gt = std::ranges::iota_view(30, 50);
    std::unordered_set values(gt.begin(), gt.end());
    for (int i = 0; i < 50; i++) {
        pool.emplace_back(i);
    }

    while (not pool.empty()) {
        EXPECT_TRUE(values.contains(pool.peekRandom()));
        int res = pool.popRandom();
        ASSERT_TRUE(values.contains(res));
        values.erase(res);
    }

    EXPECT_EQ(pool.size(), 0);
    EXPECT_TRUE(pool.empty());
    EXPECT_TRUE(values.empty());
}
