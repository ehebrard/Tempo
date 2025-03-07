/**
* @author Tim Luchterhand
* @date 24.08.24
* @brief
*/

#include <gtest/gtest.h>
#include <vector>

#include "testing.hpp"
#include "heuristics/RankingHeuristic.hpp"
#include "heuristics/Tightest.hpp"

struct TestRanker : public tempo::heuristics::RankingHeuristic<TestRanker> {
    template<typename T>
    [[nodiscard]] auto chooseBest(tempo::var_t x, tempo::var_t y, T &&) const {
        using namespace tempo;
        if (x == Constant::NoVar) {
            return y;
        }

        if (y == Constant::NoVar) {
            return x;
        }

        return x > y ? x : y;
    }
};

TEST(variable_heuristics, tightest) {
    using namespace tempo;
    using namespace heuristics;
    using tempo::testing::heuristics::LitProvider;
    Tightest<int> tightest;
    LitProvider provider(makeBooleanLiteral<int>(true, 0, 2), {10, 8, 11, 6}, {5, 8, 9, 2});
    EXPECT_EQ(tightest.getCost(0, provider), 6);
    provider.boolean.lit = makeBooleanLiteral<int>(true, 3, 2);
    EXPECT_EQ(tightest.getCost(3, provider), 9);
    provider.boolean.lit = makeBooleanLiteral<int>(true, 0, Constant::NoSemantic);
    EXPECT_EQ(tightest.getCost(0, provider), 1);
}

TEST(variable_heuristics, tightest_ranking) {
    using namespace tempo;
    using namespace heuristics;
    using tempo::testing::heuristics::LitProvider;
    Tightest<int> tightest;
    std::vector<var_t> variables{0, 1, 2, 3, 4};
    LitProvider provider(makeBooleanLiteral<int>(true, 0, 2), {10, 5, 7, 6, 6}, {5, 5, 4, 3, 1});
    EXPECT_EQ(tightest.bestVariable(variables, provider), 1);
}