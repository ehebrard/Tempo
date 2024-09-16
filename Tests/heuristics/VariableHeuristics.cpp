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
    Tightest tightest;
    LitProvider provider(makeBooleanLiteral<int>(true, 0, 2), {10, 8, 11, 6}, {5, 8, 9, 2});
    EXPECT_EQ(tightest.getCost(0, provider), 6);
    provider.boolean.lit = makeBooleanLiteral<int>(true, 3, 2);
    EXPECT_EQ(tightest.getCost(0, provider), 9);
    provider.boolean.lit = makeBooleanLiteral<int>(true, 0, Constant::NoSemantic);
    EXPECT_EQ(tightest.getCost(0, provider), 1);
}