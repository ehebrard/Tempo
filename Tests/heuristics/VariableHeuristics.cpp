/**
* @author Tim Luchterhand
* @date 24.08.24
* @brief
*/

#include <gtest/gtest.h>
#include <vector>

#include "Solver.hpp"
#include "src/header/heuristics/RankingHeuristic.hpp"
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
/*
TEST(variable_heuristics, ranking_base) {
    using namespace tempo;
    TestRanker ranker;
    std::vector<var_t> vars{5, 2, 7, 4, 6, 3};
    Solver s;
    EXPECT_EQ(ranker.bestVariable(vars, s), 7);
}*/

TEST(variable_heuristics, tightest) {
    using namespace tempo;
    using namespace heuristics;
    Tightest tightest;

}