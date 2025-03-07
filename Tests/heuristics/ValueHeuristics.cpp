/**
 * @author Tim Luchterhand
 * @date 13.05.24
 * @brief
 */

#include <gtest/gtest.h>
#include <vector>

#include "heuristics/SolutionGuided.hpp"
#include "heuristics/TightestValue.hpp"
#include "heuristics/PerfectValueOracle.hpp"
#include "util/serialization.hpp"
#include "Solver.hpp"
#include "testing.hpp"

struct TestValueHeuristic : tempo::heuristics::BaseBooleanHeuristic<TestValueHeuristic, int> {
    bool called = false;

    explicit TestValueHeuristic(double epsilon): BaseBooleanHeuristic(epsilon) {}

    template<typename Sched>
    auto choose(tempo::var_t, const Sched &) {
        called = true;
        return tempo::makeBooleanLiteral<int>(true, 0, 0);
    }

    using BaseBooleanHeuristic::valueDecisionImpl;
};

struct TestBaseHeuristic {
    bool &called;
    tempo::Literal<int> response;

    TestBaseHeuristic(bool &called, tempo::Literal<int> response) : called(called), response(response) {}

    template<typename S>
    auto valueDecision(const tempo::heuristics::VariableSelection &, const S &) {
        called = true;
        return response;
    }
};

TEST(value_heuristics, base_value_heuristic) {
    using namespace tempo;
    using tempo::testing::heuristics::LitProvider;
    using enum tempo::heuristics::VariableType;
    TestValueHeuristic h(0);
    auto lit = makeBooleanLiteral<int>(true, 0, 0);
    LitProvider provider(lit);
    EXPECT_EQ(h.valueDecisionImpl({0, Boolean}, provider), lit);
    EXPECT_TRUE(h.called);
    h = TestValueHeuristic(1);
    h.valueDecisionImpl({0, Boolean}, provider);
    EXPECT_FALSE(h.called);
    EXPECT_THROW(TestValueHeuristic(2), std::runtime_error);
}

TEST(value_heuristics, TightestValue) {
    using namespace tempo;
    using namespace tempo::heuristics;
    using tempo::testing::heuristics::LitProvider;
    EXPECT_TRUE((value_heuristic<TightestValue<int>, Solver<int>>));
    auto lit = makeBooleanLiteral<int>(true, 0, 1);
    LitProvider provider(lit, {7, 5}, {3, 4});
    EXPECT_EQ(TightestValue<int>::choose(0, provider), lit);
    provider = LitProvider(lit, {7, 8}, {3, 4});
    EXPECT_EQ(TightestValue<int>::choose(0, provider), ~lit);
}

//TODO fix test
/*
TEST(value_heuristics, SolutionGuided) {
    using namespace tempo;
    using namespace tempo::heuristics;
    using tempo::testing::heuristics::LitProvider;
    bool called = false;
    const auto baseResponse = makeBooleanLiteral<int>(true, 24, 19);
    SolutionGuided<TestBaseHeuristic> h(0, called, baseResponse);
    LitProvider solver(Literal<int>(3, 17, tempo::detail::Boolean{}));
    auto res = h.valueDecision({3, VariableType::Boolean}, solver);
    EXPECT_TRUE(called);
    EXPECT_EQ(baseResponse, res);
    EXPECT_EQ(res.semantic(), baseResponse.semantic());
    called = false;
    solver.boolean.solution = {true, false, false, true};
    res = h.valueDecision({0, VariableType::Boolean}, solver);
    EXPECT_FALSE(called);
    EXPECT_EQ(res.id(), 3);
    solver.boolean.lit = {1, 15, tempo::detail::Boolean{}};
    res = h.valueDecision({0, VariableType::Boolean}, solver);
    EXPECT_FALSE(called);
    EXPECT_EQ(res.id(), 0);
}*/

TEST(value_heuristics, oracle) {
    using namespace tempo;
    using namespace tempo::heuristics;
    using tempo::testing::heuristics::LitProvider;
    const serialization::Solution solution(0, 0, {{0, true}, {1, true}, {2, false}, {3, true}, {4, false}});
    PerfectValueHeuristic oracle(0, solution);
    const auto lit = makeBooleanLiteral<int>(true, 0, 0);
    LitProvider provider(lit);
    for (auto [var, val] : solution.decisions) {
        EXPECT_EQ(oracle.choose(var, provider), (val ? lit : ~lit));
    }
}

TEST(value_heuristics, oracle_determinism) {
    using namespace tempo;
    using namespace tempo::heuristics;
    using tempo::testing::heuristics::LitProvider;
    const serialization::Solution solution(0, 0, {{0, true}, {1, true}, {2, false}, {3, true}, {4, false}});
    const auto lit = makeBooleanLiteral<int>(true, 0, 0);
    LitProvider provider(lit);
    std::vector<bool> deviations(solution.decisions.size());
    for (int c = 0; c < 100; c++) {
        PerfectValueHeuristic oracle(0.5, solution);
        for (auto [decision, deviation] : iterators::zip(solution.decisions, deviations)) {
            deviation = (decision.second ? lit : ~lit) == oracle.choose(decision.first, provider);
        }

        for (int i = 0; i < 100; ++i) {
            for (auto [decision, deviation] : iterators::const_zip(solution.decisions, deviations)) {
                EXPECT_EQ(deviation, (decision.second ? lit : ~lit) == oracle.choose(decision.first, provider));
            }
        }
    }
}

TEST(value_heuristics, oracle_polarity_oob) {
    using namespace tempo;
    using namespace tempo::heuristics;
    using tempo::testing::heuristics::LitProvider;
    const serialization::Solution solution(0, 0, {{0, true}});
    PerfectValueHeuristic oracle(0, solution);
    const auto lit = makeBooleanLiteral<int>(true, 0, 0);
    LitProvider provider(lit);
    EXPECT_EQ(oracle.choose(0, provider), lit);
}
