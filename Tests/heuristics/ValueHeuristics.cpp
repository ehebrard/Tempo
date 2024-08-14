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
#include "util/traits.hpp"
#include "Solver.hpp"

struct TestValueHeuristic
        : public tempo::heuristics::BaseBooleanHeuristic<TestValueHeuristic> {
    bool called = false;

    explicit TestValueHeuristic(double epsilon)
            : tempo::heuristics::BaseBooleanHeuristic<TestValueHeuristic>(epsilon) {}

    template<typename Sched>
    auto choose(tempo::var_t, const Sched &) {
        called = true;
        return tempo::makeBooleanLiteral<int>(true, 0, 0);
    }
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

struct LitProvider {
    struct Storage {
        explicit Storage(tempo::Literal<int> lit) : lit(lit) {}

        tempo::Literal<int> lit;
        std::vector<bool> solution{};

        auto getLiteral(bool sign, tempo::var_t) const { return sign ? lit : ~lit; }

        auto getEdge(bool sign, tempo::var_t x) const -> tempo::DistanceConstraint<int> {
            return sign ? tempo::DistanceConstraint{x, lit.semantic(), 0} :
                   tempo::DistanceConstraint{lit.semantic(), x, 0};
        }

        bool hasSemantic(tempo::var_t) const {
            return lit.hasSemantic();
        }

        bool hasSolution() const {
            return not solution.empty();
        }

        const auto &bestSolution() const {
            return solution;
        }
    };

    struct Numeric {
        int upper(tempo::var_t x) const {
            return static_cast<int>(x);
        }

        int lower(tempo::var_t x) const {
            return -static_cast<int>(x / 2);
        }
    };

    Storage boolean;
    Numeric numeric{};

    explicit LitProvider(tempo::Literal<int> lit) : boolean(lit) {}
};

TEST(value_heuristics, base_value_heuristic) {
    using namespace tempo;
    using enum tempo::heuristics::VariableType;
    TestValueHeuristic h(0);
    auto lit = makeBooleanLiteral<int>(true, 0, 0);
    LitProvider provider(lit);
    EXPECT_EQ(h.valueDecision({0, Boolean}, provider), lit);
    EXPECT_TRUE(h.called);
    h = TestValueHeuristic(1);
    h.valueDecision({0, Boolean}, provider);
    EXPECT_FALSE(h.called);
    EXPECT_THROW(TestValueHeuristic(2), std::runtime_error);
}

TEST(value_heuristics, TightestValue) {
    using namespace tempo;
    using namespace tempo::heuristics;
    EXPECT_TRUE((value_heuristic<TightestValue, Solver<int>>));
    auto lit = makeBooleanLiteral<int>(true, 0, 4);
    LitProvider provider(lit);
    EXPECT_EQ(TightestValue::choose(5, provider), lit);
    lit = makeBooleanLiteral<int>(true, 0, 4);
    provider = LitProvider{lit};
    EXPECT_EQ(TightestValue::choose(2, provider), ~lit);
}


TEST(value_heuristics, SolutionGuided) {
    using namespace tempo;
    using namespace tempo::heuristics;
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
}

TEST(value_heuristics, oracle) {
    using namespace tempo;
    using namespace tempo::heuristics;
    const serialization::Solution solution(0, 0, {{0, true}, {1, true}, {2, false}, {3, true}, {4, false}});
    PerfectValueHeuristic oracle(0, solution);
    const auto lit = makeBooleanLiteral<int>(true, 0, 0);
    LitProvider provider(lit);
    for (auto [var, val] : solution.decisions) {
        EXPECT_EQ(oracle.choose(var, provider), (val ? lit : ~lit));
    }
}
