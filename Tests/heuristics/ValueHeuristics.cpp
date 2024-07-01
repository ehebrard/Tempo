/**
 * @author Tim Luchterhand
 * @date 13.05.24
 * @brief
 */

#include <gtest/gtest.h>
#include <vector>

#include "heuristics/SolutionGuided.hpp"
#include "heuristics/TightestValue.hpp"
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

struct LitProvider {
    struct Storage {
        explicit Storage(tempo::Literal<int> lit) : lit(lit) {}

        tempo::Literal<int> lit;

        auto getLiteral(bool sign, tempo::var_t) const { return sign ? lit : ~lit; }

        auto getEdge(bool sign, tempo::var_t x) const -> tempo::DistanceConstraint<int> {
            return sign ? tempo::DistanceConstraint{x, lit.semantic(), 0} :
                   tempo::DistanceConstraint{lit.semantic(), x, 0};
        }

        bool hasSemantic(tempo::var_t) const {
            return lit.hasSemantic();
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
    using
    enum tempo::VariableType;
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
/*
TEST(value_heuristics, SolutionGuided) {
    using namespace tempo;
    using namespace tempo::heuristics;
    EXPECT_TRUE((value_heuristic<SolutionGuided, Scheduler<int>>));
    DummyScheduler sched;
    SolutionGuided h(0);
    EXPECT_EQ(h.choose(1, sched), TightestValue::choose(1, sched));
    sched.invert = true;
    EXPECT_EQ(h.choose(1, sched), TightestValue::choose(1, sched));
    sched.solution = {true, false, false};
    sched.hasSol = true;
    EXPECT_EQ(h.choose(0, sched), POS(0));
    EXPECT_EQ(h.choose(1, sched), NEG(1));
    EXPECT_EQ(h.choose(2, sched), NEG(2));
}*/