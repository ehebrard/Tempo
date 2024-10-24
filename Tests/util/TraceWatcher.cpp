/**
* @author Tim Luchterhand
* @date 25.07.24
* @brief
*/


#include <gtest/gtest.h>
#include <algorithm>
#include <vector>

#include "util/TraceWatcher.hpp"

class TestSolver : std::vector<int> {
public:
    using std::vector<int>::vector;

    auto getTruthFunction() {
        return [this](tempo::var_t var) {
            if (this->at(var) == 0) {
                return tempo::TraceWatcher::TruthVal::Undefined;
            }

            return static_cast<tempo::TraceWatcher::TruthVal>(this->at(var) == 1);
        };
    }
};

struct TestWatcher : public tempo::TraceWatcher {
    using tempo::TraceWatcher::isAligned;
    using tempo::TraceWatcher::isEqual;
};


TEST(util, TraceWatcher_ctor) {
    using namespace tempo;
    TraceWatcher traceWatcher(std::views::iota(0, 5));
    EXPECT_EQ(traceWatcher.getOffset(), 0);
    EXPECT_FALSE(traceWatcher.isOnTrack());
}

TEST(util, TraceWatcher_ctor_empty) {
    using namespace tempo;
    EXPECT_THROW(TraceWatcher traceWatcher(std::vector<var_t>{}), std::runtime_error);
}

TEST(util, TraceWatcher_ctor_non_continuous) {
    using namespace tempo;
    EXPECT_THROW(TraceWatcher(std::vector{1, 0, 4, 3}), std::runtime_error);
}

TEST(util, TraceWatcher_ctor_offset) {
    using namespace tempo;
    TraceWatcher traceWatcher(std::views::iota(4, 9));
    EXPECT_FALSE(traceWatcher.isOnTrack());
    EXPECT_EQ(traceWatcher.getOffset(), 4);
}

TEST(util, TraceWatcher_register_soution) {
    using namespace tempo;
    TraceWatcher traceWatcher(std::views::iota(0, 5));
    TestSolver solver{1, 1, 1, -1, -1};
    traceWatcher.registerSolution(solver.getTruthFunction());
    EXPECT_TRUE(traceWatcher.isOnTrack());
    EXPECT_EQ(traceWatcher.getLastSolution(), (std::vector{true, true, true, false, false}));
}

TEST(util, TraceWatcher_register_solution_offset) {
    using namespace tempo;
    TraceWatcher traceWatcher(std::views::iota(2, 5));
    TestSolver solver{0, 0, 1, -1, -1};
    traceWatcher.registerSolution(solver.getTruthFunction());
    EXPECT_TRUE(traceWatcher.isOnTrack());
    EXPECT_EQ(traceWatcher.getLastSolution(), (std::vector{true, false, false}));
}

TEST(util, TraceWatcher_step) {
    using namespace tempo;
    TraceWatcher traceWatcher(std::views::iota(0, 5));
    TestSolver solver{1, 1, 1, -1, -1};
    traceWatcher.registerSolution(solver.getTruthFunction());
    traceWatcher.step(makeBooleanLiteral<int>(true, 2, 0));
    EXPECT_TRUE(traceWatcher.isOnTrack());
    traceWatcher.step(makeBooleanLiteral<int>(false, 3, 0));
    EXPECT_TRUE(traceWatcher.isOnTrack());
    traceWatcher.step(makeBooleanLiteral<int>(false, 1, 0));
    EXPECT_FALSE(traceWatcher.isOnTrack());
    traceWatcher.step(makeBooleanLiteral<int>(false, 4, 0));
    EXPECT_FALSE(traceWatcher.isOnTrack());
}

TEST(util, TraceWatcher_step_offset) {
    using namespace tempo;
    TraceWatcher traceWatcher(std::views::iota(1, 5));
    TestSolver solver{0, 1, 1, -1, -1};
    traceWatcher.registerSolution(solver.getTruthFunction());
    traceWatcher.step(makeBooleanLiteral<int>(true, 2, 0));
    EXPECT_TRUE(traceWatcher.isOnTrack());
    traceWatcher.step(makeBooleanLiteral<int>(false, 3, 0));
    EXPECT_TRUE(traceWatcher.isOnTrack());
    traceWatcher.step(makeBooleanLiteral<int>(false, 1, 0));
    EXPECT_FALSE(traceWatcher.isOnTrack());
    traceWatcher.step(makeBooleanLiteral<int>(false, 4, 0));
    EXPECT_FALSE(traceWatcher.isOnTrack());
}

TEST(util, TraceWatcher_update_on_track) {
    using namespace tempo;
    TraceWatcher traceWatcher(std::views::iota(0, 5));
    TestSolver solver{1, 1, 1, -1, -1};
    traceWatcher.registerSolution(solver.getTruthFunction());
    solver = {0, 1, 0, 0, -1};
    auto conflicting = traceWatcher.updateOnTrack(solver.getTruthFunction());
    EXPECT_TRUE(conflicting.empty());
    EXPECT_TRUE(traceWatcher.isOnTrack());
    solver = {1, -1, 0, 0, 1};
    conflicting = traceWatcher.updateOnTrack(solver.getTruthFunction());
    ASSERT_EQ(conflicting.size(), 2);
    std::ranges::sort(conflicting);
    EXPECT_EQ(conflicting.front(), (std::pair{1u, true}));
    EXPECT_EQ(conflicting.back(), (std::pair{4u, false}));
    EXPECT_FALSE(traceWatcher.isOnTrack());
}

TEST(util, TraceWatcher_update_on_track_offset) {
    using namespace tempo;
    TraceWatcher traceWatcher(std::views::iota(2, 7));
    TestSolver solver{0, 0, 1, 1, 1, -1, -1};
    traceWatcher.registerSolution(solver.getTruthFunction());
    solver = {0, 0, 0, 1, 0, 0, -1};
    auto conflicting = traceWatcher.updateOnTrack(solver.getTruthFunction());
    EXPECT_TRUE(conflicting.empty());
    EXPECT_TRUE(traceWatcher.isOnTrack());
    solver = {0, 0, 1, -1, 0, 0, 1};
    conflicting = traceWatcher.updateOnTrack(solver.getTruthFunction());
    ASSERT_EQ(conflicting.size(), 2);
    std::ranges::sort(conflicting);
    EXPECT_EQ(conflicting.front(), (std::pair{3u, true}));
    EXPECT_EQ(conflicting.back(), (std::pair{6u, false}));
    EXPECT_FALSE(traceWatcher.isOnTrack());
}

TEST(util, TraceWatcher_set_on_track) {
    using namespace tempo;
    TraceWatcher traceWatcher(std::views::iota(0, 5));
    traceWatcher.setOnTrack(true);
    EXPECT_TRUE(traceWatcher.isOnTrack());
    traceWatcher.setOnTrack(false);
    EXPECT_FALSE(traceWatcher.isOnTrack());
}

TEST(util, TraceWatcher_getVariablesOnTrack_step) {
    using namespace tempo;
    TraceWatcher traceWatcher(std::views::iota(0, 5));
    EXPECT_TRUE(traceWatcher.getVariablesOnTrack().empty());
    TestSolver solver{1, 1, 1, -1, -1};
    traceWatcher.registerSolution(solver.getTruthFunction());
    EXPECT_TRUE(traceWatcher.getVariablesOnTrack().empty());
    traceWatcher.step(makeBooleanLiteral<int>(true, 1));
    ASSERT_EQ(traceWatcher.getVariablesOnTrack().size(), 1);
    EXPECT_EQ(traceWatcher.getVariablesOnTrack().front(), std::pair(1u, true));
    traceWatcher.step(makeBooleanLiteral<int>(false, 3));
    ASSERT_EQ(traceWatcher.getVariablesOnTrack().size(), 2);
    EXPECT_EQ(traceWatcher.getVariablesOnTrack().front(), std::pair(1u, true));
    EXPECT_EQ(traceWatcher.getVariablesOnTrack().back(), std::pair(3u, false));
    traceWatcher.step(makeBooleanLiteral<int>(false, 0));
    EXPECT_EQ(traceWatcher.getVariablesOnTrack(), (serialization::Branch{{1, true}, {3, false}}));
}

TEST(util, TraceWatcher_getVariablesOnTrack_register_solution) {
    using namespace tempo;
    TraceWatcher traceWatcher(std::views::iota(0, 5));
    EXPECT_TRUE(traceWatcher.getVariablesOnTrack().empty());
    TestSolver solver{1, 1, 1, -1, -1};
    traceWatcher.registerSolution(solver.getTruthFunction());
    EXPECT_TRUE(traceWatcher.getVariablesOnTrack().empty());
    traceWatcher.step(makeBooleanLiteral<int>(true, 1));
    traceWatcher.step(makeBooleanLiteral<int>(false, 3));
    traceWatcher.registerSolution(solver.getTruthFunction());
    EXPECT_TRUE(traceWatcher.getVariablesOnTrack().empty());
}

TEST(util, TraceWatcher_getVariablesOnTrack_update_on_track) {
    using namespace tempo;
    TraceWatcher traceWatcher(std::views::iota(0, 5));
    EXPECT_TRUE(traceWatcher.getVariablesOnTrack().empty());
    TestSolver solver{1, 1, 1, -1, -1};
    traceWatcher.registerSolution(solver.getTruthFunction());
    solver = {0, 1, 0, 0, -1};
    traceWatcher.updateOnTrack(solver.getTruthFunction());
    ASSERT_TRUE(traceWatcher.isOnTrack());
    EXPECT_EQ(traceWatcher.getVariablesOnTrack(), (serialization::Branch{{1, true}, {4, false}}));
    traceWatcher.step(makeBooleanLiteral<int>(true, 3));
    solver = {1, 0, 0, -1, -1};
    traceWatcher.updateOnTrack(solver.getTruthFunction());
    EXPECT_EQ(traceWatcher.getVariablesOnTrack(), (serialization::Branch{{0, true}, {3, false}, {4, false}}));
}

TEST(util, TraceWatcher_equal) {
    using enum tempo::TraceWatcher::TruthVal;
    EXPECT_TRUE(TestWatcher::isEqual(True, true));
    EXPECT_TRUE(TestWatcher::isEqual(False, false));
    EXPECT_FALSE(TestWatcher::isEqual(True, false));
    EXPECT_FALSE(TestWatcher::isEqual(False, true));
    EXPECT_FALSE(TestWatcher::isEqual(Undefined, true));
    EXPECT_FALSE(TestWatcher::isEqual(Undefined, false));
}

TEST(util, TraceWatcher_aligned) {
    using enum tempo::TraceWatcher::TruthVal;
    EXPECT_TRUE(TestWatcher::isAligned(True, true));
    EXPECT_TRUE(TestWatcher::isAligned(False, false));
    EXPECT_FALSE(TestWatcher::isAligned(True, false));
    EXPECT_FALSE(TestWatcher::isAligned(False, true));
    EXPECT_TRUE(TestWatcher::isAligned(Undefined, true));
    EXPECT_TRUE(TestWatcher::isAligned(Undefined, false));
}
