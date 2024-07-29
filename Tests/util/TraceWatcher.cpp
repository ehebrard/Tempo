/**
* @author Tim Luchterhand
* @date 25.07.24
* @brief
*/


#include <gtest/gtest.h>
#include <algorithm>

#include "util/TraceWatcher.hpp"

class TestSolver : std::vector<int> {
public:
    using std::vector<int>::vector;

    auto getTruthFunction() {
        return [this](tempo::var_t var) {
            if (this->at(var) == 0) {
                ADD_FAILURE() << "variable " << var << " is undefined";
            }

            return this->at(var) == 1;
        };
    }

    auto getAlignmentFunction() {
        return [this](tempo::var_t var, bool val) {
            return this->at(var) == 0 or (this->at(var) == 1) == val;
        };
    }
};


TEST(util, TraceWatcher_ctor) {
    using namespace tempo;
    TraceWatcher traceWatcher(10);
    EXPECT_FALSE(traceWatcher.isOnTrack());
}

TEST(util, TraceWatcher_register_soution) {
    using namespace tempo;
    TraceWatcher traceWatcher(5);
    TestSolver solver{1, 1, 1, -1, -1};
    traceWatcher.registerSolution(solver.getTruthFunction());
    EXPECT_TRUE(traceWatcher.isOnTrack());
    EXPECT_EQ(traceWatcher.getLastSolution(), (std::vector{true, true, true, false, false}));
}

TEST(util, TraceWatcher_step) {
    using namespace tempo;
    TraceWatcher traceWatcher(5);
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

TEST(util, TraceWatcher_update_on_track) {
    using namespace tempo;
    TraceWatcher traceWatcher(5);
    TestSolver solver{1, 1, 1, -1, -1};
    traceWatcher.registerSolution(solver.getTruthFunction());
    solver = {0, 1, 0, 0, -1};
    auto conflicting = traceWatcher.updateOnTrack(solver.getAlignmentFunction());
    EXPECT_TRUE(conflicting.empty());
    EXPECT_TRUE(traceWatcher.isOnTrack());
    solver = {1, -1, 0, 0, 1};
    conflicting = traceWatcher.updateOnTrack(solver.getAlignmentFunction());
    ASSERT_EQ(conflicting.size(), 2);
    std::ranges::sort(conflicting);
    EXPECT_EQ(conflicting.front(), (std::pair{1u, true}));
    EXPECT_EQ(conflicting.back(), (std::pair{4u, false}));
    EXPECT_FALSE(traceWatcher.isOnTrack());
}

TEST(util, TraceWatcher_set_on_track) {
    using namespace tempo;
    TraceWatcher traceWatcher(5);
    traceWatcher.setOnTrack(true);
    EXPECT_TRUE(traceWatcher.isOnTrack());
    traceWatcher.setOnTrack(false);
    EXPECT_FALSE(traceWatcher.isOnTrack());
}
