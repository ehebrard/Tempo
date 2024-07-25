/**
* @author Tim Luchterhand
* @date 01.07.24
* @brief
*/

#include <gtest/gtest.h>

#include "util/SchedulingProblemHelper.hpp"
#include "Model.hpp"
#include "testing.hpp"

TEST(util, VarTaskMapping) {
    using namespace tempo;
    std::vector<Interval<>> tasks(3);
    tasks[0].start = {4, 0};
    tasks[0].end = {4, 7};
    tasks[1].start = {6, 0};
    tasks[1].end = {7, 0};
    tasks[2].start = {5, 0};
    tasks[2].end = {5, 2};
    VarTaskMapping mapping(tasks);
    EXPECT_TRUE(mapping.contains(6));
    EXPECT_TRUE(mapping.contains(4));
    EXPECT_FALSE(mapping.contains(-2));
    EXPECT_FALSE(mapping.contains(8));
    EXPECT_EQ(mapping(5), 2);
    EXPECT_EQ(mapping(7), mapping(6));
    EXPECT_EQ(mapping(4), 0);
    std::vector<Interval<>> empty;
    VarTaskMapping emptyMapping(empty);
    EXPECT_FALSE(emptyMapping.contains(0));
    EXPECT_FALSE(emptyMapping.contains(3));
}

TEST(util, VarTaskMapping_non_continuous) {
    using namespace tempo;
    std::vector<Interval<>> tasks(3);
    tasks[0].start = {4, 0};
    tasks[0].end = {5, 0};
    tasks[1].start = {18, 0};
    tasks[1].end = {19, 0};
    tasks[2].start = {9, 0};
    tasks[2].end = {10, 2};
    EXPECT_THROW(VarTaskMapping{tasks}, std::runtime_error);
}

TEST(util, SchedulingProblemHelper_basic) {
    using namespace tempo;
    using T = tempo::testing::TaskSpec;
    auto [tasks, _] = tempo::testing::createTasks(
            {T{.minDur = 3, .maxDur = 6}, T{.minDur = 5, .maxDur = 5}, T{.minDur = 1, .maxDur = 7},
             T{.minDur = 0, .maxDur = 100}});
    auto sched = tasks.back();
    tasks.pop_back();
    SchedulingProblemHelper<int, tempo::testing::Resource> schedulingProb(tasks, {}, {}, sched);
    EXPECT_EQ(tasks, schedulingProb.tasks());
    EXPECT_TRUE(schedulingProb.hasTask(0));
    EXPECT_TRUE(schedulingProb.hasTask(1));
    EXPECT_TRUE(schedulingProb.hasTask(2));
    EXPECT_FALSE(schedulingProb.hasTask(3));

    EXPECT_TRUE(schedulingProb.hasVariable(tasks[0].start.id()));
    EXPECT_TRUE(schedulingProb.hasVariable(tasks[0].end.id()));
    EXPECT_TRUE(schedulingProb.hasVariable(tasks[1].start.id()));
    EXPECT_TRUE(schedulingProb.hasVariable(tasks[1].end.id()));
    EXPECT_TRUE(schedulingProb.hasVariable(tasks[2].start.id()));
    EXPECT_TRUE(schedulingProb.hasVariable(tasks[2].end.id()));
    EXPECT_FALSE(schedulingProb.hasVariable(tasks[1].duration.id()));
    EXPECT_FALSE(schedulingProb.hasVariable(120));

    EXPECT_EQ(schedulingProb.getTask(tasks.back().end.id()), tasks.back());
    EXPECT_EQ(schedulingProb.getTask(tasks.front().id()), tasks.front());
}

TEST(util, SchedulingProblemView_task_distances) {
    using namespace tempo;
    std::vector<Interval<int>> tasks{{{0, 0}, {1, 0}, NumericVar{}},
                                     {{2, 0}, {2, 5}, NumericVar{}},
                                     {{3, 0}, {4, 0}, NumericVar{}}};
    SchedulingProblemHelper<int, tempo::testing::Resource> schedulingProb(tasks, {}, {}, Interval<int>{});
    BacktrackEnvironment env;
    DirectedGraph<LabeledEdge<int>> graph(5, &env);
    graph.emplace_edge(2, 4, -2);
    graph.emplace_edge(0, 2, -1);
    graph.emplace_edge(0, 3, 4);
    tempo::testing::BoundProvider bounds({6, 5, 3, 8, 4}, {3, 2, 0, 4, 1});
    auto distance = schedulingProb.getTaskDistances(graph, bounds);
    EXPECT_EQ(distance(0, 1), 4);
    EXPECT_EQ(distance(1, 0), 1);
    EXPECT_EQ(distance(0, 2), 5);
    EXPECT_EQ(distance(2, 0), 6);
    EXPECT_EQ(distance(2, 1), 13);
    EXPECT_EQ(distance(1, 2), -2);
}