/**
* @author Tim Luchterhand
* @date 01.07.24
* @brief
*/

#include <gtest/gtest.h>

#include "util/SchedulingProblemHelper.hpp"
#include "Model.hpp"

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
    EXPECT_EQ(mapping.size(), 3);
    EXPECT_TRUE(mapping.contains(6));
    EXPECT_TRUE(mapping.contains(4));
    EXPECT_FALSE(mapping.contains(-2));
    EXPECT_FALSE(mapping.contains(8));
    EXPECT_EQ(mapping(5), 2);
    EXPECT_EQ(mapping(7), mapping(6));
    EXPECT_EQ(mapping(4), 0);
    std::vector<Interval<>> empty;
    VarTaskMapping emptyMapping(empty);
    EXPECT_EQ(emptyMapping.size(), 0);
    EXPECT_FALSE(emptyMapping.contains(0));
    EXPECT_FALSE(emptyMapping.contains(3));
}

TEST(util, VarTaskMapping_non_continuous) {
    using namespace tempo;
    std::vector<Interval<>> tasks(3);
    tasks[0].start = {4, 0};
    tasks[0].end = {5, 0};
    tasks[1].start = {6, 0};
    tasks[1].end = {7, 0};
    tasks[2].start = {9, 0};
    tasks[2].end = {10, 2};
    EXPECT_THROW(VarTaskMapping{tasks}, std::runtime_error);
}


class BoundProvider {
    std::vector<int> u, l;
public:
    BoundProvider(std::vector<int> upper, std::vector<int> lower) : u(std::move(upper)), l(std::move(lower)) {}
    [[nodiscard]] int upper(tempo::var_t var) const { return u.at(var); }
    [[nodiscard]] int lower(tempo::var_t var) const { return l.at(var); }
};


TEST(util, SchedulingView_basic) {
    using namespace tempo;
    std::vector<Interval<int>> tasks{{{0, 0}, {1, 0}, 3, 6}, {{2, 0}, {2, 5}, 5, 5}, {{3, 0}, {4, 0}, 1, 7}};
    SchedulingProblemHelper<int, DisjunctiveResource<>> schedulingProb(tasks, {}, {});
    EXPECT_EQ(tasks, schedulingProb.tasks());
    EXPECT_TRUE(schedulingProb.hasTask(2));
    EXPECT_FALSE(schedulingProb.hasTask(3));

    EXPECT_TRUE(schedulingProb.hasVariable(0));
    EXPECT_TRUE(schedulingProb.hasVariable(1));
    EXPECT_TRUE(schedulingProb.hasVariable(2));
    EXPECT_TRUE(schedulingProb.hasVariable(3));
    EXPECT_TRUE(schedulingProb.hasVariable(4));
    EXPECT_FALSE(schedulingProb.hasVariable(5));

    EXPECT_EQ(schedulingProb.getTask(3), tasks.back());
    EXPECT_EQ(schedulingProb.getTask(1), tasks.front());
}

TEST(util, SchedulingView_task_distances) {
    using namespace tempo;
    std::vector<Interval<int>> tasks{{{0, 0}, {1, 0}, 3, 6}, {{2, 0}, {2, 5}, 5, 5}, {{3, 0}, {4, 0}, 1, 7}};
    SchedulingProblemHelper<int, DisjunctiveResource<>> schedulingProb(tasks, {}, {});
    BacktrackEnvironment env;
    DirectedGraph<LabeledEdge<int>> graph(5, &env);
    graph.emplace_edge(2, 4, -2);
    graph.emplace_edge(0, 2, -1);
    graph.emplace_edge(0, 3, 4);
    BoundProvider bounds({6, 5, 3, 8, 4}, {3, 2, 0, 4, 1});
    auto distance = schedulingProb.getTaskDistances(graph, bounds);
    EXPECT_EQ(distance(0, 1), 4);
    EXPECT_EQ(distance(1, 0), 1);
    EXPECT_EQ(distance(0, 2), 5);
    EXPECT_EQ(distance(2, 0), 6);
    EXPECT_EQ(distance(2, 1), 13);
    EXPECT_EQ(distance(1, 2), -2);
}