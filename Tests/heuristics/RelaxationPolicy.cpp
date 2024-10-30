/**
* @author Tim Luchterhand
* @date 24.10.24
* @brief
*/

#include <gtest/gtest.h>
#include <array>
#include <ranges>

#include "testing.hpp"
#include "heuristics/RelaxationPolicy.hpp"

template<typename ...Args>
auto Lit(Args ...args) {
    return tempo::makeBooleanLiteral<int>(args...);
}

auto makeTestTaskMap() {
    using namespace tempo;
    using namespace std::views;
    namespace h = heuristics;
    namespace t = tempo::testing;
    using std::ranges::subrange;
    using M = Matrix<Literal<int>>;
    constexpr auto C = Contradiction<int>;
    auto tasks = t::createDummyTasks(5);
    tasks[0]._id_ = 4; tasks[1]._id_ = 2; tasks[2]._id_ = 6;
    tasks[3]._id_ = 3; tasks[4]._id_ = 1;
    std::array resources{
            t::DummyResourceExpression(std::array{tasks[0], tasks[1]},
                                       M(2, 2, {C, Lit(true, 4), Lit(false, 4), C})),
            t::DummyResourceExpression(std::array{tasks[1], tasks[3]},
                                       M(2, 2, {C, Lit(true, 4), Lit(false, 9), C})),
            t::DummyResourceExpression(std::array{tasks[2], tasks[1], tasks[4]},
                                       M(3, 3, {C, Lit(true, 17), Lit(true, 8),
                                                Lit(false, 9), C, Lit(true, 18),
                                                Lit(true, 21), Lit(false, 16), C}))
    };

    return std::make_pair(
            h::detail::TaskVarMap<int>(tasks | transform([](const auto &t) { return Interval<int>(t); }), resources),
            std::move(tasks));
}

TEST(relaxation_policies, TaskVarMap_mapping) {
    using namespace tempo;
    using B = BooleanVar<int>;
    auto [map, tasks] = makeTestTaskMap();
    EXPECT_EQ(map(tasks[0]), std::vector{B(4)});
    EXPECT_EQ(map(tasks[1]), (std::vector{B(4), B(9), B(18)}));
    EXPECT_EQ(map(tasks[2]), (std::vector{B(8), B(17)}));
    EXPECT_EQ(map(tasks[3]), std::vector{B(9)});
    EXPECT_EQ(map(tasks[4]), (std::vector{B(16), B(21)}));
}

TEST(relaxation_policies, TaskVarMap_getTaskLiterals) {
    using namespace tempo;
    using B = BooleanVar<int>;
    auto [map, tasks] = makeTestTaskMap();
    auto variables = map.getTaskLiterals(std::array{tasks[2], tasks[3]});
    EXPECT_EQ(variables, (std::vector{B(8), B(9), B(17)}));
    variables = map.getTaskLiterals(std::array{tasks[0], tasks[3], tasks[4]});
    EXPECT_EQ(variables, (std::vector{B(4), B(9), B(16), B(21)}));
    variables = map.getTaskLiterals(std::array{tasks[0], tasks[1], tasks[3]});
    EXPECT_EQ(variables, (std::vector{B(4), B(9), B(18)}));
}