/**
* @author Tim Luchterhand
* @date 13.09.24
* @brief
*/

#include <gtest/gtest.h>
#include <vector>

#include "testing.hpp"
#include "util/edge_distance.hpp"


TEST(util, distance_bound_estimation) {
    using namespace tempo;
    tempo::testing::heuristics::LitProvider provider(makeBooleanLiteral<int>(true, 0, 3), {9, 7, 7, 4}, {8, 7, 3, 3});
    EXPECT_EQ(boundEstimation(DistanceConstraint{2, 3, 0}, provider), boundEstimation(true, 2, provider));
    EXPECT_EQ(boundEstimation(DistanceConstraint{2, 0, 0}, provider), -1);
    EXPECT_EQ(boundEstimation(DistanceConstraint{0, 1, 0}, provider), 2);
}

TEST(util, distance_bound_estimation_null_edge) {
    using namespace tempo;
    using namespace tempo::testing;
    DummyScheduler scheduler(std::vector<int>{9, 7, 7, 4}, std::vector<int>{8, 7, 3, 3});
    EXPECT_FALSE(boundEstimation(DistanceConstraint{Constant::NoVar, 2, 0}, scheduler).has_value());
    EXPECT_FALSE(boundEstimation(DistanceConstraint{Constant::NoVar, 2, 0.3}, scheduler).has_value());
    EXPECT_FALSE(boundEstimation(DistanceConstraint{Constant::NoVar, 2, 0.3f}, scheduler).has_value());
    tempo::testing::heuristics::LitProvider provider(makeBooleanLiteral<int>(true, -1));
    EXPECT_FALSE(boundEstimation(true, 0, provider).has_value());
    EXPECT_FALSE(boundEstimation(false, 0, provider).has_value());
}
