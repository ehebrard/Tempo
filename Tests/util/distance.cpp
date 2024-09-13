/**
* @author Tim Luchterhand
* @date 13.09.24
* @brief
*/

#include <gtest/gtest.h>

#include "testing.hpp"
#include "util/distance.hpp"


TEST(util, distance_bound_estimation) {
    using namespace tempo;
    using namespace tempo::testing;
    heuristics::LitProvider provider(makeBooleanLiteral<int>(true, 0, 3), {9, 7, 7, 4}, {8, 7, 3, 3});
    EXPECT_EQ(boundEstimation(DistanceConstraint{2, 3, 0}, provider), boundEstimation(true, 2, provider));
    EXPECT_EQ(boundEstimation(DistanceConstraint{2, 0, 0}, provider), -1);
    EXPECT_EQ(boundEstimation(DistanceConstraint{0, 1, 0}, provider), 2);
}

TEST(util, distance_bound_estimation_null_edge) {
    using namespace tempo;
    using namespace tempo::testing;
    DummyScheduler scheduler(std::vector{9, 7, 7, 4}, std::vector{8, 7, 3, 3});
    EXPECT_TRUE(std::isnan(boundEstimation(DistanceConstraint{Constant::NoVar, 2, 0}, scheduler)));
    EXPECT_TRUE(std::isnan(boundEstimation(DistanceConstraint{Constant::NoVar, 2, 0.3}, scheduler)));
    EXPECT_TRUE(std::isnan(boundEstimation(DistanceConstraint{Constant::NoVar, 2, 0.3f}, scheduler)));
    heuristics::LitProvider provider(makeBooleanLiteral<int>(true, -1));
    EXPECT_TRUE(std::isnan(boundEstimation(true, 0, provider)));
    EXPECT_TRUE(std::isnan(boundEstimation(false, 0, provider)));
}
