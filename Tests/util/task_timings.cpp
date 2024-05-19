/**
 * @author Tim Luchterhand
 * @date 10.11.23.
 */

#include <gtest/gtest.h>

#include "util/task_timings.hpp"
#include "util/Matrix.hpp"
#include "Tests/testing.hpp"

TEST(task_timings, minDuration) {
    using namespace tempo;
    Matrix<int> distances(10, 10);
    tempo::testing::setTaskDurations(2, 4, 8, 2, 9, distances);
    EXPECT_EQ(minDuration(2, distances), 4);
}

TEST(task_timings, maxDuration) {
    using namespace tempo;
    Matrix<int> distances(10, 10);
    tempo::testing::setTaskDurations(2, 4, 8, 2, 9, distances);
    EXPECT_EQ(maxDuration(2, distances), 8);
}

TEST(task_timings, earliestStart) {
    using namespace tempo;
    Matrix<int> distances(10, 10);
    tempo::testing::setTaskDurations(2, 4, 8, 2, 9, distances);
    EXPECT_EQ(earliestStartTime(2, distances), 2);
}
TEST(task_timings, latestCompletion) {
    using namespace tempo;
    Matrix<int> distances(10, 10);
    tempo::testing::setTaskDurations(2, 4, 8, 2, 9, distances);
    EXPECT_EQ(latestCompletion(2, distances), 9);
}

TEST(task_timings, taskDistance) {
    using namespace tempo;
    Matrix<int> distances(10, 10);
    tempo::testing::setTaskDistance(1, 3, -4, distances);
    EXPECT_EQ(taskDistance(1, 3, distances), -4);
}

TEST(task_timings, upperBound) {
    using namespace tempo;
    Matrix<int> distances(10, 10);
    tempo::testing::setUpperBound(8, distances);
    EXPECT_EQ(upperBound(distances), 8);
}
TEST(task_timings, lowerBound) {
    using namespace tempo;
    Matrix<int> distances(10, 10);
    tempo::testing::setLowerBound(8, distances);
    EXPECT_EQ(lowerBound(distances), 8);
}
