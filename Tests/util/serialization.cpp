/**
* @author Tim Luchterhand
* @date 02.05.24
* @brief
*/

#include <gtest/gtest.h>
#include <Iterators.hpp>

#include "util/Matrix.hpp"
#include "util/parsing/format.hpp"

TEST(serialization, Matrix_rt) {
    using namespace tempo;
    auto data = {1, 2, 3, 4, 5, 6};
    Matrix<int> matrix(2, 3, data);
    nlohmann::json j = matrix;
    auto recovered = j.get<Matrix<int>>();
    EXPECT_EQ(recovered.numRows(), matrix.numRows());
    EXPECT_EQ(recovered.numColumns(), matrix.numColumns());
    EXPECT_EQ(recovered.storageLayout(), matrix.storageLayout());
    EXPECT_EQ(recovered.rawData(), matrix.rawData());
}

TEST(serialization, ProblemInstance_rt) {
    ProblemInstance instance{.lowerBound = 3, .optimalSolution = 11, .durations = {3, 1}, .constraints = {{0, 1, 5}},
                             .resources = {{{1}, {6}, {{1, 2, 3}, {4, 5, 6}}, 14}, {{0, 1}, {2, 3}, {{4, 2, 1}}, 3}}};
    nlohmann::json j = instance;
    auto recovered = j.get<ProblemInstance>();
    EXPECT_EQ(recovered.lowerBound, instance.lowerBound);
    EXPECT_EQ(recovered.optimalSolution, instance.optimalSolution);
    EXPECT_EQ(recovered.durations, instance.durations);
    EXPECT_EQ(recovered.constraints, instance.constraints);
    ASSERT_EQ(recovered.resources.size(), recovered.resources.size());
    for (auto [rec, orig] : iterators::const_zip(recovered.resources, instance.resources)) {
        EXPECT_EQ(rec.capacity, orig.capacity);
        EXPECT_EQ(rec.demand, orig.demand);
        EXPECT_EQ(rec.transition, orig.transition);
    }
}

