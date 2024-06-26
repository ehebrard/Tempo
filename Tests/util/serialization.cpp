/**
* @author Tim Luchterhand
* @date 02.05.24
* @brief
*/

#include <gtest/gtest.h>
#include <Iterators.hpp>

#include "util/Matrix.hpp"

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
