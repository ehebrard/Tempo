/**
* @author Tim Luchterhand
* @date 11.09.24
* @brief
*/

#include <gtest/gtest.h>

#include "nn/heat_map_utils.hpp"
#include "Tests/testing.hpp"
#include "util/Matrix.hpp"


TEST(nn_utils, dstToBayes) {
    using namespace tempo::nn;
    using tempo::testing::random_float;
    for (int i = 0; i < 100; ++i) {
        auto a = random_float(0.0, 100.0);
        auto p = random_float(0.0, 1.0);
        EXPECT_NEAR(dstToBayes(a, a), 0.5, 0.00001);
        EXPECT_NEAR(dstToBayes(p, 1 - p), p, 0.00001);
    }
}

TEST(nn_utils, prob_mass_extraction) {
    using namespace tempo;
    Matrix<nn::DataType> heatMap(3, 3, nn::GNN::NoValue);
    heatMap(0, 1) = 0.9;
    heatMap(1, 0) = 0.1;
    heatMap(2, 1) = 1;
    heatMap(1, 2) = 1;

    EXPECT_EQ(nn::probabilityMass(0, 1, heatMap), nn::DataType(0.9));
    EXPECT_EQ(nn::probabilityMass(2, 1, heatMap), 1);
    EXPECT_THROW(nn::probabilityMass(2, 2, heatMap), std::runtime_error);
}
