/**
 * @author Tim Luchterhand
 * @data 26.07.23.
 */

#include <gtest/gtest.h>
#include <torch/nn/module.h>
#include <Iterators.hpp>

#include "nn/GNN.hpp"
#include "nn/tensor_utils.hpp"
#include "testing.hpp"

class TestGNN: public tempo::nn::GNN {
public:
    explicit TestGNN(const std::filesystem::path &path) : tempo::nn::GNN(path) {}
    auto getModel() const noexcept -> const torch::jit::Module& {
        return this->model;
    }
};

class TestEdgeRegressor : public tempo::nn::EdgeRegressor {
public:
    using tempo::nn::EdgeRegressor::extractHeatMap;
};

TEST(nn_GNN, base_gnn_ctor) {
    TestGNN gnn(tempo::testing::TestData::TestNN);
    EXPECT_FALSE(gnn.getModel().is_training());
}

TEST(nn_GNN, edge_regressor_heat_map){
    using namespace tempo::nn;
    EdgeVector edges{{1, 2}, {2, 4}, {4, 2}, {3, 1}, {1, 3}};
    const auto edgeTensor = util::makeIndexTensor(edges);
    std::vector<DataType> probsData{0.3, 0.8, 0.2, 0.4, 0.6};
    auto probs = torch::from_blob(probsData.data(), static_cast<long>(probsData.size()), torch::dtype<DataType>());
    auto heatMap = TestEdgeRegressor::extractHeatMap(probs, edgeTensor, 5);
    ASSERT_EQ(heatMap.numRows(), heatMap.numColumns());
    ASSERT_EQ(heatMap.numRows(), 5);
    for (auto [edge, prob] : iterators::const_zip(edges, probsData)) {
        EXPECT_EQ(heatMap.at(edge.first, edge.second), prob);
    }

    EXPECT_EQ(heatMap.at(0, 1), TestEdgeRegressor::NoValue);
    EXPECT_EQ(heatMap.at(1, 0), TestEdgeRegressor::NoValue);
    EXPECT_EQ(heatMap.at(4, 1), TestEdgeRegressor::NoValue);
}

TEST(nn_GNN, gnn_heat_map_choose_polarity) {
    using namespace tempo;
    Matrix<nn::DataType> heatMap(3, 3, nn::GNN::NoValue);
    heatMap(0, 1) = 0.9;
    heatMap(1, 0) = 0.1;
    heatMap(2, 1) = 1;
    heatMap(1, 2) = 1;
    EXPECT_NEAR(nn::EdgeRegressor::dstEdgeProbability(0, 1, heatMap), 0.9, 0.00001);
    EXPECT_NEAR(nn::EdgeRegressor::dstEdgeProbability(1, 0, heatMap), 0.1, 0.00001);
    EXPECT_DOUBLE_EQ(nn::EdgeRegressor::dstEdgeProbability(2, 1, heatMap), 0.5);
    EXPECT_DOUBLE_EQ(nn::EdgeRegressor::dstEdgeProbability(1, 2, heatMap), 0.5);
    EXPECT_THROW(nn::EdgeRegressor::dstEdgeProbability(0, 2, heatMap), std::runtime_error);
    EXPECT_THROW(nn::EdgeRegressor::dstEdgeProbability(2, 0, heatMap), std::runtime_error);
}
