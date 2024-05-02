/**
 * @author Tim Luchterhand
 * @date 21.07.23.
 */

#include <gtest/gtest.h>
#include <ranges>
#include <filesystem>
#include <unordered_map>
#include "nn/tensor_utils.hpp"


using tempo::nn::indexTensorOptions;
using tempo::nn::IndexType;

class FileTest : public testing::Test {
protected:
    static constexpr auto FileName = "./__test_file__";
    static constexpr auto Directory = "./__test_directory__";
    void SetUp() override {
        std::filesystem::create_directory(Directory);
    };

    void TearDown() override {
        if (std::filesystem::exists(FileName)) {
            std::filesystem::remove(FileName);
        }

        if (std::filesystem::is_directory(Directory)) {
            std::filesystem::remove_all(Directory);
        }
    }
};

TEST(nn_utils, index_read) {
    using namespace tempo::nn::util;
    auto tensor = torch::tensor({{{1,  2,  3,  4},  {5,  6,  7,  8}},
                                 {{10, 20, 30, 40}, {50, 60, 70, 80}}}, indexTensorOptions());
    EXPECT_EQ((tensor[0][0][3].item<IndexType>()), (index<IndexType>(tensor, {0, 0, 3})));
    EXPECT_EQ((tensor[1][0][1].item<IndexType>()), (index<IndexType>(tensor, {1, 0, 1})));
    EXPECT_EQ((tensor[0][1][1].item<IndexType>()), (index<IndexType>(tensor, {0, 1, 1})));
    EXPECT_EQ((tensor[1][1][0].item<IndexType>()), (index<IndexType>(tensor, {1, 1, 0})));
    EXPECT_EQ((tensor[0][0][3].item<IndexType>()), (c_index<IndexType>(tensor, {0, 0, 3})));
    EXPECT_EQ((tensor[1][0][1].item<IndexType>()), (c_index<IndexType>(tensor, {1, 0, 1})));
    EXPECT_EQ((tensor[0][1][1].item<IndexType>()), (c_index<IndexType>(tensor, {0, 1, 1})));
    EXPECT_EQ((tensor[1][1][0].item<IndexType>()), (c_index<IndexType>(tensor, {1, 1, 0})));
}

TEST(nn_utils, index_write) {
    using namespace tempo::nn::util;
    auto tensor = torch::tensor({{{1,  2,  3,  4},  {5,  6,  7,  8}},
                                 {{10, 20, 30, 40}, {50, 60, 70, 80}}}, indexTensorOptions());
    index<IndexType>(tensor, {0, 1, 2}) = 17;
    EXPECT_EQ(tensor[0][1][2].item<IndexType>(), 17);
    index<IndexType>(tensor, {1, 1, 1}) = 12;
    EXPECT_EQ(tensor[1][1][1].item<IndexType>(), 12);
    index<IndexType>(tensor, {1, 0, 3}) = 103;
    EXPECT_EQ(tensor[1][0][3].item<IndexType>(), 103);
    index<IndexType>(tensor, {0, 0, 2}) = 145;
    EXPECT_EQ(tensor[0][0][2].item<IndexType>(), 145);
}

TEST(nn_utils, index_exceptions) {
    using namespace tempo::nn::util;
    auto tensor = torch::tensor({{{1,  2,  3,  4},  {5,  6,  7,  8}},
                                 {{10, 20, 30, 40}, {50, 60, 70, 80}}}, indexTensorOptions());
    EXPECT_THROW(index<IndexType>(tensor, {0, 0, 0, 0}), std::out_of_range);
    EXPECT_THROW(index<IndexType>(tensor, {2, 2, 4}), std::out_of_range);

    auto sparseTensor = torch::empty({2, 2}, torch::layout(torch::kSparse));
    EXPECT_THROW(index<IndexType>(sparseTensor, {0, 0}), std::runtime_error);
}

TEST(nn_utils, slice_assign) {
    using namespace tempo::nn::util;
    auto tensor = torch::tensor({{{1,  2,  3,  4},  {5,  6,  7,  8}},
                                 {{10, 20, 30, 40}, {50, 60, 70, 80}}}, indexTensorOptions());
    sliceAssign<IndexType>(tensor, {1, -1, 2}, {17, 18});
    auto gtTensor = torch::tensor({{{1,  2,  3,  4},  {5,  6,  7,  8}},
                                   {{10, 20, 17, 40}, {50, 60, 18, 80}}}, indexTensorOptions());
    EXPECT_TRUE(torch::all(gtTensor == tensor).item<bool>());

    sliceAssign<IndexType>(tensor, {0, 0, -1}, {11, 12, 13, 14});
    gtTensor = torch::tensor({{{11, 12, 13, 14}, {5,  6,  7,  8}},
                              {{10, 20, 17, 40}, {50, 60, 18, 80}}}, indexTensorOptions());
    EXPECT_TRUE(torch::all(gtTensor == tensor).item<bool>());

    sliceAssign<IndexType>(tensor, {0, 1, 3}, {102});
    gtTensor = torch::tensor({{{11, 12, 13, 14}, {5,  6,  7,  102}},
                              {{10, 20, 17, 40}, {50, 60, 18, 80}}}, indexTensorOptions());
    EXPECT_TRUE(torch::all(gtTensor == tensor).item<bool>());
}

TEST(nn_utils, slice_assign_exceptions) {
    using namespace tempo::nn::util;
    auto tensor = torch::empty({2, 3, 4}, indexTensorOptions());
    EXPECT_THROW(sliceAssign(tensor, {0, 0, 0, 0}, {1}), std::out_of_range);
    EXPECT_THROW(sliceAssign(tensor, {1, 5, -1}, {1, 2, 3, 4}), std::out_of_range);
    tensor = torch::empty({2, 2}, torch::layout(torch::kSparse));
    EXPECT_THROW(sliceAssign(tensor, {1, -1}, {1, 2}), std::runtime_error);
}

TEST(nn_utils, makeIndexTensor) {
    std::vector<IndexType> from{1, 2, 3, 4};
    std::vector<IndexType> to{4, 3, 2, 1};
    auto indexTensor = tempo::nn::util::makeIndexTensor(from, to);
    EXPECT_TRUE(indexTensor.is_contiguous());
    ASSERT_EQ(indexTensor.sizes().size(), 2);
    ASSERT_EQ(indexTensor.size(0), 2);
    ASSERT_EQ(indexTensor.size(1), 4);
    for (int i = 0; i < 4; ++i) {
        auto sum = (indexTensor[0][i] + indexTensor[1][i]).item<IndexType>();
        EXPECT_EQ(sum, 5);
    }
}

TEST(nn_utils, makeIndexTensor_EdgeVector) {
    using namespace tempo::nn;
    EdgeVector edges{{1, 0}, {4, 2}, {2, 1}};
    auto indexTensor = util::makeIndexTensor(edges);
    EXPECT_TRUE(indexTensor.is_contiguous());
    ASSERT_EQ(indexTensor.sizes().size(), 2);
    ASSERT_EQ(indexTensor.size(0), 2);
    ASSERT_EQ(indexTensor.size(1), edges.size());
    for (auto [idx, edge] : iterators::enumerate(edges)) {
        EXPECT_EQ(indexTensor[0][idx].item<IndexType>(), edge.first);
        EXPECT_EQ(indexTensor[1][idx].item<IndexType>(), edge.second);
    }
}

TEST(nn_utils, getIndexSlice) {
    using namespace tempo::nn::util;
    std::vector<IndexType> from{1, 2, 3, 4};
    std::vector<IndexType> to{4, 3, 2, 1};
    auto indexTensor = makeIndexTensor(from, to);
    auto fromSlice = getIndexSlice(indexTensor, 0);
    auto toSlice = getIndexSlice(indexTensor, 1);
    EXPECT_TRUE(std::ranges::equal(from, fromSlice));
    EXPECT_TRUE(std::ranges::equal(to, toSlice));
}

TEST(nn_utils, getEdgeView) {
    using namespace tempo::nn::util;
    tempo::nn::EdgeVector edges{{1, 2}, {2, 3}, {3, 4}};
    auto tensor = makeIndexTensor(edges);
    auto rt = getEdgeView(tensor);
    EXPECT_TRUE(std::ranges::equal(edges, rt));
}
TEST_F(FileTest, load_save_tensor) {
    using namespace tempo::nn;
    auto tensor = torch::randn({3, 4, 2});
    util::saveTensor(tensor, FileName);
    auto rt = util::loadTensor(FileName);
    EXPECT_TRUE((tensor == rt).all().item<bool>());
}

struct TensorSurrogate {
    TensorSurrogate(torch::Tensor tensor) noexcept: t(std::move(tensor)) {}
    bool operator==(const TensorSurrogate &other) const {
        return (t == other.t).all().item<bool>();
    }
private:
    torch::Tensor t;
};

TEST_F(FileTest, load_save_graph) {
    using namespace tempo::nn;
    InputGraph graph;
    std::unordered_map<std::string, TensorSurrogate> gt;
    for (auto key : GraphKeys::AllKeys) {
        auto t = torch::randn({3, 4, 2});
        graph.insert(key, t);
        gt.emplace(key, t);
    }

    util::saveGraph(graph, Directory);
    auto rt = util::loadGraph(Directory);
    std::unordered_map<std::string, TensorSurrogate> rtMap;
    for (const auto &entry : rt) {
        rtMap.emplace(entry.key(), entry.value());
    }

    EXPECT_EQ(gt, rtMap);
}

auto createSimpleGraph() -> tempo::nn::InputGraph {
    using namespace tempo::nn;
    InputGraph g;
    auto taskFeats = torch::tensor({{0.1, 0.4}, {-1.4, 0.4}, {0.5, -0.2}}, dataTensorOptions());
    auto edges = torch::tensor({{0, 1, 2, 1},
                                {1, 2, 0, 0}}, indexTensorOptions());
    auto edgeFeats = torch::tensor({{0.2, 0.2}, {0.4, -0.1}, {0.9, 0.45}, {0.9, -0.4}}, dataTensorOptions());
    auto edgePairMask = torch::tensor({0, 1, 2, 0}, indexTensorOptions());
    auto resourceFeats = torch::tensor({{0.5}, {0.1}}, dataTensorOptions());
    auto resourceCons = torch::tensor({{0.1}, {0.4}, {0.7}, {0.5}}, dataTensorOptions());
    auto resourceDeps = torch::tensor({{0, 0, 1, 2}, {0, 1, 1, 1}}, indexTensorOptions());
    auto edgeResRel = torch::tensor({{0, 1, 2, 3}, {1, 1, 1, 1}}, indexTensorOptions());
    g.insert(GraphKeys::TaskFeatures, std::move(taskFeats));
    g.insert(GraphKeys::EdgeFeatures, std::move(edgeFeats));
    g.insert(GraphKeys::ResourceFeatures, std::move(resourceFeats));
    g.insert(GraphKeys::ResourceConsumptions, std::move(resourceCons));
    g.insert(GraphKeys::EdgeIdx, std::move(edges));
    g.insert(GraphKeys::ResourceDependencies, std::move(resourceDeps));
    g.insert(GraphKeys::EdgeResourceRelations, std::move(edgeResRel));
    g.insert(GraphKeys::EdgePairMask, std::move(edgePairMask));
    return g;
}

TEST(nn_utils, graphs_equivalent_same) {
    auto graph = createSimpleGraph();
    EXPECT_TRUE(tempo::nn::util::graphsEquivalent(graph, graph));
    EXPECT_TRUE(tempo::nn::util::compareGraphs(graph, graph));
}

TEST(nn_utils, graphs_equivalent_index_permutations) {
    using namespace tempo::nn;
    auto graph = createSimpleGraph();
    auto graphPerm = createSimpleGraph();
    graphPerm.at(GraphKeys::EdgeIdx) = torch::tensor({{1, 0, 1, 2}, {0, 1, 2, 0}}, indexTensorOptions());
    graphPerm.at(GraphKeys::EdgePairMask) = torch::tensor({0, 0, 1, 2}, indexTensorOptions());
    graphPerm.at(GraphKeys::EdgeFeatures) = torch::tensor({{0.9, -0.4}, {0.2, 0.2}, {0.4, -0.1}, {0.9, 0.45}},
                                                          dataTensorOptions());
    graphPerm.at(GraphKeys::EdgeResourceRelations) = torch::tensor({{3, 0, 1, 2}, {1, 1, 1, 1}}, indexTensorOptions());
    graphPerm.at(GraphKeys::ResourceDependencies) = torch::tensor({{0, 1, 2, 0}, {1, 1, 1, 0}}, indexTensorOptions());
    graphPerm.at(GraphKeys::ResourceConsumptions) = torch::tensor({{0.4}, {0.7}, {0.5}, {0.1}}, dataTensorOptions());
    EXPECT_TRUE(util::graphsEquivalent(graph, graphPerm));
    EXPECT_TRUE(util::compareGraphs(graph, graphPerm));
}

void swapTestRestore(const tempo::nn::InputGraph &gt, tempo::nn::InputGraph &graph,
                     const char *field, torch::Tensor t) {
    using namespace tempo::nn;
    auto orig = graph.at(field).clone();
    graph.insert_or_assign(field, std::move(t));
    EXPECT_FALSE(util::graphsEquivalent(gt, graph));
    EXPECT_FALSE(util::compareGraphs(gt, graph));
    graph.insert_or_assign(field, std::move(orig));
}

TEST(nn_utils, graphs_equivalent_not_equiv_feat) {
    using namespace tempo::nn;
    auto graph = createSimpleGraph();
    auto graph1 = createSimpleGraph();
    swapTestRestore(graph, graph1, GraphKeys::TaskFeatures,
                    torch::tensor({{0.1, 0.4}, {-1.4, 0.4}, {0.5, -0.1}}, dataTensorOptions()));
    swapTestRestore(graph, graph1, GraphKeys::EdgeFeatures,
                    torch::tensor({{0.9, -0.4}, {0.2, 0.2}, {0.9, 0.45}}, dataTensorOptions()));
    swapTestRestore(graph, graph1, GraphKeys::ResourceFeatures,
                    torch::tensor({{0.2}, {0.1}}, dataTensorOptions()));
    swapTestRestore(graph, graph1, GraphKeys::ResourceConsumptions,
                    torch::tensor({{0.1}, {0.7}, {0.5}}, dataTensorOptions()));

}

TEST(nn_utils, graphs_equivalent_not_equiv_indices) {
    using namespace tempo::nn;
    auto graph = createSimpleGraph();
    auto graph1 = createSimpleGraph();
    swapTestRestore(graph, graph1, GraphKeys::EdgeIdx,
                    torch::tensor({{0, 1, 2, 1}, {1, 0, 0, 0}}, indexTensorOptions()));
    swapTestRestore(graph, graph1, GraphKeys::ResourceDependencies,
                    torch::tensor({{0, 0, 1, 2}, {0, 1, 1, 0}}, indexTensorOptions()));
    swapTestRestore(graph, graph1, GraphKeys::EdgeResourceRelations,
                    torch::tensor({{0, 1, 2, 3}, {1, 0, 1, 1}}, indexTensorOptions()));
    swapTestRestore(graph, graph1, GraphKeys::EdgePairMask,
                    torch::tensor({0, 0, 1, 2}, indexTensorOptions()));
}