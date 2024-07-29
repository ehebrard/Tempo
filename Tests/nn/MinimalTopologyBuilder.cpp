/**
 * @author Tim Luchterhand
 * @date 17.11.23.
 */

#include <gtest/gtest.h>
#include <span>
#include <ranges>
#include <unordered_set>
#include <boost/container_hash/hash.hpp>


#include "testing.hpp"
#include "TopologyBuilderTest.hpp"
#include "nn/torch_types.hpp"
#include "nn/tensor_utils.hpp"
#include "nn/topology_extractors.hpp"

struct TestMinimalTopologyBuilder : public tempo::nn::MinimalTopologyBuilder {
    using tempo::nn::MinimalTopologyBuilder::MinimalTopologyBuilder;
    using tempo::nn::MinimalTopologyBuilder::getTopology;
    using tempo::nn::MinimalTopologyBuilder::addPrecedenceEdges;
    using tempo::nn::MinimalTopologyBuilder::addEdge;
    using tempo::nn::MinimalTopologyBuilder::completeSubGraph;
};

using tempo::nn::Edge;
using tempo::nn::EdgeVector;
using tempo::nn::IndexType;
using EdgeSet = std::unordered_set<Edge, boost::hash<Edge>>;
using EdgeMap = tempo::nn::impl::EdgeLookup;

TEST_F(TopologyBuilderTest, MinimalBuilder_resource_dependencies) {
    using namespace tempo::nn;
    MinimalTopologyBuilder topologyBuilder(instance());
    const auto& topology = topologyBuilder.getTopology();
    testResourceDependencies(topology);
}

TEST_F(TopologyBuilderTest, MinimalBuilder_edges) {
    using namespace tempo::nn;
    MinimalTopologyBuilder topologyBuilder(instance());
    const auto& topology = topologyBuilder.getTopology();
    testEdges(topology);
}

TEST_F(TopologyBuilderTest, MinimalBuilder_edge_mask) {
    using namespace tempo::nn;
    MinimalTopologyBuilder topologyBuilder(instance());
    const auto& topology = topologyBuilder.getTopology();
    testEdgePairMask(topology);
}

TEST_F(TopologyBuilderTest, MinimalBuilder_edge_resource_relations) {
    tempo::nn::MinimalTopologyBuilder topologyBuilder(instance());
    auto topology = topologyBuilder.getTopology();
    auto edgeView = tempo::nn::util::getEdgeView(topology.edgeIndices);
    std::vector edges(edgeView.begin(), edgeView.end());
    testEdgeResourceRelations(topology.edgeResourceRelations, getResourceMatrix(instance()), edges);
    topologyBuilder = tempo::nn::MinimalTopologyBuilder(extendedInstance());
    topology = topologyBuilder.getTopology();
    testEdgeResourceRelations(topology.edgeResourceRelations, getResourceMatrix(extendedInstance()), edges);
}

TEST(nn_MinimalTopologyBuilder, MinimalBuilder_addEdge) {
    using tempo::nn::Edge;
    tempo::nn::impl::TopologyData data{.edgeLookup = tempo::nn::impl::EdgeLookup(100)};
    TestMinimalTopologyBuilder::addEdge({3, 2}, false, 19, data);
    ASSERT_EQ(data.edges.size(), 1);
    ASSERT_EQ(data.edgePairMask.size(), 1);
    EXPECT_EQ(data.edges.back(), Edge(3, 2));
    EXPECT_EQ(data.edgePairMask.back(), 19);
    EXPECT_TRUE(data.edgeIdx.empty());
    ASSERT_EQ(data.edgeLookup.getNumEdges(), 1);
    EXPECT_EQ(data.edgeLookup.at(3, 2), 0);

    TestMinimalTopologyBuilder::addEdge({1, 4}, true, 17, data);
    ASSERT_EQ(data.edges.size(), 2);
    ASSERT_EQ(data.edgePairMask.size(), 2);
    EXPECT_EQ(data.edges.back(), Edge(1, 4));
    EXPECT_EQ(data.edgePairMask.back(), 17);
    EXPECT_EQ(data.edgeIdx.size(), 1);
    EXPECT_EQ(data.edgeIdx.back(), 1);
    EXPECT_EQ(data.edgeLookup.getNumEdges(), 2);
    EXPECT_EQ(data.edgeLookup.at(1, 4), 1);
}

TEST(nn_MinimalTopologyBuilder, MinimalBuilder_addPrecedenceEdges) {
    using tempo::nn::Edge;
    using iterators::const_zip_enumerate;
    using namespace tempo;
    using std::views::transform;
    tempo::nn::impl::TopologyData data{.edgeLookup = tempo::nn::impl::EdgeLookup(100)};
    auto tasks = tempo::testing::createDummyTasks(5);
    VarTaskMapping mapping(tasks);
    auto evtViewer = transform(
            [&tasks](const auto &dc) {
                return DistanceConstraint<int>(tasks.at(dc.from).start.id(), tasks.at(dc.to).end.id(), dc.distance);
            });
    std::vector<DistanceConstraint<int>> edges{{1, 2, 0}, {2, 3, 2}, {3, 4, -1}};
    TestMinimalTopologyBuilder::addEdge({2, 5}, true, -1, data);
    TestMinimalTopologyBuilder::addPrecedenceEdges(edges | evtViewer, mapping, data);
    ASSERT_EQ(data.edges.size(), 4);
    EXPECT_EQ(data.edgeIdx.size(), 1);
    EXPECT_EQ(data.edgeIdx.front(), 0);
    EdgeSet gtEdges{{2, 5}, {1, 2}, {2, 3}, {3, 4}};
    EXPECT_EQ(gtEdges, EdgeSet(data.edges.begin(), data.edges.end()));
    for (auto [idx, e, m] : const_zip_enumerate(data.edges, data.edgePairMask)) {
        EXPECT_EQ(data.edgeLookup.at(e.first, e.second), idx);
        EXPECT_EQ(m, idx - 1);
    }
}


auto edgeFromTensor(const torch::Tensor &tensor) -> tempo::nn::Edge {
    using tempo::nn::Edge;
    using tempo::nn::IndexType;
    using namespace torch::indexing;
    assert(tensor.size(0) == 2);
    assert(tensor.sizes().size() == 1);
    return {tensor[0].item<IndexType>(), tensor[1].item<IndexType>()};
}

auto edgesFromTensor(const torch::Tensor &tensor) -> std::vector<tempo::nn::Edge> {
    using tempo::nn::Edge;
    using tempo::nn::IndexType;
    using namespace torch::indexing;
    assert(tensor.size(0) == 2);
    std::vector<Edge> ret;
    ret.reserve(tensor.size(1));
    for (long i = 0; i < tensor.size(1); ++i) {
        ret.emplace_back(tensor[0][i].item<IndexType>(), tensor[1][i].item<IndexType>());
    }

    return ret;
}

template<typename T>
auto tensorFromVec(std::vector<T> &vec) -> torch::Tensor {
    return torch::from_blob(vec.data(), static_cast<long>(vec.size()), torch::dtype<T>());
}

void testEdges(const EdgeVector &edges, const EdgeSet &gtEdges, const EdgeMap &edgeLookup,
               std::vector<IndexType> &edgePair) {
    using namespace torch::indexing;
    using namespace tempo::nn;
    // test edges
    EdgeSet res(edges.begin(), edges.end());
    auto edgeTensor = util::makeIndexTensor(edges);
    auto maskTensor = tensorFromVec(edgePair);
    EXPECT_EQ(res, gtEdges);
    for (const auto &e : gtEdges) {
        ASSERT_TRUE(edgeLookup.contains(e));
        auto edgeIdx = edgeLookup.at(e.first, e.second);
        auto edge = edgeFromTensor(edgeTensor.index({Slice{None}, edgeIdx}));
        EXPECT_EQ(edge, e);

    }
    // test pair mask
    const auto maxMask = maskTensor.max().item<IndexType>();
    for (auto m = 0; m <= maxMask; ++m) {
        auto slice = edgeTensor.index({Slice{None}, maskTensor == m});
        auto resEdges = edgesFromTensor(slice);
        ASSERT_EQ(resEdges.size(), 2);
        EXPECT_EQ(resEdges.front().first, resEdges.back().second);
        EXPECT_EQ(resEdges.front().second, resEdges.back().first);
    }
}

TEST(nn_MinimalTopologyBuilder, MinimalBuilder_completeSubGraphOneResource) {
    using namespace tempo;
    std::vector<unsigned> taskIds{1, 2, 3};
    std::vector<int> demands{2, 1, 4};
    constexpr int Capacity = 4;
    tempo::testing::Resource resource(Capacity, taskIds, demands);
    tempo::nn::impl::TopologyData data{.edgeLookup = tempo::nn::impl::EdgeLookup(100)};
    TestMinimalTopologyBuilder::completeSubGraph(resource, 17, data);
    ASSERT_EQ(data.taskIdx.size(), 3);
    ASSERT_EQ(data.resIdx.size(), 3);
    ASSERT_EQ(data.resDemands.size(), 3);
    ASSERT_EQ(data.edges.size(), 6);
    ASSERT_EQ(data.edgePairMask.size(), 6);
    ASSERT_EQ(data.edgeRelResIdx.size(), 6);
    ASSERT_EQ(data.edgeIdx.size(), 6);
    ASSERT_EQ(data.edgeLookup.getNumEdges(), 6);
    EdgeSet gtEdges {{1, 2}, {2, 1}, {1, 3}, {3, 1}, {2, 3}, {3, 2}};
    testEdges(data.edges, gtEdges, data.edgeLookup, data.edgePairMask);

    // test task resource relation
    for (auto r : data.resIdx) {
        EXPECT_EQ(r, 17);
    }

    for (auto [t, gtDemand] : iterators::const_zip(taskIds, demands)) {
        auto res = std::ranges::find(data.taskIdx, t);
        ASSERT_NE(res, data.taskIdx.end());
        auto idx = static_cast<std::size_t>(res - data.taskIdx.begin());
        EXPECT_EQ(data.resDemands.at(idx), gtDemand / static_cast<tempo::nn::DataType>(Capacity));
    }

    // test edge resource relation
    for (auto r : data.edgeRelResIdx) {
        EXPECT_EQ(r, 17);
    }

    using iterators::const_enumerate;
    for (auto [idx, edgeIdx] : const_enumerate(data.edgeIdx)) {
        EXPECT_EQ(idx, edgeIdx);
    }
}

void checkIdxTensor(std::vector<IndexType> &sourceIdx, std::vector<IndexType> &targetIdx, IndexType idxVal,
                    std::initializer_list<IndexType> expectedIndices) {

    const auto sourceIdxTensor = tensorFromVec(sourceIdx);
    const auto targetIdxTensor = tensorFromVec(targetIdx);
    auto targets = targetIdxTensor.index({sourceIdxTensor == idxVal});
    std::unordered_set<IndexType> tSet(targets.data_ptr<IndexType>(), targets.data_ptr<IndexType>() + targets.size(0));
    EXPECT_EQ(tSet, decltype(tSet)(expectedIndices));
}

void checkEdge(const Edge &edge, const EdgeVector &allEdges, std::vector<IndexType> &edgeIdx,
               std::vector<IndexType> &resourceIdx, std::initializer_list<IndexType> expectedIndices) {
    auto idx = static_cast<IndexType>(std::ranges::find(allEdges, edge) - allEdges.begin());
    checkIdxTensor(edgeIdx, resourceIdx, idx, expectedIndices);
}

TEST(nn_MinimalTopologyBuilder, MinimalBuilder_completeSubGraph_multiple_resources) {
    using namespace tempo;
    std::vector<unsigned> tasks17{1, 2};
    std::vector<int> demands17{2, 1};
    constexpr int Capacity17 = 2;
    tempo::nn::impl::TopologyData data{.edgeLookup = tempo::nn::impl::EdgeLookup(100)};
    TestMinimalTopologyBuilder::completeSubGraph(tempo::testing::Resource(Capacity17, tasks17, demands17), 17, data);
    decltype(tasks17) tasks18{1, 2, 3};
    decltype(demands17) demands18{2, 1, 1};
    constexpr int Capacity18 = 3;
    TestMinimalTopologyBuilder::completeSubGraph(
            tempo::testing::Resource(Capacity18, tasks18, demands18), 18, data);
    ASSERT_EQ(data.taskIdx.size(), 5);
    ASSERT_EQ(data.resIdx.size(), 5);
    ASSERT_EQ(data.resDemands.size(), 5);
    ASSERT_EQ(data.edges.size(), 6);
    ASSERT_EQ(data.edgePairMask.size(), 6);
    ASSERT_EQ(data.edgeIdx.size(), 8);
    ASSERT_EQ(data.edgeRelResIdx.size(), 8);
    EdgeSet gtEdges{{1, 2}, {2, 1}, {1, 3}, {3, 1}, {2, 3}, {3, 2}};
    testEdges(data.edges, gtEdges, data.edgeLookup, data.edgePairMask);

    checkIdxTensor(data.resIdx, data.taskIdx, 17, {1, 2});
    checkIdxTensor(data.resIdx, data.taskIdx, 18, {1, 2, 3});
    for (auto [t, r, d]: iterators::zip(data.taskIdx, data.resIdx,data.resDemands)) {
        const auto &gtT = r == 17 ? tasks17 : tasks18;
        const auto &gtD = r == 17 ? demands17 : demands18;
        auto cap = r == 17 ? Capacity17 : Capacity18;
        auto res = std::ranges::find(gtT, t);
        ASSERT_NE(res, gtT.end());
        auto idx = static_cast<std::size_t>(res - gtT.begin());
        EXPECT_EQ(d, gtD.at(idx) / static_cast<tempo::nn::DataType>(cap));
    }

    checkEdge({1, 2}, data.edges, data.edgeIdx, data.edgeRelResIdx, {17, 18});
    checkEdge({2, 1}, data.edges, data.edgeIdx, data.edgeRelResIdx, {17, 18});
    checkEdge({1, 3}, data.edges, data.edgeIdx, data.edgeRelResIdx, {18});
    checkEdge({3, 1}, data.edges, data.edgeIdx, data.edgeRelResIdx, {18});
    checkEdge({2, 3}, data.edges, data.edgeIdx, data.edgeRelResIdx, {18});
    checkEdge({3, 2}, data.edges, data.edgeIdx, data.edgeRelResIdx, {18});
}