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
    TestMinimalTopologyBuilder() : tempo::nn::MinimalTopologyBuilder(100) {}
    using tempo::nn::MinimalTopologyBuilder::MinimalTopologyBuilder;
    using tempo::nn::MinimalTopologyBuilder::getTopology;
    using tempo::nn::MinimalTopologyBuilder::addPrecedenceEdges;
    using tempo::nn::MinimalTopologyBuilder::resetToImmutableTopology;
    using tempo::nn::MinimalTopologyBuilder::addEdge;
    using tempo::nn::MinimalTopologyBuilder::completeSubGraph;
    auto &getTaskIdx() noexcept { return this->taskIdx; }
    auto &getResIdx() noexcept { return this->resIdx; }
    auto &getResDemands() noexcept { return this->resDemands; }
    auto &getEdgeIdx() noexcept { return this->edgeIdx; }
    auto &getEdgeRelResIdx() noexcept { return this->edgeRelResIdx; }
    auto &getEdges() noexcept { return this->edges; }
    auto &getEdgePairMask() noexcept { return this->edgePairMask; }
    auto &getEdgeLookup() noexcept { return this->edgeLookup; }
    std::size_t edgeLookupSize() {
        std::size_t ret = 0;
        getEdgeLookup().for_each([&ret](auto val) { ret += (val != tempo::nn::impl::EdgeLookup::NoValue); });
        return ret;
    }
};

using tempo::nn::Edge;
using tempo::nn::EdgeVector;
using tempo::nn::IndexType;
using EdgeSet = std::unordered_set<Edge, boost::hash<Edge>>;
using EdgeMap = std::remove_cvref_t<decltype(std::declval<TestMinimalTopologyBuilder>().getEdgeLookup())>;

TEST_F(TopologyBuilderTest, MinimalBuilder_resource_dependencies) {
    using namespace tempo::nn;
    MinimalTopologyBuilder topologyBuilder(instance());
    auto topology = topologyBuilder.getTopology();
    testResourceDependencies(topology);
}

TEST_F(TopologyBuilderTest, MinimalBuilder_edges) {
    using namespace tempo::nn;
    MinimalTopologyBuilder topologyBuilder(instance());
    auto topology = topologyBuilder.getTopology();
    testEdges(topology);
}

TEST_F(TopologyBuilderTest, MinimalBuilder_edge_mask) {
    using namespace tempo::nn;
    MinimalTopologyBuilder topologyBuilder(instance());
    auto topology = topologyBuilder.getTopology();
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
    TestMinimalTopologyBuilder topologyBuilder;
    topologyBuilder.addEdge({3, 2}, true, false, 19);
    ASSERT_EQ(topologyBuilder.getEdges().size(), 1);
    ASSERT_EQ(topologyBuilder.getEdgePairMask().size(), 1);
    EXPECT_EQ(topologyBuilder.getEdges().back(), Edge(3, 2));
    EXPECT_EQ(topologyBuilder.getEdgePairMask().back(), 19);
    EXPECT_TRUE(topologyBuilder.getEdgeIdx().empty());
    ASSERT_EQ(topologyBuilder.edgeLookupSize(), 1);
    EXPECT_EQ(topologyBuilder.getEdgeLookup().at(3, 2), 0);

    topologyBuilder.addEdge({1, 4}, true, true, 17);
    ASSERT_EQ(topologyBuilder.getEdges().size(), 2);
    ASSERT_EQ(topologyBuilder.getEdgePairMask().size(), 2);
    EXPECT_EQ(topologyBuilder.getEdges().back(), Edge(1, 4));
    EXPECT_EQ(topologyBuilder.getEdgePairMask().back(), 17);
    EXPECT_EQ(topologyBuilder.getEdgeIdx().size(), 1);
    EXPECT_EQ(topologyBuilder.getEdgeIdx().back(), 1);
    EXPECT_EQ(topologyBuilder.edgeLookupSize(), 2);
    EXPECT_EQ(topologyBuilder.getEdgeLookup().at(1, 4), 1);

    topologyBuilder.addEdge({8, 3}, false, true, 7);
    ASSERT_EQ(topologyBuilder.getEdges().size(), 3);
    ASSERT_EQ(topologyBuilder.getEdgePairMask().size(), 3);
    EXPECT_EQ(topologyBuilder.getEdges().back(), Edge(8, 3));
    EXPECT_EQ(topologyBuilder.getEdgePairMask().back(), 7);
    ASSERT_EQ(topologyBuilder.getEdgeIdx().size(), 2);
    EXPECT_EQ(topologyBuilder.getEdgeIdx().back(), 2);
    EXPECT_EQ(topologyBuilder.edgeLookupSize(), 3);
    EXPECT_EQ(topologyBuilder.getEdgeLookup().at(8, 3), 2);

    topologyBuilder.addEdge({2, 3}, false, false, 5);
    ASSERT_EQ(topologyBuilder.getEdges().size(), 4);
    ASSERT_EQ(topologyBuilder.getEdgePairMask().size(), 4);
    EXPECT_EQ(topologyBuilder.getEdges().back(), Edge(2, 3));
    EXPECT_EQ(topologyBuilder.getEdgePairMask().back(), 5);
    ASSERT_EQ(topologyBuilder.getEdgeIdx().size(), 2);
    EXPECT_EQ(topologyBuilder.edgeLookupSize(), 3);

}

TEST(nn_MinimalTopologyBuilder, MinimalBuilder_addPrecedenceEdges) {
    using tempo::nn::Edge;
    using iterators::const_zip_enumerate;
    using namespace tempo;
    using std::views::transform;
    auto evtViewer = transform(
            [](const auto &dc) { return DistanceConstraint<int>(START(dc.from), END(dc.to), dc.distance); });
    std::vector<DistanceConstraint<int>> edges{{1, 2, 0}, {2, 3, 2}, {3, 4, -1}};
    TestMinimalTopologyBuilder topologyBuilder;
    topologyBuilder.addEdge({2, 5}, true, true, -1);
    topologyBuilder.addPrecedenceEdges(edges | evtViewer, true);
    ASSERT_EQ(topologyBuilder.getEdges().size(), 4);
    EXPECT_EQ(topologyBuilder.getEdgeIdx().size(), 1);
    EXPECT_EQ(topologyBuilder.getEdgeIdx().front(), 0);
    EdgeSet gtEdges{{2, 5}, {1, 2}, {2, 3}, {3, 4}};
    EXPECT_EQ(gtEdges, EdgeSet(topologyBuilder.getEdges().begin(), topologyBuilder.getEdges().end()));
    for (auto [idx, e, m] : const_zip_enumerate(topologyBuilder.getEdges(), topologyBuilder.getEdgePairMask())) {
        EXPECT_EQ(topologyBuilder.getEdgeLookup().at(e.first, e.second), idx);
        EXPECT_EQ(m, idx - 1);
    }

    edges = decltype(edges){{4, 5, -2}, {5, 6, -1}};
    topologyBuilder.addPrecedenceEdges(edges | evtViewer, false);
    ASSERT_EQ(topologyBuilder.getEdges().size(), 6);
    EXPECT_EQ(topologyBuilder.getEdgeIdx().size(), 1);
    EXPECT_EQ(topologyBuilder.edgeLookupSize(), 4);
    EXPECT_EQ(topologyBuilder.getEdges()[4], Edge(4, 5));
    EXPECT_EQ(topologyBuilder.getEdges()[5], Edge(5, 6));
}

TEST(nn_MinimalTopologyBuilder, MinimalBuilder_resetToImmutableTopology) {
    using tempo::nn::Edge;
    using iterators::const_zip_enumerate;
    using namespace tempo;
    TestMinimalTopologyBuilder topologyBuilder;
    topologyBuilder.addEdge({1, 2}, true, true, 1);
    topologyBuilder.addEdge({2, 3}, true, false, 2);
    std::vector<DistanceConstraint<int>> precedences{{3, 4, -1}, {4, 5, 0}};
    topologyBuilder.addPrecedenceEdges(precedences, false);
    topologyBuilder.addEdge({17, 4}, false, false, 9);
    topologyBuilder.resetToImmutableTopology();
    ASSERT_EQ(topologyBuilder.getEdges().size(), 2);
    ASSERT_EQ(topologyBuilder.getEdgePairMask().size(), 2);
    for (auto [idx, e, m] : const_zip_enumerate(topologyBuilder.getEdges(), topologyBuilder.getEdgePairMask(), 1)) {
        EXPECT_EQ(e, Edge(idx, idx + 1));
        EXPECT_EQ(m, idx);
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
    TestMinimalTopologyBuilder topologyBuilder;
    std::vector<int> tasks{1, 2, 3};
    std::vector<int> demands{2, 1, 4};
    constexpr int Capacity = 4;
    Resource<int> resource(tasks, demands, {}, Capacity);
    topologyBuilder.completeSubGraph(resource, 17);
    topologyBuilder.resetToImmutableTopology();
    ASSERT_EQ(topologyBuilder.getTaskIdx().size(), 3);
    ASSERT_EQ(topologyBuilder.getResIdx().size(), 3);
    ASSERT_EQ(topologyBuilder.getResDemands().size(), 3);
    ASSERT_EQ(topologyBuilder.getEdges().size(), 6);
    ASSERT_EQ(topologyBuilder.getEdgePairMask().size(), 6);
    ASSERT_EQ(topologyBuilder.getEdgeRelResIdx().size(), 6);
    ASSERT_EQ(topologyBuilder.getEdgeIdx().size(), 6);
    ASSERT_EQ(topologyBuilder.edgeLookupSize(), 6);
    EdgeSet gtEdges {{1, 2}, {2, 1}, {1, 3}, {3, 1}, {2, 3}, {3, 2}};
    testEdges(topologyBuilder.getEdges(), gtEdges, topologyBuilder.getEdgeLookup(), topologyBuilder.getEdgePairMask());

    // test task resource relation
    for (auto r : topologyBuilder.getResIdx()) {
        EXPECT_EQ(r, 17);
    }

    for (auto [t, gtDemand] : iterators::const_zip(tasks, demands)) {
        auto res = std::ranges::find(topologyBuilder.getTaskIdx(), t);
        ASSERT_NE(res, topologyBuilder.getTaskIdx().end());
        auto idx = static_cast<std::size_t>(res - topologyBuilder.getTaskIdx().begin());
        EXPECT_EQ(topologyBuilder.getResDemands().at(idx), gtDemand / static_cast<tempo::nn::DataType>(Capacity));
    }

    // test edge resource relation
    for (auto r : topologyBuilder.getEdgeRelResIdx()) {
        EXPECT_EQ(r, 17);
    }

    using iterators::const_enumerate;
    for (auto [idx, edgeIdx] : const_enumerate(topologyBuilder.getEdgeIdx())) {
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
    TestMinimalTopologyBuilder topologyBuilder;
    std::vector<int> tasks17{1, 2};
    std::vector<int> demands17{2, 1};
    constexpr int Capacity17 = 2;
    topologyBuilder.completeSubGraph(Resource<int>(tasks17, demands17, {}, Capacity17), 17);
    decltype(tasks17) tasks18{1, 2, 3};
    decltype(demands17) demands18{2, 1, 1};
    constexpr int Capacity18 = 3;
    topologyBuilder.completeSubGraph(Resource<int>(tasks18, demands18, {}, Capacity18), 18);
    ASSERT_EQ(topologyBuilder.getTaskIdx().size(), 5);
    ASSERT_EQ(topologyBuilder.getResIdx().size(), 5);
    ASSERT_EQ(topologyBuilder.getResDemands().size(), 5);
    ASSERT_EQ(topologyBuilder.getEdges().size(), 6);
    ASSERT_EQ(topologyBuilder.getEdgePairMask().size(), 6);
    ASSERT_EQ(topologyBuilder.getEdgeIdx().size(), 8);
    ASSERT_EQ(topologyBuilder.getEdgeRelResIdx().size(), 8);
    EdgeSet gtEdges{{1, 2}, {2, 1}, {1, 3}, {3, 1}, {2, 3}, {3, 2}};
    testEdges(topologyBuilder.getEdges(), gtEdges, topologyBuilder.getEdgeLookup(), topologyBuilder.getEdgePairMask());

    checkIdxTensor(topologyBuilder.getResIdx(), topologyBuilder.getTaskIdx(), 17, {1, 2});
    checkIdxTensor(topologyBuilder.getResIdx(), topologyBuilder.getTaskIdx(), 18, {1, 2, 3});
    for (auto [t, r, d]: iterators::zip(topologyBuilder.getTaskIdx(), topologyBuilder.getResIdx(),
                                        topologyBuilder.getResDemands())) {
        const auto &gtT = r == 17 ? tasks17 : tasks18;
        const auto &gtD = r == 17 ? demands17 : demands18;
        auto cap = r == 17 ? Capacity17 : Capacity18;
        auto res = std::ranges::find(gtT, t);
        ASSERT_NE(res, gtT.end());
        auto idx = static_cast<std::size_t>(res - gtT.begin());
        EXPECT_EQ(d, gtD.at(idx) / static_cast<tempo::nn::DataType>(cap));
    }

    checkEdge({1, 2}, topologyBuilder.getEdges(), topologyBuilder.getEdgeIdx(), topologyBuilder.getEdgeRelResIdx(),
              {17, 18});
    checkEdge({2, 1}, topologyBuilder.getEdges(), topologyBuilder.getEdgeIdx(), topologyBuilder.getEdgeRelResIdx(),
              {17, 18});
    checkEdge({1, 3}, topologyBuilder.getEdges(), topologyBuilder.getEdgeIdx(), topologyBuilder.getEdgeRelResIdx(),
              {18});
    checkEdge({3, 1}, topologyBuilder.getEdges(), topologyBuilder.getEdgeIdx(), topologyBuilder.getEdgeRelResIdx(),
              {18});
    checkEdge({2, 3}, topologyBuilder.getEdges(), topologyBuilder.getEdgeIdx(), topologyBuilder.getEdgeRelResIdx(),
              {18});
    checkEdge({3, 2}, topologyBuilder.getEdges(), topologyBuilder.getEdgeIdx(), topologyBuilder.getEdgeRelResIdx(),
              {18});
}