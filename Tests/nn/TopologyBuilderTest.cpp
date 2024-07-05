/**
 * @author Tim Luchterhand
 * @data 17.11.23.
 */

#include <span>
#include <ranges>
#include <Iterators.hpp>

#include "TopologyBuilderTest.hpp"
#include "nn/tensor_utils.hpp"
#include "util/serialization.hpp"
#include "Tests/testing.hpp"

void TopologyBuilderTest::SetUp() {
    inst = tempo::testing::createTestProblem();
    extInst = tempo::testing::createTestProblem(); //@TODO
}

auto TopologyBuilderTest::instance() const noexcept -> const tempo::testing::ProblemInstance & {
    return inst;
}

auto TopologyBuilderTest::extendedInstance() const noexcept -> const tempo::testing::ProblemInstance & {
    return extInst;
}

void TopologyBuilderTest::testEdgeResourceRelations(const torch::Tensor &edgeResourceRelations,
                                                    const tempo::Matrix<int> &resConsumptions,
                                                    const tempo::nn::EdgeVector &edges) {
    ASSERT_EQ(edgeResourceRelations.sizes().size(), 2);
    ASSERT_EQ(edgeResourceRelations.size(0), 2);
    for (int i = 0; i < edgeResourceRelations.size(1); ++i) {
        const auto e = edgeResourceRelations[0][i].item<long>();
        const auto r = edgeResourceRelations[1][i].item<long>();
        const auto [src, dest] = edges.at(e);
        EXPECT_GT(resConsumptions.at(r, src), 0);
        EXPECT_GT(resConsumptions.at(r, dest), 0);
    }

    const auto edgeView = edgeResourceRelations[0];
    const auto resView = edgeResourceRelations[1];
    for (auto [idx, edge] : iterators::enumerate(edges, 0l)) {
        auto resources = resView.index({edgeView == idx});
        for (auto r = 0ul; r < resConsumptions.numRows(); ++r) {
            if (resConsumptions.at(r, edge.first) > 0 and resConsumptions.at(r, edge.second) > 0) {
                EXPECT_EQ(torch::sum(resources == static_cast<long>(r)).item<long>(), 1);
            }
        }
    }
}

void TopologyBuilderTest::testResourceDependencies(const tempo::nn::Topology &topology) const {
    using namespace tempo::nn;
    ASSERT_EQ(topology.resourceDependencies.size(0), 2);
    ASSERT_EQ(topology.resourceDemands.sizes().size(), 2);
    ASSERT_EQ(topology.resourceDependencies.size(1), topology.resourceDemands.size(0));
    EXPECT_EQ(topology.resourceDemands.size(1), 1);
    const auto resMat = getResourceMatrix(instance());
    const auto tasks = topology.resourceDependencies.slice(0, 0, 1).flatten().contiguous();
    const auto resources = topology.resourceDependencies.slice(0, 1, 2).flatten().contiguous();
    ASSERT_TRUE(topology.resourceDemands.is_contiguous());
    auto zRange = iterators::const_zip(std::span(tasks.data_ptr<IndexType>(), tasks.size(0)),
                                       std::span(resources.data_ptr<IndexType>(), tasks.size(0)),
                                       std::span(topology.resourceDemands.data_ptr<DataType>(), tasks.size(0)));
    for (auto [tIdx, rIdx, d] : zRange) {
        EXPECT_GT(resMat.at(rIdx, tIdx), 0);
        EXPECT_EQ(static_cast<DataType>(resMat.at(rIdx, tIdx)) / instance().resources()[rIdx].resourceCapacity(), d);
    }

    using std::views::iota;
    using std::views::filter;
    for (auto t = 0ul; t < instance().tasks().size(); ++t) {
        auto gtRange =
                iota(0ul, resMat.numRows()) | filter([&resMat, t](auto r) { return resMat.at(r, t) > 0; });
        std::unordered_set<std::size_t> gtResources(gtRange.begin(), gtRange.end());
        const auto resTensor = resources.index({tasks == static_cast<IndexType>(t)});
        ASSERT_TRUE(resTensor.is_contiguous());
        const auto *dataPtr = resTensor.data_ptr<IndexType>();
        std::unordered_set<std::size_t> res(dataPtr, dataPtr + resTensor.numel());
        EXPECT_EQ(gtResources, res);
    }
}

void TopologyBuilderTest::testEdges(const tempo::nn::Topology &topology) const {
    using namespace tempo::nn;
    const auto consumptions = getResourceMatrix(instance());
    auto edges = util::getEdgeView(topology.edgeIndices);
    EXPECT_EQ(edges.size(), 9); // only for this particular test problem
    for (auto [from, to] : edges) {
        auto res = std::ranges::find_if(instance().precedences(), [from, to, this](const auto &prec) {
            return instance().getMapping()(prec.from) == from and instance().getMapping()(prec.to) == to;
        });
        if (res != instance().precedences().end()) {
            continue;
        }

        bool edgeFound = false;
        for (auto r = 0ul; r < consumptions.numRows(); ++r) {
            if (consumptions.at(r, from) > 0 and consumptions(r, to) > 0) {
                edgeFound = true;
                break;
            }
        }

        EXPECT_TRUE(edgeFound);
    }
}

void TopologyBuilderTest::testEdgePairMask(const tempo::nn::Topology &topology) const {
    using namespace tempo::nn;
    using namespace tempo::testing;
    auto edges = util::getEdgeView(topology.edgeIndices);
    ASSERT_EQ(edges.size(), topology.edgePairMask.numel());
    ASSERT_TRUE(topology.edgePairMask.is_contiguous());
    const std::span pairMaskVec(topology.edgePairMask.data_ptr<IndexType>(), edges.size());
    const auto maxId = topology.edgePairMask.max().item<IndexType>();
    for (auto pairId = 0; pairId <= maxId; ++pairId) {
        auto numEdges = (topology.edgePairMask == pairId).sum().item<int>();
        ASSERT_TRUE(numEdges == 1 or numEdges == 2);
        if (numEdges == 1) {
            continue;
        }

        auto loc1 = std::ranges::find(pairMaskVec, pairId);
        auto loc2 = std::ranges::find(pairMaskVec | std::views::reverse, pairId);
        ASSERT_NE(loc1, pairMaskVec.end());
        ASSERT_NE(loc2, pairMaskVec.rend());
        auto edge1 = edges[loc1 - pairMaskVec.begin()];
        auto edge2 = edges[edges.size() - 1 - (loc2 - pairMaskVec.rbegin())];
        EXPECT_EQ(edge1.first, edge2.second);
        EXPECT_EQ(edge2.first, edge1.second);
    }

    for (auto [idx1, edge1] : iterators::enumerate(edges)) {
        for (auto [idx2, edge2] : iterators::enumerate(edges)) {
            if (edge1.first == edge2.second and edge2.first == edge1.second) {
                ASSERT_LT(idx1, pairMaskVec.size());
                ASSERT_LT(idx2, pairMaskVec.size());
                EXPECT_EQ(pairMaskVec[idx1], pairMaskVec[idx2]);
                break;
            }
        }
    }
}

auto TopologyBuilderTest::getResourceMatrix(const tempo::testing::ProblemInstance &probInstance) -> tempo::Matrix<int> {
    tempo::Matrix<int> ret(probInstance.resources().size(), probInstance.tasks().size());
    for (const auto [r, res] : iterators::const_enumerate(probInstance.resources())) {
        for (auto [idx, t] : iterators::enumerate(res)) {
            ret.at(r, probInstance.getMapping()(t.id())) = res.getDemand(idx);
        }
    }

    return ret;
}
