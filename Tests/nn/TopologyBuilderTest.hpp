/**
 * @author Tim Luchterhand
 * @date 17.11.23.
 */

#ifndef TEMPO_TOPOLOGYBUILDERTEST_HPP
#define TEMPO_TOPOLOGYBUILDERTEST_HPP

#include <gtest/gtest.h>

#include "util/parsing/format.hpp"
#include "nn/torch_types.hpp"
#include "util/Matrix.hpp"

class TopologyBuilderTest : public testing::Test {
    ProblemInstance inst;
    ProblemInstance extInst;

protected:
    void SetUp() override;

    [[nodiscard]] auto instance() const noexcept -> const ProblemInstance&;

    [[nodiscard]] auto extendedInstance() const noexcept -> const ProblemInstance&;

    void testResourceDependencies(const tempo::nn::Topology &topology) const;

    void testEdges(const tempo::nn::Topology &topology) const;

    void testEdgePairMask(const tempo::nn::Topology &topology) const;

    [[nodiscard]] static auto getResourceMatrix(const ProblemInstance &probInstance) -> tempo::Matrix<int>;

public:
    static void testEdgeResourceRelations(const torch::Tensor &edgeResourceRelations,
                                          const tempo::Matrix<int> &resConsumptions,
                                          const tempo::nn::EdgeVector &edges);

};


#endif //TEMPO_TOPOLOGYBUILDERTEST_HPP
