/**
 * @author Tim Luchterhand
 * @date 17.11.23.
 */

#ifndef TEMPO_TOPOLOGYBUILDERTEST_HPP
#define TEMPO_TOPOLOGYBUILDERTEST_HPP

#include <gtest/gtest.h>

#include "testing.hpp"
#include "nn/torch_types.hpp"
#include "util/Matrix.hpp"
#include "Model.hpp"
#include "util/SchedulingProblemHelper.hpp"

class TopologyBuilderTest : public testing::Test {
    tempo::testing::ProblemInstance inst;
    tempo::testing::ProblemInstance extInst;

protected:
    void SetUp() override;

    [[nodiscard]] auto instance() const noexcept -> const tempo::testing::ProblemInstance&;

    [[nodiscard]] auto extendedInstance() const noexcept -> const tempo::testing::ProblemInstance&;

    void testResourceDependencies(const tempo::nn::Topology &topology) const;

    void testEdges(const tempo::nn::Topology &topology) const;

    void testEdgePairMask(const tempo::nn::Topology &topology) const;

    [[nodiscard]] static auto getResourceMatrix(const tempo::testing::ProblemInstance &probInstance) -> tempo::Matrix<int>;

public:
    static void testEdgeResourceRelations(const torch::Tensor &edgeResourceRelations,
                                          const tempo::Matrix<int> &resConsumptions,
                                          const tempo::nn::EdgeVector &edges);

};


#endif //TEMPO_TOPOLOGYBUILDERTEST_HPP
