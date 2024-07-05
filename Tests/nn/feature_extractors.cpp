/**
 * @author Tim Luchterhand
 * @date 27.06.23.
 */

#include <gtest/gtest.h>
#include <Iterators.hpp>

#include "testing.hpp"
//#include "nn/feature_extractors.hpp"
#include "util/Matrix.hpp"
#include "nn/GraphBuilder.hpp"
#include "util/SchedulingProblemHelper.hpp"
#include "Model.hpp"


auto getInput() -> std::pair<tempo::testing::ProblemInstance, tempo::nn::Topology> {
    using namespace tempo;
    auto problem = tempo::testing::createTestProblem();
    nn::MinimalTopologyBuilder topologyBuilder(problem);
    auto topology = topologyBuilder.getTopology();
    return {std::move(problem), std::move(topology)};
}
/*
TEST(nn_feature_extractors, TaskTimingExtractor) {
    using namespace tempo::nn;
    using namespace tempo;
    using tempo::testing::random_int;
    using tempo::testing::setTaskDurations;
    using tempo::testing::setUpperBound;
    constexpr bool SatisfiesConcept = feature_extractor<TaskTimingFeatureExtractor, Matrix<int>>;
    constexpr DataType Ub = 10;
    EXPECT_TRUE(SatisfiesConcept);
    TaskTimingFeatureExtractor extractor;
    const auto [instance, topology] = getInput();
    const auto numEvents = instance.durations.size() * 2 + 2;
    Matrix<int> eventNet(numEvents, numEvents);
    std::vector<tempo::testing::TaskSpec> taskSpecs;
    for (std::size_t t = 0; t < instance.durations.size(); ++t) {
        taskSpecs.emplace_back(random_int(1, 5), random_int(5, 10), random_int(-10, 0), random_int(-10, 0));
        setTaskDurations(static_cast<tempo::task>(t), taskSpecs.back().minDur, taskSpecs.back().maxDur,
                         taskSpecs.back().release, taskSpecs.back().deadline, eventNet);
    }

    setUpperBound<int>(Ub, eventNet);
    auto taskFeatures = extractor(topology, eventNet);
    for (auto [idx, task] : iterators::const_enumerate(taskSpecs)) {
        auto features = taskFeatures.slice(0, idx, idx + 1);
        EXPECT_EQ(features.numel(), 4);
        EXPECT_EQ(features[0][0].item<DataType>(), task.minDur/ Ub);
        EXPECT_EQ(features[0][1].item<DataType>(), task.maxDur / Ub);
        EXPECT_EQ(features[0][2].item<DataType>(), task.release / Ub);
        EXPECT_EQ(features[0][3].item<DataType>(), task.deadline / Ub);
    }
}

TEST(nn_feature_extractors, ResourceEnergyExtractor) {
    using namespace tempo::nn;
    using namespace tempo::serialization;
    using namespace tempo;
    using namespace tempo::testing;
    constexpr bool SatisfiesConcept = feature_extractor<ResourceEnergyExtractor, Matrix<int>>;
    EXPECT_TRUE(SatisfiesConcept);
    ProblemInstance instance{.lowerBound = 0, .optimalSolution = 0, .durations = {1, 1, 1},
                             .constraints = {}, .resources = {{{0, 1, 2}, {1, 1, 1}, {}, 1}}};
    const auto topology = MinimalTopologyBuilder(instance).getTopology();
    Matrix<int> eventNet(8, 8);
    setTaskDurations(0, 6, 6, 0, 0, eventNet);
    setTaskDurations(1, 2, 4, 0, 0, eventNet);
    setTaskDurations(2, 1, 5, 0, 0, eventNet);
    setUpperBound(12, eventNet);
    ResourceEnergyExtractor extractor;
    auto resEnergies = extractor(topology, eventNet);
    ASSERT_EQ(resEnergies.size(0), 1);
    ASSERT_EQ(resEnergies.size(1), 1);
    EXPECT_EQ(resEnergies.item<DataType>(), 1);
}

TEST(nn_features_extractors, ResourceEnergyExtractor_complex_consumptions) {
    using namespace tempo::nn;
    using namespace tempo::serialization;
    using namespace tempo;
    using namespace tempo::testing;
    ProblemInstance instance{.lowerBound = 0, .optimalSolution = 0, .durations = {1, 1, 1},
            .constraints = {}, .resources = {{{0, 1, 2}, {3, 2, 2}, {}, 4}}};
    const auto topology = MinimalTopologyBuilder(instance).getTopology();
    Matrix<int> eventNet(8, 8);
    setTaskDurations(0, 4, 4, 0, 0, eventNet);
    setTaskDurations(1, 2, 4, 0, 0, eventNet);
    setTaskDurations(2, 8, 10, 0, 0, eventNet);
    setUpperBound(9, eventNet);
    ResourceEnergyExtractor extractor;
    auto resEnergies = extractor(topology, eventNet);
    ASSERT_EQ(resEnergies.size(0), 1);
    ASSERT_EQ(resEnergies.size(1), 1);
    EXPECT_EQ(resEnergies.item<DataType>(), 1);
}

TEST(nn_feature_extractors, TimingEdgeExtractor) {
    using namespace tempo::nn;
    using namespace tempo::testing;
    using namespace tempo;
    constexpr bool SatisfiesConcept = feature_extractor<TimingEdgeExtractor, Matrix<int>>;
    EXPECT_TRUE(SatisfiesConcept);
    tempo::Matrix<int> timings(8, 8);
    setTaskDistance(0, 1, 4, timings);
    setTaskDistance(1, 0, 3, timings);
    setTaskDistance(2, 1, -2, timings);
    setTaskDistance(1, 2, 5, timings);
    setUpperBound(4, timings);
    timings.at(1, 2) = 5;
    Topology topology;
    EdgeVector edges{{0, 1}, {2, 1}};
    topology.edgeIndices = util::makeIndexTensor(edges);
    torch::Tensor features = TimingEdgeExtractor()(topology, timings);
    ASSERT_EQ(features.size(0), edges.size());
    ASSERT_EQ(features.size(1), 2);
    ASSERT_EQ(features.sizes().size(), 2);
    EXPECT_EQ(features[0][0].item<float>(), 0);
    EXPECT_EQ(features[0][1].item<float>(), 1);
    EXPECT_EQ(features[1][0].item<float>(), 1);
    EXPECT_EQ(features[1][1].item<float>(), -0.5);
}*/