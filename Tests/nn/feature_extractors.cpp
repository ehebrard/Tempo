/**
 * @author Tim Luchterhand
 * @date 27.06.23.
 */

#include <gtest/gtest.h>
#include <Iterators.hpp>

#include "testing.hpp"
#include "nn/feature_extractors.hpp"
#include "util/Matrix.hpp"
#include "nn/GraphBuilder.hpp"
#include "util/SchedulingProblemHelper.hpp"


auto getInput() -> std::pair<tempo::testing::ProblemInstance, tempo::nn::Topology> {
    using namespace tempo;
    auto [problem, _] = tempo::testing::createTestProblem();
    nn::MinimalTopologyBuilder topologyBuilder(problem);
    auto topology = topologyBuilder.getTopology();
    return {std::move(problem), std::move(topology)};
}

TEST(nn_feature_extractors, TaskTimingExtractor) {
    using namespace tempo::nn;
    using namespace tempo;
    using namespace tempo::testing;
    using T = tempo::testing::TaskSpec;
    using tempo::testing::random_int;
    constexpr DataType Ub = 10;
    TaskTimingFeatureExtractor extractor;
    auto [tasks, scheduler] = createTasks(
            {T{.minDur = 2, .maxDur = 4, .earliestStart = 4, .latestDeadline = static_cast<int>(Ub)},
             T{.minDur = 6, .maxDur = 6, .earliestStart = 0, .latestDeadline = 5},
             T{.minDur = 4, .maxDur = 6, .earliestStart = 2, .latestDeadline = 8},
             T{.minDur = 0, .maxDur = static_cast<int>(Ub), .latestDeadline = static_cast<int>(Ub)}});
    auto sched = tasks.back();
    tasks.pop_back();
    ProblemInstance instance(std::move(tasks), {}, {}, sched);
    auto topology = MinimalTopologyBuilder(instance).getTopology();
    auto taskFeatures = extractor(topology, makeSolverState(Matrix<int>{}, scheduler), instance);
    for (auto [idx, task] : iterators::const_enumerate(instance.tasks())) {
        auto features = taskFeatures.slice(0, idx, idx + 1);
        EXPECT_EQ(features.numel(), 4);
        EXPECT_EQ(features[0][0].item<DataType>(), task.minDuration(scheduler) / Ub);
        EXPECT_EQ(features[0][1].item<DataType>(), task.maxDuration(scheduler) / Ub);
        EXPECT_EQ(features[0][2].item<DataType>(), task.getEarliestStart(scheduler) / Ub);
        EXPECT_EQ(features[0][3].item<DataType>(), task.getLatestEnd(scheduler) / Ub);
    }
}

TEST(nn_feature_extractors, ResourceEnergyExtractor) {
    using namespace tempo::nn;
    using namespace tempo::serialization;
    using namespace tempo;
    using namespace tempo::testing;
    auto [tasks, scheduler] = tempo::testing::createTasks({{6, 6}, {2, 4}, {1, 5}, {0, 12, 0, 12}});
    auto sched = tasks.back();
    tasks.pop_back();
    ProblemInstance instance(tasks, {{1, tasks, {1, 1, 1}}}, {}, sched);
    const auto topology = MinimalTopologyBuilder(instance).getTopology();
    ResourceEnergyExtractor extractor;
    auto resEnergies = extractor(topology, makeSolverState(Matrix<int>{}, scheduler), instance);
    ASSERT_EQ(resEnergies.size(0), 1);
    ASSERT_EQ(resEnergies.size(1), 1);
    EXPECT_EQ(resEnergies.item<DataType>(), 1);
}

TEST(nn_features_extractors, ResourceEnergyExtractor_complex_consumptions) {
    using namespace tempo::nn;
    using namespace tempo::serialization;
    using namespace tempo;
    using namespace tempo::testing;
    auto [tasks, scheduler] = tempo::testing::createTasks({{4, 4}, {2, 4}, {8, 10}, {0, 9, 0, 9}});
    auto sched = tasks.back();
    tasks.pop_back();
    ProblemInstance instance(tasks, {{4, tasks, {3, 2, 2}}}, {}, sched);
    const auto topology = MinimalTopologyBuilder(instance).getTopology();
    ResourceEnergyExtractor extractor;
    auto resEnergies = extractor(topology, makeSolverState(Matrix<int>{}, scheduler), instance);
    ASSERT_EQ(resEnergies.size(0), 1);
    ASSERT_EQ(resEnergies.size(1), 1);
    EXPECT_EQ(resEnergies.item<DataType>(), 1);
}

TEST(nn_feature_extractors, TimingEdgeExtractor) {
    using namespace tempo::nn;
    using namespace tempo::testing;
    using namespace tempo;
    tempo::Matrix<int> timings(3, 3, {0, 4, 0,
                                      3, 0, 5,
                                      0, -2, 0});
    auto [schedule, scheduler] = createTasks({{0, 4, 0, 4}});
    Topology topology;
    EdgeVector edges{{0, 1}, {2, 1}};
    topology.edgeIndices = util::makeIndexTensor(edges);
    torch::Tensor features = TimingEdgeExtractor()(topology, makeSolverState(timings, scheduler),
                                                   ProblemInstance({}, {}, {}, schedule.front()));
    ASSERT_EQ(features.size(0), edges.size());
    ASSERT_EQ(features.size(1), 2);
    ASSERT_EQ(features.sizes().size(), 2);
    EXPECT_EQ(features[0][0].item<float>(), 0);
    EXPECT_EQ(features[0][1].item<float>(), 1);
    EXPECT_EQ(features[1][0].item<float>(), 1);
    EXPECT_EQ(features[1][1].item<float>(), -0.5);
}