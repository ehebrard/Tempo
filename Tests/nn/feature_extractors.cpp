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
    auto problem = tempo::testing::createTestProblem();
    nn::MinimalTopologyBuilder topologyBuilder(problem);
    auto topology = topologyBuilder.getTopology();
    return {std::move(problem), std::move(topology)};
}

struct DummyScheduler {
    tempo::testing::BoundProvider numeric;
    template<typename ...Args>
    explicit DummyScheduler(Args &&...args): numeric(std::forward<Args>(args)...) {}
};

TEST(nn_feature_extractors, TaskTimingExtractor) {
    using namespace tempo::nn;
    using namespace tempo;
    using tempo::testing::random_int;
    constexpr DataType Ub = 10;
    TaskTimingFeatureExtractor extractor;
    auto [instance, topology] = getInput();
    std::vector<int> upper(instance.tasks().size() * 2 + 2, -1);
    std::vector<int> lower(instance.tasks().size() * 2 + 2, -1);
    for (const auto &t : instance.tasks()) {
        lower.at(t.start) = random_int(-10, 0);
        upper.at(t.end) = random_int(-10, 0);
    }

    lower.at(instance.schedule().start) = 0;
    upper.at(instance.schedule().end) = Ub;
    DummyScheduler scheduler(std::move(upper), std::move(lower));

    auto taskFeatures = extractor(topology, makeSolverState(Matrix<int>{}, scheduler), instance);
    for (auto [idx, task] : iterators::const_enumerate(instance.tasks())) {
        auto features = taskFeatures.slice(0, idx, idx + 1);
        EXPECT_EQ(features.numel(), 4);
        EXPECT_EQ(features[0][0].item<DataType>(), task.minDuration() / Ub);
        EXPECT_EQ(features[0][1].item<DataType>(), task.maxDuration() / Ub);
        EXPECT_EQ(features[0][2].item<DataType>(), task.getEarliestStart(scheduler) / Ub);
        EXPECT_EQ(features[0][3].item<DataType>(), task.getLatestEnd(scheduler) / Ub);
    }
}

TEST(nn_feature_extractors, ResourceEnergyExtractor) {
    using namespace tempo::nn;
    using namespace tempo::serialization;
    using namespace tempo;
    using namespace tempo::testing;
    const auto tasks = tempo::testing::createTasks({{6, 6}, {2, 4}, {1, 5}});
    ProblemInstance instance(tasks, {{1, tasks, {1, 1, 1}}}, {}, {{3, 0}, {3, 0}, 0, 0});
    DummyScheduler scheduler(std::vector<int>{-1, -1, -1, 12}, std::vector<int>{-1, -1, -1, 12});
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
    const auto tasks = tempo::testing::createTasks({{4, 4}, {2, 4}, {8, 10}});
    ProblemInstance instance(tasks, {{4, tasks, {3, 2, 2}}}, {}, {{3, 0}, {3, 0}, 0, 0});
    DummyScheduler scheduler(std::vector<int>{-1, -1, -1, 9}, std::vector<int>{-1, -1, -1, 9});
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
    DummyScheduler scheduler(std::vector<int>{4}, std::vector<int>{4});
    Topology topology;
    EdgeVector edges{{0, 1}, {2, 1}};
    topology.edgeIndices = util::makeIndexTensor(edges);
    torch::Tensor features = TimingEdgeExtractor()(topology, makeSolverState(timings, scheduler),
                                                   ProblemInstance({}, {}, {}, {{0, 0}, {0, 0}, 0, 0}));
    ASSERT_EQ(features.size(0), edges.size());
    ASSERT_EQ(features.size(1), 2);
    ASSERT_EQ(features.sizes().size(), 2);
    EXPECT_EQ(features[0][0].item<float>(), 0);
    EXPECT_EQ(features[0][1].item<float>(), 1);
    EXPECT_EQ(features[1][0].item<float>(), 1);
    EXPECT_EQ(features[1][1].item<float>(), -0.5);
}