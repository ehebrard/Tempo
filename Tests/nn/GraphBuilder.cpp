/**
 * @author Tim Luchterhand
 * @date 27.06.23.
 */

#include <gtest/gtest.h>
#include <ranges>

#include "nn/GraphBuilder.hpp"
#include "testing.hpp"


auto getGtTopology(const ProblemInstance &problemInstance) -> tempo::nn::Topology {
    using namespace tempo::nn;
    using namespace tempo::testing;
    std::ifstream configFile(TestData::GraphBuilderConfig);
    auto config = nlohmann::json::parse(configFile).get<GraphBuilderConfig>();
    auto extractor = TopologyBuilderFactory::getInstance().create(config.topologyExtractor.extractorName,
                                                                  problemInstance);
    return std::visit([](auto &e) {
        return e.getTopology(makeSolverState(tempo::Matrix<int>{}));
    }, extractor);
}

void testFeatureExtraction(const std::string &featureName, const std::string &graphKey) {
    using namespace tempo::nn;
    using namespace tempo::testing;
    const auto [problem, evtNetwork] = createRandomProblem(100, 50);
    const auto topology = getGtTopology(problem);
    GraphBuilder graphBuilder(TestData::GraphBuilderConfig, problem);
    auto graph = graphBuilder.getGraph(makeSolverState(evtNetwork));
    std::ifstream configFile(TestData::GraphBuilderConfig);
    auto config = nlohmann::json::parse(configFile).at(featureName).get<ExtractorConfig>();
    auto extractor = FeatureExtractorFactory::getInstance().create(config.extractorName, config.arguments);
    auto gt = std::visit([&evtNetwork, &topology](auto &ext) { return ext(topology, evtNetwork); }, extractor);
    EXPECT_TRUE(torch::all(gt == graph.at(graphKey)).item<bool>());
}

TEST(nn_GraphBuilder, get_graph_basic) {
    using namespace tempo::nn;
    using namespace tempo::testing;
    const auto [problem, evtNetwork] = createRandomProblem(3, 2);
    GraphBuilder graphBuilder(TestData::GraphBuilderConfig, problem);
    auto graph = graphBuilder.getGraph(makeSolverState(evtNetwork));
    EXPECT_TRUE(graph.contains(GraphKeys::TaskFeatures));
    EXPECT_TRUE(graph.contains(GraphKeys::EdgeFeatures));
    EXPECT_TRUE(graph.contains(GraphKeys::ResourceFeatures));
    EXPECT_TRUE(graph.contains(GraphKeys::EdgeIdx));
    EXPECT_TRUE(graph.contains(GraphKeys::EdgePairMask));
    EXPECT_TRUE(graph.contains(GraphKeys::EdgeResourceRelations));
    EXPECT_TRUE(graph.contains(GraphKeys::ResourceDependencies));
    EXPECT_TRUE(graph.contains(GraphKeys::ResourceConsumptions));
}


TEST(nn_GraphBuilder, get_graph_task_features) {
    using namespace tempo::nn;
    testFeatureExtraction("taskFeatureExtractor", tempo::nn::GraphKeys::TaskFeatures);
}

TEST(nn_GraphBuilder, get_graph_edge_features) {
    using namespace tempo::nn;
    testFeatureExtraction("edgeFeatureExtractor", tempo::nn::GraphKeys::EdgeFeatures);
}

TEST(nn_GraphBuilder, get_graph_resource_features) {
    testFeatureExtraction("resourceFeatureExtractor", tempo::nn::GraphKeys::ResourceFeatures);
}

TEST(nn_GraphBuilder, get_graph_resource_deps) {
    using namespace tempo::nn;
    using namespace tempo::testing;
    auto [problem, evtNetwork] = createRandomProblem(10, 5);
    GraphBuilder graphBuilder(TestData::GraphBuilderConfig, problem);
    auto graph = graphBuilder.getGraph(makeSolverState(evtNetwork));
    const auto topology = getGtTopology(problem);
    EXPECT_TRUE(torch::all(topology.resourceDependencies == graph.at(GraphKeys::ResourceDependencies)).item<bool>());
    EXPECT_TRUE(torch::all(topology.resourceDemands == graph.at(GraphKeys::ResourceConsumptions)).item<bool>());
}

TEST(nn_GraphBuilder, get_graph_topology) {
    using namespace tempo::nn;
    using namespace tempo::testing;
    const auto [problem, evtNetwork] = createRandomProblem(10, 5);
    GraphBuilder graphBuilder(TestData::GraphBuilderConfig, problem);
    const auto graph = graphBuilder.getGraph(makeSolverState(evtNetwork));
    const auto topology = getGtTopology(problem);
    EXPECT_TRUE(torch::all(topology.edgePairMask == graph.at(GraphKeys::EdgePairMask)).item<bool>());
    EXPECT_TRUE(torch::all(topology.edgeResourceRelations == graph.at(GraphKeys::EdgeResourceRelations)).item<bool>());
    EXPECT_TRUE(torch::all(topology.edgeIndices == graph.at(GraphKeys::EdgeIdx)).item<bool>());
}
