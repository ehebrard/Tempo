/**
 * @author Tim Luchterhand
 * @date 13.11.23.
 */

#include "nn/GraphBuilder.hpp"

namespace tempo::nn {

    GraphBuilder::GraphBuilder(const std::filesystem::path &configPath, const ProblemInstance &problemInstance) {
        auto data = serialization::deserializeFromFile<GraphBuilderConfig>(configPath);
        topologyExtractor = TopologyBuilderFactory::getInstance().create(data.topologyExtractor.extractorName,
                                                                         problemInstance);
        taskFeatureExtractor = FeatureExtractorFactory::getInstance().create(
                data.taskFeatureExtractor.extractorName, data.taskFeatureExtractor.arguments);
        resourceFeatureExtractor = FeatureExtractorFactory::getInstance().create(
                data.resourceFeatureExtractor.extractorName, data.resourceFeatureExtractor.arguments);
        edgeFeatureExtractor = FeatureExtractorFactory::getInstance().create(
                data.edgeFeatureExtractor.extractorName, data.edgeFeatureExtractor.arguments);
    }
}