/**
 * @author Tim Luchterhand
 * @date 21.06.23.
 */

#ifndef TEMPO_GRAPHBUILDER_HPP
#define TEMPO_GRAPHBUILDER_HPP
#include <string>
#include <variant>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <fstream>
#include <ranges>
#include <algorithm>
#include <optional>
#include <Iterators.hpp>

#include "feature_extractors.hpp"
#include "topology_extractors.hpp"
#include "torch_types.hpp"
#include "util/traits.hpp"
#include "util/parsing/format.hpp"

namespace tempo::nn {

    MAKE_FACTORY_PATTERN(FeatureExtractor, nlohmann::json, TaskTimingFeatureExtractor, TimingEdgeExtractor,
                         ResourceEnergyExtractor)

    MAKE_FACTORY_PATTERN(TopologyBuilder, ProblemInstance, MinimalTopologyBuilder)

    /**
     * Serializable type that holds the configuration of a feature extractor type
     */
    struct ExtractorConfig {
        std::string extractorName; ///< name of the feature extractor (same as the class name)
        nlohmann::json arguments; ///< arguments to the constructor of the feature extractor
    };

    /**
     * Serializable type that holds the full configuration of a GraphBuilder
     */
    struct GraphBuilderConfig {
        ExtractorConfig taskFeatureExtractor, edgeFeatureExtractor, resourceFeatureExtractor, topologyExtractor;
    };

    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ExtractorConfig, extractorName, arguments)
    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(GraphBuilderConfig, taskFeatureExtractor, edgeFeatureExtractor,
                                       resourceFeatureExtractor, topologyExtractor)

    /**
     * Class that extracts all tensors necessary for running the inference of a GNN from a given problem and temporal
     * network
     */
    class GraphBuilder {
        template<typename EvtFun>
        struct FactoryChecker {
            HOLDS_FOR_ALL(FeatureExtractor, feature_extractor, EvtFun)
            HOLDS_FOR_ALL(TopologyExtractor, topology_extractor, EvtFun)
            static constexpr bool value = __FeatureExtractor_tester__<FeatureExtractor>::value and
                                          __TopologyExtractor_tester__<TopologyBuilder>::value;
        };
    public:

        /**
         * CTor
         * @param configPath path to the feature extractor configurations
         * @param problemInstance Initial problem containing immutable information about tasks and resources
         */
        GraphBuilder(const std::filesystem::path &configPath, const ProblemInstance &problemInstance);

        /**
         * Extracts all features and topology information from a given event distance network
         * @return pair(InputGraph containing topology and feature information, vector containing all edges of the
         * problem)
         */
        template<concepts::arbitrary_event_dist_fun EvtFun>
        auto getGraph(const EvtFun& distances) -> InputGraph {
            static_assert(FactoryChecker<EvtFun>::value,
                          "Not all feature extractors or topology extractors have the correct signature");
            auto topology = std::visit(
                    [&](auto &extractor) { return extractor.getTopology(makeSolverState(distances)); },
                    *topologyExtractor);
            auto taskFeats = std::visit(
                    [t = std::cref(topology), &distances](const auto &extractor) { return extractor(t, distances); },
                    taskFeatureExtractor);
            auto edgeFeats = std::visit(
                    [t = std::cref(topology), &distances](const auto &extractor) { return extractor(t, distances); },
                    edgeFeatureExtractor);
            auto resourceFeats = std::visit(
                    [t = std::cref(topology), &distances](const auto &extractor) { return extractor(t, distances); },
                    resourceFeatureExtractor);
            InputGraph ret;
            using K = GraphKeys;
            ret.insert(K::TaskFeatures, std::move(taskFeats));
            ret.insert(K::EdgeFeatures, std::move(edgeFeats));
            ret.insert(K::ResourceFeatures, std::move(resourceFeats));
            ret.insert(K::ResourceConsumptions, std::move(topology.resourceDemands));
            ret.insert(K::EdgeIdx, std::move(topology.edgeIndices));
            ret.insert(K::ResourceDependencies, std::move(topology.resourceDependencies));
            ret.insert(K::EdgeResourceRelations, std::move(topology.edgeResourceRelations));
            ret.insert(K::EdgePairMask, std::move(topology.edgePairMask));
            return ret;
        }


    private:
        std::optional<TopologyBuilder> topologyExtractor;
        FeatureExtractor taskFeatureExtractor;
        FeatureExtractor edgeFeatureExtractor;
        FeatureExtractor resourceFeatureExtractor;
    };
}

#endif //TEMPO_GRAPHBUILDER_HPP
