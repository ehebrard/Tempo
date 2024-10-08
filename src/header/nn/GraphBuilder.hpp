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
#include "util/SchedulingProblemHelper.hpp"
#include "util/serialization.hpp"

namespace tempo::nn {

    MAKE_POLYMORPHIC_TYPE(FeatureExtractor, TaskTimingFeatureExtractor, TimingEdgeExtractor, ResourceEnergyExtractor,
                          CumulativeTimingEdgeExtractor)
    MAKE_FACTORY_PATTERN(FeatureExtractor, const nlohmann::json&, TaskTimingFeatureExtractor, TimingEdgeExtractor,
                         ResourceEnergyExtractor, CumulativeTimingEdgeExtractor)

    MAKE_POLYMORPHIC_TYPE(TopologyBuilder, MinimalTopologyBuilder)
    MAKE_T_FACTORY_PATTERN_RAW(TopologyBuilder, ESCAPE(template<typename T, typename R>),
                               ESCAPE(const SchedulingProblemHelper<T, R> &p, const nlohmann::json &params),
                               ESCAPE(p, params), MinimalTopologyBuilder)

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
    template<concepts::scalar T, SchedulingResource R>
    class GraphBuilder {
    public:

        /**
         * CTor
         * @param configPath path to the feature extractor configurations
         * @param problemInstance Initial problem containing immutable information about tasks and resources
         */
        GraphBuilder(const std::filesystem::path &configPath, SchedulingProblemHelper<T, R> problemInstance)
                : problemDefinition(std::move(problemInstance)) {
            auto data = serialization::deserializeFromFile<GraphBuilderConfig>(configPath);
            topologyExtractor = TopologyBuilderFactory::getInstance().create(data.topologyExtractor.extractorName,
                                                                             problemDefinition,
                                                                             data.topologyExtractor.arguments);
            taskFeatureExtractor = FeatureExtractorFactory::getInstance().create(
                    data.taskFeatureExtractor.extractorName, data.taskFeatureExtractor.arguments);
            resourceFeatureExtractor = FeatureExtractorFactory::getInstance().create(
                    data.resourceFeatureExtractor.extractorName, data.resourceFeatureExtractor.arguments);
            edgeFeatureExtractor = FeatureExtractorFactory::getInstance().create(
                    data.edgeFeatureExtractor.extractorName, data.edgeFeatureExtractor.arguments);
        }

        /**
         * Extracts all features and topology information from a given solver state
         * @tparam Args Types of solver state fields
         * @param state current state of the solver
         * @return InputGraph containing topology and feature information
         */
        template<typename ...Args>
        auto getGraph(const SolverState<Args...> &state) -> InputGraph {
            auto topology = std::visit(
                    [&](auto &extractor) { return extractor.getTopology(state); },
                    *topologyExtractor);
            const auto &p = problemDefinition;
            auto taskFeats = std::visit([t = std::cref(topology), &p, &state](const auto &extractor) {
                return extractor(t, state, p);
            }, taskFeatureExtractor);
            auto edgeFeats = std::visit([t = std::cref(topology), &p, &state](const auto &extractor) {
                return extractor(t, state, p);
            }, edgeFeatureExtractor);
            auto resourceFeats = std::visit([t = std::cref(topology), &p, &state](const auto &extractor) {
                return extractor(t, state, p);
            }, resourceFeatureExtractor);
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

        /**
         * gets the stored problem definition
         * @return problem definition
         */
        auto getProblem() const noexcept -> const SchedulingProblemHelper<T, R>& {
            return problemDefinition;
        }

        const auto &getTaskFeatureExtractor() const noexcept {
            return taskFeatureExtractor;
        }

        const auto &getEdgeFeatureExtractor() const noexcept {
            return edgeFeatureExtractor;
        }

        const auto &getResourceFeatureExtractor() const noexcept {
            return resourceFeatureExtractor;
        }

    private:
        SchedulingProblemHelper<T, R> problemDefinition;
        std::optional<TopologyBuilder> topologyExtractor;
        FeatureExtractor taskFeatureExtractor;
        FeatureExtractor edgeFeatureExtractor;
        FeatureExtractor resourceFeatureExtractor;
    };
}

#endif //TEMPO_GRAPHBUILDER_HPP
