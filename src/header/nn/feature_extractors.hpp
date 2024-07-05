/**
 * @brief File containing different methods of extracting graph features
 * @author Tim Luchterhnad
 */

#ifndef TEMPO_FEATURE_EXTRACTORS_HPP
#define TEMPO_FEATURE_EXTRACTORS_HPP
#include <unordered_map>
#include <torch/torch.h>
#include <variant>
#include <vector>
#include <span>
#include <Iterators.hpp>
#include <nlohmann/json.hpp>

#include "util/traits.hpp"
#include "util/serialization.hpp"
#include "util/SchedulingProblemHelper.hpp"
#include "tensor_utils.hpp"
#include "torch_types.hpp"
#include "util/factory_pattern.hpp"

/**
 * @brief namespace containing neural network specific code
 */
namespace tempo::nn {
    /**
     * @brief Requirement for a type that can be used as resource or task feature extractor
     * @details @copybrief
     * Requires a call operator that takes a graph topology and an event network and returns a tensor containing all
     * features
     * @tparam Extractor
     */
    template<typename Extractor, typename EvtFun>
    concept feature_extractor = concepts::callable_r<Extractor, torch::Tensor, const Topology, const EvtFun>
            and concepts::arbitrary_event_dist_fun<EvtFun>;

    /**
     * @brief Extracts features from timing information of tasks.
     * @details @copybrief
     * These are the normalized shortest and longest execution time as well as the normalized minimal offset from
     * the origin and horizon.
     */
    struct TaskTimingFeatureExtractor {
        /**
         * Returns a tensor containing timing features for all tasks.
         * @param topology Graph topology for which to extract the features
         * @param eventDistances current event network from which the features are generated
         * @return torch::Tensor containing task features
         */
        template<concepts::arbitrary_event_dist_fun EvtFun>
        auto operator()(const Topology &topology, const EvtFun &eventDistances) const -> torch::Tensor {
            using namespace torch::indexing;
            const auto ub = upperBound(eventDistances);
            torch::Tensor ret = torch::empty({static_cast<long>(topology.numTasks), 4}, dataTensorOptions());
            for (long t = 0; t < topology.numTasks; ++t) {
                util::sliceAssign(ret, {static_cast<long>(t), util::SliceHere},
                                  {minDuration(t, eventDistances)/ static_cast<DataType>(ub),
                                   maxDuration(t, eventDistances) / static_cast<DataType>(ub),
                                   earliestStartTime(t, eventDistances) / static_cast<DataType>(ub),
                                   latestCompletion(t, eventDistances) / static_cast<DataType>(ub)});
            }

            return ret;
        }
    };

    /**
     * @brief Calculates the scalar resource energy feature.
     * @details @copybrief
     * \f$e_r = \sum_{t \in T} t.t_\text{avg} \cdot t.consumption_r \f$
     *
     * where \f$ t.t_\text{avg} \f$ is the normalized average duration of task \f$t\f$ and
     * \f$ t.consumption_r \f$ is the normalized consumption of resource \f$r\f$ by \f$t\f$
     */
    struct ResourceEnergyExtractor {
        /**
         * Returns a tensor containing resource energy feature for all resources.
         * @param topology Graph topology for which to extract the features
         * @param eventDistances current event network from which the features are generated
         * @return torch::Tensor containing resource features
         */
        template<concepts::arbitrary_event_dist_fun EvtFun>
        auto operator()(const Topology &topology, const EvtFun &eventDistances) const -> torch::Tensor {
            using namespace torch::indexing;
            const auto ub = upperBound(eventDistances);
            auto ret = torch::zeros({static_cast<long>(topology.numResources), 1}, dataTensorOptions());
            const std::span demands(topology.resourceDemands.data_ptr<DataType>(), topology.resourceDemands.size(0));
            const auto tasks = util::getIndexSlice(topology.resourceDependencies, 0);
            const auto resources = util::getIndexSlice(topology.resourceDependencies, 1);
            std::span energy(ret.data_ptr<DataType>(), topology.numResources);
            for (auto [task, res, demand] : iterators::zip(tasks, resources, demands)) {
                auto avgDuration = static_cast<DataType>(minDuration(task, eventDistances)
                        + maxDuration(task, eventDistances)) / DataType(2) / ub;
                energy[res] += avgDuration * demand;
            }

            return ret;
        }
    };

    /**
     * @brief Calculates timing features for edges
     * @details @copybrief
     * The edge feature contains first a flag that indicates whether the edge is an already fixed precedence edge or not
     * and second the normalized temporal difference between source and destination task
     */
    struct TimingEdgeExtractor {
        /**
         * Returns a tensor containing timing features for all edges.
         * @param topology Graph topology for which to extract the features
         * @param eventDistances current event network from which the features are generated
         * @return torch::Tensor containing edge features
         */
        template<concepts::arbitrary_event_dist_fun EvtFun>
        auto operator()(const Topology &topology, const EvtFun &eventDistances) const -> torch::Tensor {
            using namespace torch::indexing;
            const auto ub = upperBound(eventDistances);
            const auto edgeFrom = util::getIndexSlice(topology.edgeIndices, 0);
            const auto edgeTo = util::getIndexSlice(topology.edgeIndices, 1);
            auto ret = torch::empty({static_cast<long>(edgeFrom.size()), 2}, dataTensorOptions());
            for (auto [idx, from, to] : iterators::zip_enumerate(edgeFrom, edgeTo, 0l)) {
                const auto tempDiff = taskDistance(from, to, eventDistances);
                const auto tempDiffRev = taskDistance(to, from, eventDistances);
                const bool isPrecedence = tempDiff <= 0 or tempDiffRev <= 0;
                util::sliceAssign(ret, {idx, util::SliceHere}, {static_cast<DataType>(isPrecedence),
                                                                tempDiff / static_cast<DataType>(ub)});
            }

            return ret;
        }
    };

    MAKE_DEFAULT_FACTORY(TaskTimingFeatureExtractor, const nlohmann::json&)
    MAKE_DEFAULT_FACTORY(ResourceEnergyExtractor, const nlohmann::json&)
    MAKE_DEFAULT_FACTORY(TimingEdgeExtractor, const nlohmann::json&)
}


#endif //TEMPO_FEATURE_EXTRACTORS_HPP
