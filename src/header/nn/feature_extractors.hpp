/**
 * @brief File containing different methods of extracting graph features
 * @author Tim Luchterhnad
 */

#ifndef TEMPO_FEATURE_EXTRACTORS_HPP
#define TEMPO_FEATURE_EXTRACTORS_HPP
#include <torch/torch.h>
#include <variant>
#include <span>
#include <optional>
#include <Iterators.hpp>
#include <nlohmann/json.hpp>

#include "util/traits.hpp"
#include "util/edge_distance.hpp"
#include "util/SchedulingProblemHelper.hpp"
#include "tensor_utils.hpp"
#include "torch_types.hpp"
#include "util/factory_pattern.hpp"
#include "util/Matrix.hpp"

/**
 * @brief namespace containing neural network specific code
 */
namespace tempo::nn {
    /**
     * @brief Extracts features from timing information of tasks.
     * @details @copybrief
     * These are the normalized shortest and longest execution time as well as the normalized minimal offset from
     * the origin and horizon.
     */
    struct TaskTimingFeatureExtractor {
        static constexpr auto LegacyKey = "legacy";
        bool legacyFeatures = false;

        /**
         * Returns a tensor containing timing features for all tasks.
         * @tparam DistFun type of distance function
         * @tparam S Solver type that has a member that provides upper and lower bounds
         * @tparam T timing type
         * @tparam R scheduling resource type
         * @param topology Graph topology for which to extract the features
         * @param state current state of the solver
         * @param problem scheduling problem description
         * @return torch::Tensor containing task features
         */
        template<typename DistFun, distance_provider S, concepts::scalar T, SchedulingResource R>
        auto operator()(const Topology &topology, const SolverState<DistFun, S> &state,
                        const SchedulingProblemHelper<T, R> &problem) const -> torch::Tensor {
            const auto ub = static_cast<DataType>(problem.schedule().getLatestEnd(state.solver));
            torch::Tensor ret = torch::empty({static_cast<long>(topology.numTasks), 4}, dataTensorOptions());
            for (auto [t, task]: iterators::enumerate(problem.tasks(), 0l)) {
                util::sliceAssign(ret, {t, util::SliceHere},
                                  {task.minDuration(state.solver) / ub,
                                   task.maxDuration(state.solver) / ub,
                                   task.getEarliestStart(state.solver) / ub *
                                   static_cast<DataType>(1 - 2 * legacyFeatures),
                                   task.getLatestEnd(state.solver) / ub - static_cast<DataType>(legacyFeatures)});
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
         * @tparam DistFun type of distance function
         * @tparam S Solver type that has a member that provides upper and lower bounds
         * @tparam T timing type
         * @tparam R scheduling resource type
         * @param topology Graph topology for which to extract the features
         * @param state current state of the solver
         * @param problem scheduling problem description
         * @return torch::Tensor containing resource features
         */
        template<typename Dist, distance_provider S, concepts::scalar T, SchedulingResource R>
        auto operator()(const Topology &topology, const SolverState<Dist, S> &state,
                        const SchedulingProblemHelper<T, R> &problem) const -> torch::Tensor {
            const auto ub = static_cast<DataType>(problem.schedule().getLatestEnd(state.solver));
            auto ret = torch::zeros({static_cast<long>(topology.numResources), 1}, dataTensorOptions());
            const std::span demands(topology.resourceDemands.data_ptr<DataType>(), topology.resourceDemands.size(0));
            const auto tasks = util::getIndexSlice(topology.resourceDependencies, 0);
            const auto resources = util::getIndexSlice(topology.resourceDependencies, 1);
            std::span energy(ret.data_ptr<DataType>(), topology.numResources);
            for (auto [task, res, demand] : iterators::zip(tasks, resources, demands)) {
                auto avgDuration = static_cast<DataType>(problem.tasks()[task].minDuration(state.solver)
                        + problem.tasks()[task].maxDuration(state.solver)) / DataType(2) / ub;
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
    class TimingEdgeExtractor {
        std::optional<Matrix<bool>> taskPrecedences;
    public:
        static constexpr auto LegacyPrecedences = "precedencesFromDistances";

        /**
         * Ctor
         * @param taskPrecedences bool matrix indicating precedences between tasks
         */
        explicit TimingEdgeExtractor(std::optional<Matrix<bool>> taskPrecedences) noexcept: taskPrecedences(
            std::move(taskPrecedences)) {}

        /**
         * Returns a tensor containing timing features for all edges.
         * @tparam DistFun type of distance function
         * @tparam S Solver type that has a member that provides upper and lower bounds
         * @tparam T timing type
         * @tparam R scheduling resource type
         * @param topology Graph topology for which to extract the features
         * @param state current state of the solver
         * @param problem scheduling problem description
         * @return torch::Tensor containing edge features
         */
        template<concepts::arbitrary_task_dist_fun Dist, typename S, concepts::scalar T, SchedulingResource R>
        auto operator()(const Topology &topology, const SolverState<Dist, S> &state,
                        const SchedulingProblemHelper<T, R> &problem) const -> torch::Tensor {
            const auto ub = static_cast<DataType>(problem.schedule().getLatestEnd(state.solver));
            const auto edgeFrom = util::getIndexSlice(topology.edgeIndices, 0);
            const auto edgeTo = util::getIndexSlice(topology.edgeIndices, 1);
            auto ret = torch::empty({static_cast<long>(edgeFrom.size()), 2}, dataTensorOptions());
            for (auto [idx, from, to] : iterators::zip_enumerate(edgeFrom, edgeTo, 0l)) {
                const auto tempDiff = state.distance(from, to);
                bool isPrecedence;
                if (taskPrecedences.has_value()) {
                    isPrecedence = taskPrecedences->at(from, to);
                } else {
                    const auto tempDiffRev = state.distance(to, from);
                    isPrecedence = tempDiff <= 0 or tempDiffRev <= 0;
                }

                util::sliceAssign(ret, {idx, util::SliceHere}, {static_cast<DataType>(isPrecedence),
                                                                static_cast<DataType>(tempDiff) / ub});
            }

            return ret;
        }
    };

    MAKE_TEMPLATE_FACTORY(TaskTimingFeatureExtractor, ESCAPE(concepts::scalar T, SchedulingResource R),
                          const nlohmann::json &config, const SchedulingProblemHelper<T, R> &) {
            return TaskTimingFeatureExtractor{.legacyFeatures = config.at(
                    TaskTimingFeatureExtractor::LegacyKey).get<bool>()};
        }
    };

    MAKE_DEFAULT_TEMPLATE_FACTORY(ResourceEnergyExtractor, ESCAPE(concepts::scalar T, SchedulingResource R),
                                  const nlohmann::json&, const SchedulingProblemHelper<T, R> &)

    MAKE_TEMPLATE_FACTORY(TimingEdgeExtractor, ESCAPE(concepts::scalar T, SchedulingResource R),
                          const nlohmann::json &config, const SchedulingProblemHelper<T, R> &problem) {
            if (config.at(TimingEdgeExtractor::LegacyPrecedences).get<bool>()) {
                return TimingEdgeExtractor(std::nullopt);
            }

            Matrix<bool> taskPrecedences(problem.tasks().size(), problem.tasks().size(), false);
            const auto &mapping = problem.getMapping();
            for (const auto &prec : problem.precedences()) {
                if (problem.hasVariable(prec.from) and problem.hasVariable(prec.to)) {
                    taskPrecedences(mapping(prec.from), mapping(prec.to)) = true;
                    taskPrecedences(mapping(prec.to), mapping(prec.from)) = true;
                }
            }

            return TimingEdgeExtractor(std::move(taskPrecedences));
        }
    };
}


#endif //TEMPO_FEATURE_EXTRACTORS_HPP
