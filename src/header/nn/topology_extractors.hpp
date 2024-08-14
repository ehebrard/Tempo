/**
 * @author Tim Luchterhand
 * @date 17.11.23.
 */

#ifndef TEMPO_TOPOLOGY_EXTRACTORS_HPP
#define TEMPO_TOPOLOGY_EXTRACTORS_HPP

#include <vector>
#include <cassert>

#include "Global.hpp"
#include "Model.hpp"
#include "util/SchedulingProblemHelper.hpp"
#include "torch_types.hpp"
#include "tensor_utils.hpp"
#include "util/Matrix.hpp"
#include "util/factory_pattern.hpp"
#include "util/traits.hpp"
#include "DistanceConstraint.hpp"

namespace tempo::nn {

    namespace impl {
        struct EdgeLookup : protected Matrix<IndexType> {
            static constexpr IndexType NoValue = std::numeric_limits<IndexType>::lowest();

            constexpr explicit EdgeLookup(std::size_t numEvents): Matrix<IndexType>(numEvents, numEvents, NoValue) {}

            [[nodiscard]] constexpr bool contains(const Edge &e) const noexcept {
                return static_cast<const Matrix<IndexType>&>(*this)(e.first, e.second) != NoValue;
            }

            [[nodiscard]] constexpr bool containsSafe(const Edge &e) const {
                return this->at(e.first, e.second) != NoValue;
            }

            constexpr auto operator()(const Edge &e) noexcept -> Matrix<IndexType>::reference {
                return static_cast<Matrix<IndexType>&>(*this)(e.first, e.second);
            }

            constexpr auto operator()(const Edge &e) const noexcept -> Matrix<IndexType>::const_reference {
                return static_cast<const Matrix<IndexType>&>(*this)(e.first, e.second);
            }

            [[nodiscard]] std::size_t getNumEdges() const;

            using Matrix<IndexType>::at;
            using Matrix<IndexType>::for_each;
        };

        struct TopologyData {
            EdgeLookup edgeLookup;
            std::vector<IndexType> taskIdx{};
            std::vector<IndexType> resIdx{};
            std::vector<DataType> resDemands{};
            std::vector<IndexType> edgeIdx{};
            std::vector<IndexType> edgeRelResIdx{};
            EdgeVector edges{};
            std::vector<IndexType> edgePairMask{};
            IndexType pairMaskVal{};
        };
    }

    /**
     * Topology extractor that creates a minimum redundancy topology from a set of decided precedences
     */
    struct MinimalTopologyBuilder {
        MinimalTopologyBuilder(const MinimalTopologyBuilder&) = delete;
        MinimalTopologyBuilder &operator=(const MinimalTopologyBuilder&) = delete;
        MinimalTopologyBuilder(MinimalTopologyBuilder&&) = default;
        MinimalTopologyBuilder &operator=(MinimalTopologyBuilder&&) = default;
        ~MinimalTopologyBuilder() = default;

        /**
         * Ctor
         * @param problem initial problem instance
         * @param bidirectionalPrecedences Whether to insert bidirectional precedence edges
         */
        template<concepts::scalar T, SchedulingResource R>
        explicit MinimalTopologyBuilder(const SchedulingProblemHelper<T, R> &problem, bool bidirectionalPrecedences) {
            using namespace iterators;
            cache.numTasks = problem.tasks().size();
            cache.numResources = problem.resources().size();
            impl::TopologyData topologyData{.edgeLookup = impl::EdgeLookup(cache.numTasks * 2 + 2)};
            for (auto [r, resourceSpec]: enumerate(problem.resources(), 0l)) {
                completeSubGraph(resourceSpec, r, topologyData);
            }

            addPrecedenceEdges(problem.precedences(), problem.getMapping(), topologyData, bidirectionalPrecedences);
            cache.edgeIndices = util::makeIndexTensor(topologyData.edges);
            cache.edgePairMask = torch::from_blob(topologyData.edgePairMask.data(),
                                                  static_cast<long>(topologyData.edgePairMask.size()),
                                                  indexTensorOptions()).clone();
            cache.edgeResourceRelations = util::makeIndexTensor(topologyData.edgeIdx, topologyData.edgeRelResIdx);
            cache.resourceDependencies = util::makeIndexTensor(topologyData.taskIdx, topologyData.resIdx);
            cache.resourceDemands = torch::from_blob(topologyData.resDemands.data(),
                                                     {static_cast<long>(topologyData.taskIdx.size()), 1},
                                                     dataTensorOptions()).clone();
        }

        /**
         * Generates a graph topology based on the initial problem instance and on the given decided precedence
         * constraints
         * @param precedences range of decided precedences
         * @return the graph topology representing the initial problem and the added precedences
         */
        template<typename ...Args>
        [[nodiscard]] auto getTopology(const SolverState<Args...> &) const -> const Topology &{
            return cache;
        }

        /**
         * Generates a graph topology based only on the initial problem instance
         * @return the graph topology representing the initial problem
         */
        [[nodiscard]] auto getTopology() const -> const Topology&;

    protected:
        template<concepts::ttyped_range<DistanceConstraint> R>
        static void addPrecedenceEdges(const R &precedences, const VarTaskMapping &vtMapping,
                                       impl::TopologyData &topologyData, bool bidirectional) {
            for (const auto &p : precedences) {
                if(not vtMapping.contains(p.from) or not vtMapping.contains(p.to)) {
                    continue;
                }

                Edge e(vtMapping(p.from), vtMapping(p.to));
                if (e.first == e.second or topologyData.edgeLookup.contains(e)) {
                    continue;
                }

                addEdge(e, false, topologyData.pairMaskVal, topologyData);
                if (bidirectional) {
                    Edge rev(e.second, e.first);
                    addEdge(rev, false, topologyData.pairMaskVal, topologyData);
                }

                ++topologyData.pairMaskVal;
            }
        }

        static void addEdge(const Edge &e, bool isResourceEdge, IndexType maskVal, impl::TopologyData &topologyData);

        template<SchedulingResource R>
        static void completeSubGraph(const R &resourceSpec, IndexType resource, impl::TopologyData &topologyData) {
            const auto &taskIds = resourceSpec;
            for (std::size_t i = 0; i < taskIds.size(); ++i) {
                topologyData.taskIdx.emplace_back(taskIds[i]);
                topologyData.resIdx.emplace_back(resource);
                topologyData.resDemands.emplace_back(static_cast<DataType>(resourceSpec.getDemand(i)) /
                                                     static_cast<DataType>(resourceSpec.resourceCapacity()));
                for (std::size_t j = i + 1; j < taskIds.size(); ++j) {
                    Edge e(taskIds[i], taskIds[j]);
                    Edge rev(taskIds[j], taskIds[i]);
                    topologyData.edgeRelResIdx.emplace_back(resource);
                    topologyData.edgeRelResIdx.emplace_back(resource);
                    if (topologyData.edgeLookup.contains(e)) {
                        topologyData.edgeIdx.emplace_back(topologyData.edgeLookup(e));
                        topologyData.edgeIdx.emplace_back(topologyData.edgeLookup(rev));
                        continue;
                    }

                    addEdge(e, true, topologyData.pairMaskVal, topologyData);
                    addEdge(rev, true, topologyData.pairMaskVal, topologyData);
                    ++topologyData.pairMaskVal;
                }
            }
        }

    private:
        Topology cache;
    };

    MAKE_TEMPLATE_FACTORY(MinimalTopologyBuilder, ESCAPE(concepts::scalar T, SchedulingResource R),
                          const ESCAPE(SchedulingProblemHelper<T, R> &problemInstance, const nlohmann::json &params)) {
            auto bidir = params.at("bidirectionalPrecedences").get<bool>();
            return MinimalTopologyBuilder(problemInstance, bidir);
        }
    };
}

#endif //TEMPO_TOPOLOGY_EXTRACTORS_HPP
