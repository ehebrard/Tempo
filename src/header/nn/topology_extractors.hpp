/**
 * @author Tim Luchterhand
 * @date 17.11.23.
 */

#ifndef TEMPO_TOPOLOGY_EXTRACTORS_HPP
#define TEMPO_TOPOLOGY_EXTRACTORS_HPP

#include <ranges>
#include <vector>

#include "Global.hpp"
#include "util/parsing/format.hpp"
#include "torch_types.hpp"
#include "tensor_utils.hpp"
#include "util/Matrix.hpp"
#include "util/factory_pattern.hpp"
#include "DistanceConstraint.hpp"

namespace tempo::nn {

    /**
     * Represents the state of a the solver that can be used to extract the graph topology
     * @tparam EvtFun type of event distance function
     */
    template<concepts::arbitrary_event_dist_fun EvtFun>
    struct SolverState {
        template<typename E>
        requires(concepts::arbitrary_event_dist_fun<std::remove_cvref_t<E>>)
        explicit constexpr SolverState(E &&evt) : eventNetwork(std::forward<E>(evt)) {}
        EvtFun eventNetwork; ///< event network functor object
    };

    /**
     * Concept that models the interface of a topology extractor
     * @tparam Extractor
     * @tparam EvtFun type of event function
     */
    template<typename Extractor, typename EvtFun>
    concept topology_extractor = requires(Extractor e, const SolverState<EvtFun> &state) {
        { e.getTopology(state) } -> std::convertible_to<Topology>;
    };


    /**
     * Helper function for SolverState
     * @tparam Args argument types to SolverState
     * @param args arguments to SolverState
     * @return solver state view of the arguments. Depending on the reference type of the arguments, this is only
     * a non owning object
     */
    template<typename ...Args>
    constexpr auto makeSolverState(Args &&...args) {
        return SolverState<Args...>(std::forward<Args>(args)...);
    }

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
         */
        explicit MinimalTopologyBuilder(const ProblemInstance &problem);

        /**
         * Generates a graph topology based on the initial problem instance and on the given decided precedence
         * constraints
         * @tparam EvtFun type of event distance function
         * @param precedences range of decided precedences
         * @return the graph topology representing the initial problem and the added precedences
         */
        template<typename EvtFun>
        [[nodiscard]] auto getTopology(const SolverState<EvtFun> &) const -> const Topology &{
            return cache;
        }

        /**
         * Generates a graph topology based only on the initial problem instance
         * @return the graph topology representing the initial problem
         */
        [[nodiscard]] auto getTopology() const -> const Topology&;

    protected:
        template<std::ranges::range R>
        static void addPrecedenceEdges(const R &precedences, impl::TopologyData &topologyData) {
            for (auto [from, to, dist] : precedences) {
                assert(from != ORIGIN and from != HORIZON and to != ORIGIN and to != HORIZON &&
                       "You probably passed a precedence between tasks instead of events!");
                Edge e(TASK(from), TASK(to));
                if (e.first == e.second or topologyData.edgeLookup.contains(e)) {
                    continue;
                }

                addEdge(e, false, topologyData.pairMaskVal++, topologyData);
            }
        }

        static void addEdge(const Edge &e, bool isResourceEdge, IndexType maskVal, impl::TopologyData &topologyData);

        static void completeSubGraph(const Resource<int> &resourceSpec, IndexType resource,
                                     impl::TopologyData &topologyData);

    private:
        Topology cache;
    };

    MAKE_FACTORY(MinimalTopologyBuilder, const ProblemInstance &problemInstance) {
        return MinimalTopologyBuilder(problemInstance);
    }};
}

#endif //TEMPO_TOPOLOGY_EXTRACTORS_HPP
