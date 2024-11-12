/**
* @author Tim Luchterhand
* @date 06.11.24
* @brief Contains a helper class that can be used to implement GNN based relaxation policies
*/

#ifndef TEMPO_GNNRELAXATIONPREDICTOR_HPP
#define TEMPO_GNNRELAXATIONPREDICTOR_HPP

#include <vector>
#include <filesystem>
#include <Iterators.hpp>

#include "GNN.hpp"
#include "Literal.hpp"
#include "util/traits.hpp"
#include "util/Matrix.hpp"
#include "GraphBuilder.hpp"
#include "util/SchedulingProblemHelper.hpp"
#include "heuristics/LiteralPredictor.hpp"
#include "Solver.hpp"

namespace tempo::nn {
    namespace fs = std::filesystem;

    template<typename P, typename T>
    concept edge_provider = requires(const P provider, Literal<T> lit) {
        { provider.getEdge(lit) } -> std::same_as<DistanceConstraint<T>>;
    };

    template<typename P>
    concept boolean_value_provider = requires(const P provider, var_t x) {
        { provider.value(x) } -> std::same_as<bool>;
    };

    /**
     * @brief GNN that predicts which literals to relax
     * @tparam T timing type
     * @tparam R resource type
     */
    template<concepts::scalar T, SchedulingResource R>
    class GNNRelaxationPredictor: public heuristics::LiteralPredictor<GNNRelaxationPredictor<T, R>, T> {
        EdgeRegressor gnn;
        GraphBuilder<T, R> graphBuilder;
        std::vector<DataType> literalConfidences;
    public:

        /**
         * Ctor
         * @param modelLocation path to gnn location
         * @param featureExtractorConfigLocation path to feature extraction config
         * @param problemInstance model of the problem
         * @param literals list with all search literals
         */
        GNNRelaxationPredictor(const fs::path &modelLocation, const fs::path &featureExtractorConfigLocation,
                               SchedulingProblemHelper<T, R> problemInstance, std::vector<Literal<T>> literals) :
                heuristics::LiteralPredictor<GNNRelaxationPredictor<T, R>, T>(std::move(literals)),
                gnn(modelLocation), graphBuilder(featureExtractorConfigLocation, std::move(problemInstance)),
                literalConfidences(this->numLiterals(), 0) {
            for (auto l: this->literals) {
                if (not l.hasSemantic() or not l.isBoolean()) {
                    throw std::runtime_error("all literals need to be boolean and have a semantic");
                }
            }
        }

        /**
         * Update confidences using the current state of the solver
         * @param solver Solver that represents the current state of the search
         * @throws std::runtime_error if no solution has been found yet
         */
        template<boolean_value_provider BVP>
        void updateConfidence(const Solver<T> &solver, const BVP &boolValProvider) {
            const auto &boolean = solver.boolean;
            if (not boolean.hasSolution()) {
                throw std::runtime_error("solver needs to have found at least one solution");
            }

            const auto graph = graphBuilder.getGraph(
                    makeSolverState(graphBuilder.getProblem().getTaskDistances(solver), solver));
            const auto edgeHeatMap = gnn.getHeatMap(graph);
            const auto &mapping = graphBuilder.getProblem().getMapping();
            for (auto [lit, confidence]: iterators::zip(this->literals, literalConfidences)) {
                const auto &edge = boolean.getEdge(lit);
                const auto tFrom = mapping(edge.from);
                const auto tTo = mapping(edge.to);
                const auto forward = probabilityMass(tFrom, tTo, edgeHeatMap);
                const auto backward = probabilityMass(tTo, tFrom, edgeHeatMap);
                confidence = (forward + backward) / 2;
                BooleanVar var(lit);
                lit = var == boolValProvider.value(var.id());
            }
        }

        /**
         * Sets all confidence values to 0
         */
        void resetConfidences() noexcept {
            literalConfidences.assign(this->literals.size(), 0);
        }

        /**
         * resets confidence values to 0 and updates possible literals
         * @param literals all valid search literals
         */
        void reinitialize(std::vector<Literal<T>> literals) noexcept {
            this->literals = std::move(literals);
            resetConfidences();
        }

        /**
         * @copydoc reinitialize
         * overload that infers search literals from solver
         */
        void reinitialize(const Solver<T> &solver) {
            reinitialize(graphBuilder.getProblem().getSearchLiterals(solver));
        }

        DataType getCertainty(Literal<T>, std::size_t idx) const noexcept {
            return literalConfidences[idx];
        }

        Literal<T> getPolarity(Literal<T> lit, std::size_t) const noexcept {
            return lit;
        }
    };
}


#endif //TEMPO_GNNRELAXATIONPREDICTOR_HPP
