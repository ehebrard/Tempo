/**
* @author Tim Luchterhand
* @date 11.09.24
* @brief
*/

#ifndef TEMPO_GNNPRECEDENCEPREDICTOR_HPP
#define TEMPO_GNNPRECEDENCEPREDICTOR_HPP

#include <vector>
#include <filesystem>
#include <Iterators.hpp>

#include "GNN.hpp"
#include "GraphBuilder.hpp"
#include "Literal.hpp"
#include "util/SchedulingProblemHelper.hpp"
#include "util/traits.hpp"
#include "torch_types.hpp"
#include "heat_map_utils.hpp"
#include "heuristics/PrecedencePredictor.hpp"

namespace tempo {
    template<typename T>
    class Solver;
}

namespace tempo::nn {
    namespace fs = std::filesystem;

    /**
     * @brief GNN based precedence prediction heuristic
     * @detail @copybrief
     * @tparam T timing type
     * @tparam R resource type
     */
    template<concepts::scalar T, SchedulingResource R>
    class GNNPrecedencePredictor : public heuristics::PrecedencePredictor<GNNPrecedencePredictor<T, R>, T> {
        EdgeRegressor gnn;
        GraphBuilder<T, R> graphBuilder;
    public:

        /**
         * Ctor
         * @param modelLocation path to gnn location
         * @param featureExtractorConfigLocation path to feature extraction config
         * @param problemInstance model of the problem
         * @param literals list with all search literals
         */
        GNNPrecedencePredictor(const fs::path &modelLocation, const fs::path &featureExtractorConfigLocation,
                               SchedulingProblemHelper<T, R> problemInstance, std::vector<Literal<T>> literals) :
                heuristics::PrecedencePredictor<GNNPrecedencePredictor<T, R>, T>(std::move(literals)), gnn(
                modelLocation), graphBuilder(featureExtractorConfigLocation, std::move(problemInstance)) {
            for (auto l: this->literals) {
                if (not l.hasSemantic() or not l.isBoolean()) {
                    throw std::runtime_error("all literals need to be boolean and have a semantic");
                }
            }
        }

        /**
         * Update confidences using the current state of the solver
         * @param solver Solver that represents the current state of the search
         */
        void updateConfidence(const Solver<T> &solver) {
            const auto taskNetwork = graphBuilder.getProblem().getTaskDistances(solver);
            const auto graph = graphBuilder.getGraph(makeSolverState(std::move(taskNetwork), solver));
            const auto edgeHeatMap = gnn.getHeatMap(graph);
            const auto &mapping = graphBuilder.getProblem().getMapping();
            for (auto [lit, pos, neg] : iterators::zip(this->literals, this->massesPos, this->massesNeg)) {
                pos += probabilityMass(lit, edgeHeatMap, solver, mapping);
                neg += probabilityMass(~lit, edgeHeatMap, solver, mapping);
            }
        }

        /**
         * Sets all confidence values to 0
         */
        void resetConfidences() noexcept {
            for (auto [pos, neg] : iterators::zip(this->massesPos, this->massesNeg)) {
                pos = neg = 0;
            }
        }

        /**
         * resets confidence values to 0 and updates possible literals
         * @param literals all valid search literals
         */
        void reinitialize(std::vector<Literal<T>> literals) {
            this->literals = std::move(literals);
            this->massesPos = std::vector(this->literals.size(), 0.0);
            this->massesNeg = std::vector(this->literals.size(), 0.0);
        }

        /**
         * @copydoc reinitialize
         * overload that infers search literals from solver
         */
        void reinitialize(const Solver<T> &solver) {
            reinitialize(graphBuilder.getProblem().getSearchLiterals(solver));
        }

        /**
         * Calculates prediction certainty from two evidence masses (not really theoretically motivated)
         * @param mPos positive evidence for a literal
         * @param mNeg negative evidence for a literal
         * @return confidence value in [0, 1]
         */
        static double getCertainty(DataType mPos, DataType mNeg) {
            // = |m1 - m2| / (m1 + m2)
            return std::abs(2.0 * mPos / (mPos + mNeg) - 1);
        }
    };
}

#endif //TEMPO_GNNPRECEDENCEPREDICTOR_HPP
