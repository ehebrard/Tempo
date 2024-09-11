/**
 * @author Tim Luchterhand
 * @date 26.07.23.
 */

#ifndef TEMPO_GNNEDGEPOLARITYPREDICTOR_HPP
#define TEMPO_GNNEDGEPOLARITYPREDICTOR_HPP

#include "Global.hpp"
#include "GNN.hpp"
#include "GraphBuilder.hpp"
#include "util/Matrix.hpp"
#include "util/traits.hpp"
#include "Literal.hpp"

namespace tempo {
    template<typename T>
    class Solver;
}

namespace tempo::nn {

    /**
     * GNN based value selection heuristic that predicts a probability for each edge and returns the polarity of
     * choicepoint accordingly
     */
    template<concepts::scalar T, SchedulingResource R>
    class GNNEdgePolarityPredictor {
    public:
        /**
         * Ctor
         * @param modelLocation location of the trained GNN model
         * @param featureExtractorConfigLocation location of the feature extractor configuration
         * (tempo::nn:GraphBuilderConfig)
         * @param problemInstance initial description of the problem
         */
        GNNEdgePolarityPredictor(const std::filesystem::path &modelLocation,
                                 const std::filesystem::path &featureExtractorConfigLocation,
                                 SchedulingProblemHelper<T, R> problemInstance) :
                gnn(modelLocation),
                graphBuilder(featureExtractorConfigLocation, std::move(problemInstance)) {}

        /**
         * Returns the choice point in the appropriate polarity according to a GNN edge regressor
         * @param x variable for which to make a choice
         * @return corresponding positive or negative literal according to the GNN
         */
        auto choose(var_t x, const Solver<T> &solver) const -> Literal<T> {
            auto edge = solver.boolean.getEdge(true, x);
            auto edgeRev = solver.boolean.getEdge(false, x);
            const auto &mapping = graphBuilder.getProblem().getMapping();
            assert(mapping.contains(edge.from));
            assert(mapping.contains(edge.to));
            assert(mapping.contains(edgeRev.from));
            assert(mapping.contains(edgeRev.to));
            unsigned from = mapping(edge.from);
            unsigned to = mapping(edge.to);
            unsigned rfrom = mapping(edgeRev.from);
            unsigned rto = mapping(edgeRev.to);
            if (from != rto or to != rfrom) {
                // @TODO not strictly necessary. simply allow distances between different tasks and use DST for probs
                throw std::runtime_error("positive and negative edges have different tasks as endpoints.");
            }

            return solver.boolean.getLiteral(EdgeRegressor::dstEdgeProbability(from, to, edgeHeatMap) > 0.5, x);
        }

        /**
         * Runs the GNN to produce an edge heat map
         * @param solver scheduler instance for which to produce the heat map
         */
        void preEvaluation(const Solver<T> &solver) {
            auto taskNetwork = graphBuilder.getProblem().getTaskDistances(solver);
            auto graph = graphBuilder.getGraph(makeSolverState(std::move(taskNetwork), solver));
            edgeHeatMap = gnn.getHeatMap(graph);
        }

    private:
        EdgeRegressor gnn;
        GraphBuilder<T, R> graphBuilder;
        Matrix<DataType> edgeHeatMap{};
    };
}

#endif //TEMPO_GNNEDGEPOLARITYPREDICTOR_HPP
