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

namespace tempo::nn::heuristics {

    namespace detail {
        bool choosePolarityFromHeatMap(unsigned taskFrom, unsigned taskTo, const Matrix<DataType> &heatMap);
    }


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
            auto edge = solver.boolean.getEdge(x, true);
            auto edgeRev = solver.boolean.getEdge(x, false);
            unsigned from = graphBuilder.getProblem().getMapping()(edge.from);
            unsigned to = graphBuilder.getProblem().getMapping()(edge.to);
            unsigned rfrom = graphBuilder.getProblem().getMapping()(edgeRev.from);
            unsigned rto = graphBuilder.getProblem().getMapping()(edgeRev.to);
            if (from != rto or to != rfrom) {
                throw std::runtime_error("positive and negative edges have different tasks as endpoints.");
            }

            return solver.boolean.getLiteral(detail::choosePolarityFromHeatMap(from, to, edgeHeatMap), x);
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
