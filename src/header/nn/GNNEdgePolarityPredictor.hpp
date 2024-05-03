/**
 * @author Tim Luchterhand
 * @date 26.07.23.
 */

#ifndef TEMPO_GNNEDGEPOLARITYPREDICTOR_HPP
#define TEMPO_GNNEDGEPOLARITYPREDICTOR_HPP

#include "GNN.hpp"
#include "GraphBuilder.hpp"
#include "util/Matrix.hpp"
#include "util/parsing/format.hpp"

namespace tempo::nn::heuristics {
    /**
     * GNN based value selection heuristic that predicts a probability for each edge and returns the polarity of
     * choicepoint accordingly
     */
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
                                 const ProblemInstance &problemInstance);

        /**
         * Returns the choice point in the appropriate polarity according to a GNN edge regressor
         * @param cp choicepoint to evaluate
         * @return cp or ~cp according to the GNN heat map
         */
        template<concepts::scalar T>
        auto choosePolarity(const DistanceConstraint<T> &cp) const -> DistanceConstraint<T> {
            return choosePolarityFromHeatMap(cp.from, cp.to, edgeHeatMap) ? cp : ~cp;
        }

        /**
         * Runs the GNN to produce an edge heat map
         * @param scheduler scheduler instance for which to produce the heat map
         */
        template<concepts::scalar T>
        void preEvaluation(const Scheduler<T> &scheduler) {
            auto [choices, implied] = scheduler.getChoices();
            std::vector<DistanceConstraint<T>> allChoices = std::move(choices);
            std::copy(implied.begin(), implied.end(), std::back_inserter(allChoices));
            auto graph = graphBuilder.getGraph(scheduler.distance, allChoices);
            edgeHeatMap = gnn.getHeatMap(graph);
        }

    protected:
        static bool choosePolarityFromHeatMap(event from, event to, const Matrix<DataType> &heatMap);

    private:
        EdgeRegressor gnn;
        GraphBuilder graphBuilder;
        Matrix<DataType> edgeHeatMap{};
    };
}

#endif //TEMPO_GNNEDGEPOLARITYPREDICTOR_HPP
