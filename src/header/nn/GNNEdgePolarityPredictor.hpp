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
        template<concepts::scalar T, SchedulingResource R>
        GNNEdgePolarityPredictor(const std::filesystem::path &modelLocation,
                                 const std::filesystem::path &featureExtractorConfigLocation,
                                 const SchedulingProblemHelper<T, R> &problemInstance) :
                gnn(modelLocation)/*,
                graphBuilder(featureExtractorConfigLocation, problemInstance)*/ {}

        /**
         * Returns the choice point in the appropriate polarity according to a GNN edge regressor
         * @param cp choicepoint to evaluate
         * @return cp or ~cp according to the GNN heat map
         */
        /*template<typename Sched>
        requires(traits::is_same_template_v<DistanceConstraint, decltype(std::declval<Sched>().getEdge(
                std::declval<lit>()))>)
        lit choose(tempo::var cp, const Sched &scheduler) const{
            // @TODO does not work anymore!!!
            auto [from, to] = scheduler.getEdge(POS(cp));
            return choosePolarityFromHeatMap(from, to, edgeHeatMap) ? POS(cp) : NEG(cp);
        }*/

        /**
         * Runs the GNN to produce an edge heat map
         * @param scheduler scheduler instance for which to produce the heat map
         */
        /*template<concepts::scalar T>
        void preEvaluation(const Scheduler<T> &scheduler) {
            auto graph = graphBuilder.getGraph(
                    makeSolverState([&](event from, event to) { return scheduler.distance(from, to); }));
            edgeHeatMap = gnn.getHeatMap(graph);
        }*/

    protected:
        static bool choosePolarityFromHeatMap(unsigned taskFrom, unsigned taskTo, const Matrix<DataType> &heatMap);

    private:
        EdgeRegressor gnn;
        //GraphBuilder graphBuilder;
        Matrix<DataType> edgeHeatMap{};
    };
}

#endif //TEMPO_GNNEDGEPOLARITYPREDICTOR_HPP
