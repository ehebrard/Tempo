/**
* @author Tim Luchterhand
* @date 08.07.24
* @brief
*/

#ifndef TEMPO_GNNVALUEHEURISTICS_HPP
#define TEMPO_GNNVALUEHEURISTICS_HPP

#include <filesystem>

#include "util/traits.hpp"
#include "Global.hpp"
#include "nn/GNNEdgePolarityPredictor.hpp"
#include "BaseBooleanHeuristic.hpp"
#include "util/SchedulingProblemHelper.hpp"

namespace fs = std::filesystem;

namespace tempo {
    template<typename T>
    class Solver;
}

namespace tempo::heuristics {

    /**
     * @brief Full guidance GNN based value branching heuristic.
     * @detail @copybrief
     * The GNN is called at every decision.
     * @tparam T timing type
     * @tparam R resource type
     * @note This Value branching heuristic is only ment to be used for scheduling problems
     */
    template<concepts::scalar T, SchedulingResource R>
    class GNNFullGuidance: public BaseBooleanHeuristic<GNNFullGuidance<T, R>> {
        nn::heuristics::GNNEdgePolarityPredictor<T, R> polarityPredictor;
    public:
        /**
         * Ctor
         * @param epsilon epsilon greedy value (see BaseBooleanHeuristic)
         * @param modelLocation location of the trained GNN model
         * @param featureExtractorConfigLocation location of the feature extractor configuration
         * (tempo::nn:GraphBuilderConfig)
         * @param problem initial description of the problem
         */
        GNNFullGuidance(double epsilon, const fs::path &modelLocation, const fs::path &featureExtractorConfigLocation,
                        SchedulingProblemHelper <T, R> problem) : BaseBooleanHeuristic<GNNFullGuidance<T, R>>(epsilon),
                                                                  polarityPredictor(modelLocation,
                                                                                    featureExtractorConfigLocation,
                                                                                    std::move(problem)) {}

        /**
         * BaseBooleanHeuristic interface
         * @param x variable id
         * @param solver solver for which to generate literal
         * @returns branching Literal
         */
        auto choose(var_t x, const Solver<T> &solver) -> Literal<T> {
            polarityPredictor.preEvaluation(solver);
            return polarityPredictor.choose(x, solver);
        }
    };
}

#endif //TEMPO_GNNVALUEHEURISTICS_HPP
