/**
* @author Tim Luchterhand
* @date 08.07.24
* @brief
*/

#ifndef TEMPO_GNNVALUEHEURISTICS_HPP
#define TEMPO_GNNVALUEHEURISTICS_HPP

#include <filesystem>
#include <iostream>
#include <ostream>

#include "util/traits.hpp"
#include "Global.hpp"
#include "Solver.hpp"
#include "nn/GNNEdgePolarityPredictor.hpp"
#include "heuristics/BaseBooleanHeuristic.hpp"
#include "util/SchedulingProblemHelper.hpp"

namespace fs = std::filesystem;

namespace tempo::nn {

    template<typename D>
    concept dispatcher = requires(D d) {
        { d.runInference() } -> std::convertible_to<bool>;
    };

    /**
     * @brief Full guidance GNN based value branching heuristic.
     * @detail @copybrief
     * The GNN is called at every decision.
     * @tparam T timing type
     * @tparam R resource type
     * @tparam D dispatcher type
     * @note This Value branching heuristic is only ment to be used for scheduling problems
     */
    template<concepts::scalar T, SchedulingResource R, dispatcher D>
    class GNNValueHeuristic: public heuristics::BaseBooleanHeuristic<GNNValueHeuristic<T, R, D>, T> {
        GNNEdgePolarityPredictor<T, R> polarityPredictor;
        D inferenceDispatcher;
    public:

        /**
         * Ctor
         * @tparam DP dispatcher type
         * @param epsilon epsilon greedy value (see BaseBooleanHeuristic)
         * @param modelLocation location of the trained GNN model
         * @param featureExtractorConfigLocation location of the feature extractor configuration
         * (tempo::nn:GraphBuilderConfig)
         * @param inferenceDispatcher dispatcher
         * @param problem initial description of the problem
         */
        template<dispatcher DP>
        GNNValueHeuristic(double epsilon, const fs::path &modelLocation, const fs::path &featureExtractorConfigLocation,
                          DP &&inferenceDispatcher, SchedulingProblemHelper<T, R> problem)
            : heuristics::BaseBooleanHeuristic<GNNValueHeuristic, T>(epsilon),
              polarityPredictor(modelLocation, featureExtractorConfigLocation, std::move(problem)),
              inferenceDispatcher(std::forward<DP>(inferenceDispatcher)) {}

        /**
         * BaseBooleanHeuristic interface
         * @param x variable id
         * @param solver solver for which to generate literal
         * @returns branching Literal
         */
        auto choose(var_t x, const Solver<T> &solver) -> Literal<T> {
            const auto verbosity = solver.getOptions().verbosity;
            if (inferenceDispatcher.runInference()) {
                if (verbosity >= Options::YACKING) {
                    std::cout << "-- rerunning GNN inference" << std::endl;
                }
                polarityPredictor.preEvaluation(solver);
            }

            return polarityPredictor.choose(x, solver);
        }
    };

    /**
     * Helper for template argument deduction
     * @copydoc GNNValueHeuristic::GNNValueHeuristic
     * @return constructed heuristic
     */
    template<concepts::scalar T, SchedulingResource R, dispatcher D>
    auto make_gnn_value_heuristic(double epsilon, const fs::path &modelLocation,
                                  const fs::path &featureExtractorConfigLocation,
                                  D &&inferenceDispatcher, SchedulingProblemHelper<T, R> problem) {
        return GNNValueHeuristic<T, R, D>(epsilon, modelLocation, featureExtractorConfigLocation,
                                          std::forward<D>(inferenceDispatcher), std::move(problem));
    }
}

#endif //TEMPO_GNNVALUEHEURISTICS_HPP
