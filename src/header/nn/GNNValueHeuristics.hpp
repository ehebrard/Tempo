/**
* @author Tim Luchterhand
* @date 08.07.24
* @brief
*/

#ifndef TEMPO_GNNVALUEHEURISTICS_HPP
#define TEMPO_GNNVALUEHEURISTICS_HPP

#include <filesystem>
#include <iostream>

#include "util/traits.hpp"
#include "util/SubscribableEvent.hpp"
#include "Global.hpp"
#include "Solver.hpp"
#include "nn/GNNEdgePolarityPredictor.hpp"
#include "heuristics/BaseBooleanHeuristic.hpp"
#include "util/SchedulingProblemHelper.hpp"

namespace fs = std::filesystem;

namespace tempo::nn {

    /**
     * @brief Full guidance GNN based value branching heuristic.
     * @detail @copybrief
     * The GNN is called at every decision.
     * @tparam T timing type
     * @tparam R resource type
     * @note This Value branching heuristic is only ment to be used for scheduling problems
     */
    template<concepts::scalar T, SchedulingResource R>
    class GNNFullGuidance: public heuristics::BaseBooleanHeuristic<GNNFullGuidance<T, R>, T> {
        GNNEdgePolarityPredictor<T, R> polarityPredictor;
        SubscriberHandle failHandler;
        SubscriberHandle successHandler;
        unsigned numFails = 0;
        unsigned numVars;
        double maxFailRatio;
        bool decisionFlag = false;
    public:

        GNNFullGuidance(const GNNFullGuidance &) = delete;
        GNNFullGuidance(GNNFullGuidance &&) = delete;
        GNNFullGuidance &operator=(const GNNFullGuidance &) = delete;
        GNNFullGuidance &operator=(GNNFullGuidance &&) = delete;
        ~GNNFullGuidance() override = default;

        /**
         * Ctor
         * @param epsilon epsilon greedy value (see BaseBooleanHeuristic)
         * @param modelLocation location of the trained GNN model
         * @param featureExtractorConfigLocation location of the feature extractor configuration
         * (tempo::nn:GraphBuilderConfig)
         * @param solver target solver
         * @param maxFailRatio maximum fail ratio (num decision fails divided by num variables) at which the GNN is
         * reevaluated
         * @param problem initial description of the problem
         */
        GNNFullGuidance(double epsilon, const fs::path &modelLocation, const fs::path &featureExtractorConfigLocation,
                        const Solver<T> &solver, float maxFailRatio, SchedulingProblemHelper<T, R> problem)
            : heuristics::BaseBooleanHeuristic<GNNFullGuidance, T>(epsilon),
              polarityPredictor(modelLocation, featureExtractorConfigLocation, std::move(problem)),
              failHandler(solver.ConflictEncountered.subscribe_handled([this](auto &&) {
                  ++numFails;
                  decisionFlag = false;
              })),
              successHandler(solver.PropagationCompleted.subscribe_handled([this](auto &&) {
                  if (decisionFlag and numFails > 0) { --numFails; }
              })),
              numVars(solver.boolean.size()), maxFailRatio(maxFailRatio) {
            polarityPredictor.preEvaluation(solver);
        }

        double failRatio() const noexcept {
            return numFails / static_cast<double>(numVars);
        }

        /**
         * BaseBooleanHeuristic interface
         * @param x variable id
         * @param solver solver for which to generate literal
         * @returns branching Literal
         */
        auto choose(var_t x, const Solver<T> &solver) -> Literal<T> {
            decisionFlag = true;
            if (failRatio() > maxFailRatio) {
                if (solver.getOptions().verbosity >= Options::YACKING) {
                    std::cout << "-- rerunning GNN inference" << std::endl;
                }
                polarityPredictor.preEvaluation(solver);
                numFails = 0;
            }

            return polarityPredictor.choose(x, solver);
        }
    };
}

#endif //TEMPO_GNNVALUEHEURISTICS_HPP
