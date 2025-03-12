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
#include "util/SubscribableEvent.hpp"
#include "Global.hpp"
#include "Solver.hpp"
#include "nn/GNNEdgePolarityPredictor.hpp"
#include "heuristics/BaseBooleanHeuristic.hpp"
#include "util/SchedulingProblemHelper.hpp"

namespace fs = std::filesystem;

namespace tempo::nn {
    /**
     * @brief Config struct for dispatcher.
     * @note all fields are in percent relative to the respective threshold except heatDecay.
     */
    struct DispatcherConfig {
        double failIncrement; ///< ratio increment on conflict
        double successDecrement; ///< ratio decrement on successful propagation
        double restartIncrement; ///< ratio increment on search restart
        double solutionDecrement; ///< ratio decrement when solution is found
        double maxFillRate; ///< maximum fill rate
        double heatIncrement; ///< heat increment on inference
        double heatDecay; ///< constant heat decay on step (in %)
        double heatLowerThreshold; ///< threshold under which the heat must fall after overheating
    };

    std::ostream &operator<<(std::ostream &os, const DispatcherConfig &config);

    class Dispatcher {
        static constexpr auto MaxTemperature = 1.0;
        static constexpr auto TriggerThreshold = 10000u;
        // config
        int failIncrement; // ratio increment on conflict
        int successDecrement; // ratio decrement on successful propagation
        int restartIncrement; // ratio increment on search restart
        int solutionDecrement; // ratio decrement when solution is found
        int maxFillRate; // maximum fill rate
        double heatIncrement; // heat increment on inference
        double heatDecay; // constant heat decay on step
        double heatLowerThreshold;
        // state
        double temperature = 0;
        unsigned waterLevel = 0;
        int fillingRate = 0;
        double isOverheated = false;

    public:
        Dispatcher(const DispatcherConfig &config) noexcept;

        void onFail() noexcept;

        void onSuccess() noexcept;

        void onRestart() noexcept;

        void onSolution() noexcept;

        void onInference() noexcept;

        void step() noexcept;

        bool inferenceAllowed() const noexcept;

        bool overheated() const noexcept;

        friend std::ostream &operator<<(std::ostream &os, const Dispatcher &dispatcher);
    };


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
        Dispatcher dispatcher;
        SubscriberHandle failHandler;
        SubscriberHandle successHandler;
        SubscriberHandle restartHandler;
        SubscriberHandle solutionHandler;
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
         * @param dispatcherConfig config for dispatching strategy
         * @param problem initial description of the problem
         */
        GNNFullGuidance(double epsilon, const fs::path &modelLocation, const fs::path &featureExtractorConfigLocation,
                        const Solver<T> &solver, const DispatcherConfig &dispatcherConfig,
                        SchedulingProblemHelper<T, R> problem)
            : heuristics::BaseBooleanHeuristic<GNNFullGuidance, T>(epsilon),
              polarityPredictor(modelLocation, featureExtractorConfigLocation, std::move(problem)),
              dispatcher(dispatcherConfig),
              failHandler(solver.ConflictEncountered.subscribe_handled([this](auto &&) {
                  dispatcher.onFail();
                  decisionFlag = false;
              })),
              successHandler(solver.PropagationCompleted.subscribe_handled([this](auto &&) {
                  if (decisionFlag) { dispatcher.onSuccess(); }
              })),
              restartHandler(solver.SearchRestarted.subscribe_handled([this](auto &&) {
                  dispatcher.onRestart();
              })),
              solutionHandler(solver.SolutionFound.subscribe_handled([this](auto &&) {
                  dispatcher.onSolution();
              })) {
            if (solver.getOptions().verbosity >= Options::NORMAL) {
                std::cout << dispatcherConfig << std::endl;
            }
            polarityPredictor.preEvaluation(solver);
        }

        /**
         * BaseBooleanHeuristic interface
         * @param x variable id
         * @param solver solver for which to generate literal
         * @returns branching Literal
         */
        auto choose(var_t x, const Solver<T> &solver) -> Literal<T> {
            decisionFlag = true;
            dispatcher.step();
            const auto verbosity = solver.getOptions().verbosity;;
            if (verbosity >= Options::SOLVERINFO) {
                std::cout << dispatcher << std::endl;
            }
            if (dispatcher.inferenceAllowed()) {
                dispatcher.onInference();
                if (verbosity >= Options::YACKING) {
                    std::cout << "-- rerunning GNN inference" << std::endl;
                }
                polarityPredictor.preEvaluation(solver);
            } else if (dispatcher.overheated() and verbosity >= Options::YACKING) {
                std::cout << "-- GNN overheated" << std::endl;
            }

            return polarityPredictor.choose(x, solver);
        }
    };
}

#endif //TEMPO_GNNVALUEHEURISTICS_HPP
