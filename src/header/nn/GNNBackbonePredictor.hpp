/**
* @author Tim Luchterhand
* @date 18.09.24
* @brief GNN based relaxation policy based on edge predictor
*/

#ifndef TEMPO_GNNBACKBONEPREDICTOR_HPP
#define TEMPO_GNNBACKBONEPREDICTOR_HPP

#include <filesystem>
#include <vector>
#include <ranges>
#include <iostream>

#include "heuristics/RelaxationInterface.hpp"
#include "util/traits.hpp"
#include "util/SchedulingProblemHelper.hpp"
#include "Literal.hpp"
#include "GNNPrecedencePredictor.hpp"
#include "util/Profiler.hpp"

namespace tempo {
    template<typename T>
    class Solver;
}

namespace tempo::nn {
    namespace fs = std::filesystem;

    struct PolicyConfig {
        double relaxationRatio, relaxationDecay, minCertainty;
        bool carefulAssumptions;
        unsigned retryLimit;
    };

    /**
     * @brief GNN based relaxation policy. Uses precedence predictor to fix most probable precedences
     * @tparam T timing type
     * @tparam R resource type
     */
    template<concepts::scalar T, SchedulingResource R>
    class GNNBackbonePredictor {
        GNNPrecedencePredictor<T, R> predictor;
        const Solver<T> &solver;
        mutable tempo::util::Profiler profiler;
        PolicyConfig config;
        std::vector<Literal<T>> assumptionCache;
        unsigned failCount = 0;

    public:
        GNNBackbonePredictor(const GNNBackbonePredictor &) = default;
        GNNBackbonePredictor(GNNBackbonePredictor &&) = default;

        /**
         * Dtor. Prints timing information if solver verbosity allows it
         */
        ~GNNBackbonePredictor() {
            if (solver.getOptions().verbosity >= Options::NORMAL) {
                profiler.printAll<std::chrono::milliseconds>(std::cout);
            }
        }

        /**
         * Ctor. Performs call to the GNN based on the current solver state
         * @param solver solver in use
         * @param modelLocation path to the model file
         * @param featureExtractorConfigLocation path to the feature extractor config
         * @param problemInstance problem instance
         * @param relaxationRatio percentage of search literals to relay (in [0, 1])
         * @param relaxationDecay decay factor applied to the relaxation ratio on each relaxation fail
         * @param minCertainty minimum GNN prediction certainty
         */
        GNNBackbonePredictor(const Solver<T> &solver, const fs::path &modelLocation,
                             const fs::path &featureExtractorConfigLocation,
                             const SchedulingProblemHelper<T, R> &problemInstance,
                             const PolicyConfig &config) :
                predictor(modelLocation, featureExtractorConfigLocation, problemInstance,
                          problemInstance.getSearchLiterals(solver)), solver(solver), config(config) {
            if (config.relaxationRatio < 0 or config.relaxationRatio > 1) {
                throw std::runtime_error("invalid relaxation ratio");
            }

            if (config.relaxationDecay < 0 or config.relaxationDecay > 1) {
                throw std::runtime_error("invalid relaxation ratio");
            }

            updateCache();
        }

        void relax(heuristics::AssumptionInterface<T> &s) {
            auto numLiterals = static_cast<std::size_t>(predictor.numLiterals() * config.relaxationRatio);
            if (numLiterals == 0) {
                return;
            }

            tempo::util::ScopeWatch sw(profiler, "make relaxation");
            if (not config.carefulAssumptions or failCount == 0) {
                s.makeAssumptions(assumptionCache);
                if (not assumptionCache.empty() and solver.getOptions().verbosity >= Options::NORMAL) {
                    std::cout << "-- fixing " << assumptionCache.size() << " literals\n";
                }
            } else {
                std::size_t litCount = 0;
                std::vector<Literal<T>> newAssumptions;
                newAssumptions.reserve(numLiterals);
                for (auto lit : assumptionCache) {
                    if (litCount == numLiterals) {
                        break;
                    }

                    bool success = s.tryMakeAssumption(lit);
                    litCount += success;
                    if (success) {
                        newAssumptions.emplace_back(lit);
                    }
                }

                std::swap(newAssumptions, assumptionCache);
                if (litCount > 0 and solver.getOptions().verbosity >= Options::NORMAL) {
                    std::cout << "-- fixing " << litCount << " literals after fail\n";
                }
            }
        }

        /**
         * Call this method when the last relaxation was a success. Currently does nothing
         */
        void notifySuccess() {
            failCount = 0;
        }

        void updateCache() {
            using namespace std::views;
            auto numLiterals = static_cast<std::size_t>(predictor.numLiterals() * config.relaxationRatio);
            if (numLiterals == 0) {
                assumptionCache.clear();
                return;
            }

            tempo::util::ScopeWatch sw(profiler, "gnn lns policy update");
            predictor.updateConfidence(solver);
            auto lits = predictor.getLiterals();
            auto selection = lits | filter([m = config.minCertainty](auto &tpl) { return std::get<1>(tpl) > m; }) |
                             take(numLiterals) | elements<0> | common;
            assumptionCache = std::vector(selection.begin(), selection.end());
        }

        /**
         * Call this method when the last relaxation lead to an UNSAT problem. Updates the GNN prediction and decreases
         * the relaxation ratio
         */
        void notifyFailure() {
            if (++failCount > config.retryLimit) {
                config.relaxationRatio *= config.relaxationDecay;
                auto numLiterals = static_cast<std::size_t>(predictor.numLiterals() * config.relaxationRatio);
                if (numLiterals > 0 and solver.getOptions().verbosity >= Options::NORMAL) {
                    std::cout << std::setprecision(2) << "-- setting relaxation ratio = "
                              << config.relaxationRatio * 100 << "%" << std::endl;
                }
            } else {
                if (solver.getOptions().verbosity >= Options::NORMAL) {
                    std::cout << "-- relaxation failed, retrying more carefully\n";
                }
            }

            updateCache();
        }
    };
}

#endif //TEMPO_GNNBACKBONEPREDICTOR_HPP
