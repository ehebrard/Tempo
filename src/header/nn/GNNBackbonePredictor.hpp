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
#include <limits>

#include "heuristics/RelaxationInterface.hpp"
#include "util/traits.hpp"
#include "util/SchedulingProblemHelper.hpp"
#include "Literal.hpp"
#include "GNNPrecedencePredictor.hpp"
#include "util/Profiler.hpp"
#include "util/enum.hpp"
#include "util/Options.hpp"

namespace tempo {
    template<typename T>
    class Solver;
}

namespace tempo::nn {
    namespace fs = std::filesystem;

    PENUM(DecayMode, Constant, Linear, Exponential);

    struct PolicyConfig {
        PolicyConfig() noexcept;

        /**
         * Ctor
         * @param relaxationRatio percentage of literals to fix [0, 1]
         * @param reactivity factor to apply to relaxation ratio on fail
         * @param minCertainty minimum GNN certainty [0, 1]
         * @param minFailRatio lower bound solver failure rate at which to increase relaxation ratio
         * @param maxFailRatio upper bound solver failure rate at which to decrease relaxation ratio
         * @param carefulAssumptions whether to propagate after each literal
         * @param retryLimit number of retries with same relaxation ration before decreasing relaxation ratio
         * @param decayMode type of decay to apply on fail or after too many solver fails
         */
        PolicyConfig(double relaxationRatio, double reactivity, double minCertainty, double minFailRatio,
                     double maxFailRatio, bool carefulAssumptions, unsigned int retryLimit,
                     DecayMode decayMode) noexcept;

        double relaxationRatio, reactivity, minCertainty, minFailRatio, maxFailRatio;
        bool carefulAssumptions;
        unsigned retryLimit;
        DecayMode decayMode;
    };

    std::ostream &operator<<(std::ostream &os, const PolicyConfig &config);

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
        unsigned solverFailCount = 0;

        auto calcFailRatio(auto fails) const noexcept {
            return static_cast<double>(fails - solverFailCount) / predictor.numLiterals();
        }

        auto decayFactor(auto failRatio) const {
            using enum DecayMode;
            switch(config.decayMode) {
                case Constant:
                    return config.reactivity;
                case Linear:
                    return std::min(config.reactivity / failRatio, config.reactivity);
                case Exponential:
                    return std::min(std::pow(config.reactivity, failRatio / config.maxFailRatio), config.reactivity);
                default:
                    throw std::runtime_error("enum out of bounds");
            }
        }

    public:
        GNNBackbonePredictor(const GNNBackbonePredictor &) = default;
        GNNBackbonePredictor(GNNBackbonePredictor &&) = default;

        /**
         * Dtor. Prints timing information if solver verbosity allows it
         */
        ~GNNBackbonePredictor() {
            if (solver.getOptions().verbosity >= Options::YACKING) {
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

            if (config.reactivity < 0 or config.reactivity >= 1) {
                throw std::runtime_error("invalid relaxation ratio");
            }

            if (solver.getOptions().verbosity >= Options::YACKING) {
                std::cout << config << std::endl;
            }

            updateCache();
        }

        template<heuristics::AssumptionInterface<T> AssumptionInterface>
        void relax(AssumptionInterface &s) {
            auto numLiterals = static_cast<std::size_t>(predictor.numLiterals() * config.relaxationRatio);
            if (numLiterals == 0) {
                return;
            }

            tempo::util::ScopeWatch sw(profiler, "make relaxation");
            auto assumptions = assumptionCache | std::views::take(numLiterals);
            if (not config.carefulAssumptions or failCount == 0) {
                s.makeAssumptions(assumptions);
                if (not assumptions.empty() and solver.getOptions().verbosity >= Options::YACKING) {
                    std::cout << "-- fixing " << assumptions.size() << " literals\n";
                }
            } else {
                std::size_t litCount = 0;
                std::vector<Literal<T>> newAssumptions;
                newAssumptions.reserve(numLiterals);
                for (auto lit : assumptions) {
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
                if (litCount > 0 and solver.getOptions().verbosity >= Options::YACKING) {
                    std::cout << "-- fixing " << litCount << " literals after fail\n";
                }
            }
        }

        /**
         * Call this method when the last relaxation was a success. Currently does nothing
         */
        void notifySuccess(unsigned fails) {
            const auto failRatio = calcFailRatio(fails);
            if (failRatio > config.maxFailRatio) {
                config.relaxationRatio *= decayFactor(failRatio);
                updateCache();
                if (solver.getOptions().verbosity >= Options::YACKING) {
                    std::cout << "-- decreasing fix ratio to " << config.relaxationRatio * 100
                              << "% after too many solver fails" << std::endl;
                }
            } else if (failRatio < config.minFailRatio) {
                config.relaxationRatio = std::min(1.0, config.relaxationRatio / config.reactivity);
                if (solver.getOptions().verbosity >= Options::YACKING) {
                    std::cout << "-- increasing fix ratio to " << config.relaxationRatio * 100
                              << "%" << std::endl;
                }
            }

            solverFailCount = fails;
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
                             elements<0> | common;
            assumptionCache = std::vector(selection.begin(), selection.end());
        }

        /**
         * Call this method when the last relaxation lead to an UNSAT problem. Updates the GNN prediction and decreases
         * the relaxation ratio
         */
        void notifyFailure(unsigned fails) {
            if (++failCount > config.retryLimit) {
                config.relaxationRatio *= decayFactor(calcFailRatio(fails));
                auto numLiterals = static_cast<std::size_t>(predictor.numLiterals() * config.relaxationRatio);
                if (numLiterals > 0 and solver.getOptions().verbosity >= Options::YACKING) {
                    std::cout << std::setprecision(2) << "-- setting fix ratio = "
                              << config.relaxationRatio * 100 << "%" << std::endl;
                }
            } else {
                if (solver.getOptions().verbosity >= Options::YACKING) {
                    std::cout << "-- backbone prediction failed, retrying more carefully\n";
                }
            }

            updateCache();
            solverFailCount = fails;
        }
    };
}

#endif //TEMPO_GNNBACKBONEPREDICTOR_HPP
