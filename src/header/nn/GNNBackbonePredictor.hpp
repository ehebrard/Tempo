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

    PENUM(DecayMode, Constant, Reciprog, Exponential);

    struct PolicyConfig {
        PolicyConfig() noexcept;

        /**
         * Ctor
         * @param fixRatio percentage of literals to fix [0, 1]
         * @param reactivity factor to apply to relaxation ratio on fail
         * @param minCertainty minimum GNN certainty [0, 1]
         * @param minFailRatio lower bound solver failure rate at which to increase relaxation ratio
         * @param maxFailRatio upper bound solver failure rate at which to decrease relaxation ratio
         * @param confidenceUpdateFailRatio fail ratio threshold at which the GNN confidence values are updated
         * @param carefulAssumptions whether to propagate after each literal
         * @param decreaseOnSuccess whether to decrease fix rate even on success
         * @param retryLimit number of retries with same relaxation ration before decreasing relaxation ratio
         * @param decayMode type of decay to apply on fail or after too many solver fails
         */
        PolicyConfig(double fixRatio, double reactivity, double minCertainty, double minFailRatio,
                     double maxFailRatio, double confidenceUpdateFailRatio, bool carefulAssumptions,
                     bool decreaseOnSuccess, DecayMode decayMode, unsigned int retryLimit) noexcept;

        double fixRatio, decay, minCertainty, minFailRatio, maxFailRatio, confidenceUpdateFailRatio;
        bool carefulAssumptions, decreaseOnSuccess;
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
        const PolicyConfig config;
        std::vector<Literal<T>> assumptionCache;
        unsigned failCount = 0;
        unsigned solverFailCount = 0;
        double fixRatio;

        auto calcFailRatio(auto fails) const noexcept {
            return static_cast<double>(fails - solverFailCount) / predictor.numLiterals();
        }

        auto decayFactor(auto failRatio) const {
            using enum DecayMode;
            switch(config.decayMode) {
                case Constant:
                    return config.decay;
                case Reciprog:
                    return std::min(config.decay / failRatio * config.maxFailRatio, config.decay);
                case Exponential:
                    return std::min(std::pow(config.decay, failRatio / config.maxFailRatio), config.decay);
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
                          problemInstance.getSearchLiterals(solver)),
                          solver(solver), config(config), fixRatio(config.fixRatio) {
            if (config.fixRatio < 0 or config.fixRatio > 1) {
                throw std::runtime_error("invalid fix ratio");
            }

            if (config.decay < 0 or config.decay >= 1) {
                throw std::runtime_error("invalid decay");
            }

            if (solver.getOptions().verbosity >= Options::YACKING) {
                std::cout << config << std::endl;
            }

            updateCache();
        }

        template<heuristics::AssumptionInterface AssumptionInterface>
        void relax(AssumptionInterface &s) {
            const auto numLits = numLiterals();
            if (numLits == 0) {
                return;
            }

            tempo::util::ScopeWatch sw(profiler, "make relaxation");
            auto assumptions = assumptionCache | std::views::take(numLits);
            if (not config.carefulAssumptions or failCount == 0) {
                s.makeAssumptions(assumptions);
                if (not assumptions.empty() and solver.getOptions().verbosity >= Options::YACKING) {
                    std::cout << "-- fixing " << assumptions.size() << " literals\n";
                }
            } else {
                std::size_t litCount = 0;
                std::vector<Literal<T>> newAssumptions;
                newAssumptions.reserve(numLits);
                for (auto lit : assumptions) {
                    if (litCount == numLits) {
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
            if (failRatio > config.maxFailRatio and config.decreaseOnSuccess) {
                fixRatio *= decayFactor(failRatio);
                if (solver.getOptions().verbosity >= Options::YACKING) {
                    std::cout << "-- decreasing fix ratio to " << fixRatio * 100
                              << "% after too many solver fails" << std::endl;
                }
            } else if (failRatio < config.minFailRatio) {
                fixRatio = std::min(1.0, fixRatio / config.decay);
                if (solver.getOptions().verbosity >= Options::YACKING) {
                    std::cout << "-- increasing fix ratio to " << fixRatio * 100
                              << "%" << std::endl;
                }
            }

            if (failRatio > config.confidenceUpdateFailRatio) {
                updateCache();
            }

            solverFailCount = fails;
            failCount = 0;
        }

        void updateCache() {
            using namespace std::views;
            if (numLiterals() == 0) {
                assumptionCache.clear();
                return;
            }

            tempo::util::ScopeWatch sw(profiler, "gnn lns policy update");
            predictor.updateConfidence(solver);
            auto lits = predictor.getLiterals();
            auto selection = lits | filter([m = config.minCertainty](auto &tpl) { return std::get<1>(tpl) > m; }) |
                             elements<0> | common;
            assumptionCache = std::vector(selection.begin(), selection.end());
            if (solver.getOptions().verbosity >= Options::YACKING) {
                std::cout << "-- Updating GNN confidence values" << std::endl;
            }
        }

        /**
         * Maximum number of literals the predictor will try to fix
         * @return upper bound on number of literals to fix
         */
        auto numLiterals() const noexcept {
           return static_cast<std::size_t>(predictor.numLiterals() * fixRatio);
        }

        /**
         * Whether the GNN is unable to fix any literals
         * @return true if no literals can be fixed, false else
         */
        bool exhausted() const noexcept {
            return numLiterals() == 0 or assumptionCache.empty();
        }

        /**
         * Call this function to reset confidences and caches. This is useful when the gnn will be applied to a new
         * region in the search space
         */
        void notifyNewRegion() noexcept {
            failCount = 0;
            predictor.resetConfidences();
            fixRatio = config.fixRatio;
        }

        /**
         * Call this method when the last relaxation lead to an UNSAT problem. Updates the GNN prediction and decreases
         * the relaxation ratio
         */
        void notifyFailure(unsigned fails) {
            if (++failCount > config.retryLimit) {
                fixRatio *= decayFactor(calcFailRatio(fails));
                if (numLiterals() > 0 and solver.getOptions().verbosity >= Options::YACKING) {
                    std::cout << std::setprecision(2) << "-- setting fix ratio = "
                              << fixRatio * 100 << "%" << std::endl;
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
