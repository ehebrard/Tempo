/**
* @author Tim Luchterhand
* @date 18.09.24
* @brief GNN based relaxation policy based on edge predictor
*/

#ifndef TEMPO_GNNREPAIR_HPP
#define TEMPO_GNNREPAIR_HPP

#include <filesystem>
#include <vector>
#include <ranges>
#include <iostream>
#include <limits>
#include <variant>
#include <Iterators.hpp>

#include "heuristics/LNS/relaxation_interface.hpp"
#include "util/traits.hpp"
#include "util/SchedulingProblemHelper.hpp"
#include "Literal.hpp"
#include "GNNPrecedencePredictor.hpp"
#include "util/Profiler.hpp"
#include "util/Options.hpp"
#include "util/factory_pattern.hpp"
#include "Solver.hpp"
#include "Constant.hpp"
#include "heuristics/LNS/fix_policies.hpp"
#include "heuristics/LNS/PolicyDecay.hpp"

namespace tempo::nn {
    namespace fs = std::filesystem;


    /**
     * @brief GNN based relaxation policy. Uses precedence predictor to fix most probable precedences
     * @tparam T timing type
     * @tparam R resource type
     */
    template<concepts::scalar T, SchedulingResource R>
    class GNNRepair {
        using FixPolicy = lns::VariantFix<lns::BestN, lns::GreedyFix<T>, lns::OptimalFix<T>>;
        GNNPrecedencePredictor<T, R> predictor;
        const Solver<T> &solver;
        mutable tempo::util::Profiler profiler{};
        const lns::PolicyDecayConfig decayConfig;
        std::vector<std::pair<Literal<T>, double>> gnnCache;
        FixPolicy fixPolicy{};
        unsigned failCount = 0;
        unsigned solverFailCount = 0;
        double fixRatio;
        std::size_t numFixed = 0;

        auto calcFailRatio(auto fails) const noexcept {
            return static_cast<double>(fails - solverFailCount) / predictor.numLiterals();
        }

        auto decayFactor(auto failRatio) const {
            using enum lns::DecayMode;
            switch(decayConfig.decayMode) {
                case Constant:
                    return decayConfig.decay;
                case Reciprog:
                    return std::min(decayConfig.decay / failRatio * decayConfig.maxFailRatio, decayConfig.decay);
                case Exponential:
                    return std::min(std::pow(decayConfig.decay, failRatio / decayConfig.maxFailRatio), decayConfig.decay);
                default:
                    throw std::runtime_error("enum out of bounds");
            }
        }

    public:
        GNNRepair(const GNNRepair &) = default;
        GNNRepair(GNNRepair &&) = default;

        /**
         * Dtor. Prints timing information if solver verbosity allows it
         */
        ~GNNRepair() {
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
         * @param decayConfig config struct for policy decay
         * @param assumptionMode how to make assumptions
         */
        GNNRepair(const Solver<T> &solver, const fs::path &modelLocation,
                  const fs::path &featureExtractorConfigLocation,
                  const SchedulingProblemHelper<T, R> &problemInstance,
                  const lns::PolicyDecayConfig &decayConfig, lns::AssumptionMode assumptionMode) :
                predictor(modelLocation, featureExtractorConfigLocation, problemInstance,
                          problemInstance.getSearchLiterals(solver)),
                solver(solver), decayConfig(decayConfig), fixRatio(decayConfig.fixRatio) {
            if (decayConfig.fixRatio < 0 or decayConfig.fixRatio > 1) {
                throw std::runtime_error("invalid fix ratio");
            }

            if (decayConfig.decay < 0 or decayConfig.decay >= 1) {
                throw std::runtime_error("invalid decay");
            }

            if (solver.getOptions().verbosity >= Options::YACKING) {
                std::cout << decayConfig << std::endl;
            }

            using enum lns::AssumptionMode;
            switch(assumptionMode) {
                case GreedySkip:
                    fixPolicy.template emplace<lns::GreedyFix<T>>(false);
                    break;
                case GreedyInverse:
                    fixPolicy.template emplace<lns::GreedyFix<T>>(true);
                    break;
                case Optimal:
                    fixPolicy.template emplace<lns::OptimalFix<T>>(100, 50, false);
                    break;
                case BestN:
                    break;
            }
        }

        /**
         * Repair interface. Fixes a number of literals based on the current fix ratio
         * @tparam AssumptionInterface
         * @param s
         */
        template<lns::assumption_interface AssumptionInterface>
        void fix(AssumptionInterface &s) {
            const auto numLits = maxNumLiterals();
            if (numLits == 0) {
                return;
            }

            tempo::util::ScopeWatch sw(profiler, "repair");
            numFixed = fixPolicy.select(s, numLits, failCount, gnnCache);
            if (solver.getOptions().verbosity >= Options::YACKING) {
                if (s.getState() == lns::AssumptionState::Fail) {
                    std::cout << "-- failed to fix literals\n";
                } else {
                    std::cout << "-- fixing " << numFixed << " / " << predictor.numLiterals()
                              << " literals" << (failCount > 0 ? " after fail\n" : "\n");
                }
            }
        }

        /**
         * Call this method when the last relaxation was a success. Currently does nothing
         */
        void notifySuccess(unsigned fails) {
            const auto failRatio = calcFailRatio(fails);
            if (failRatio > decayConfig.maxFailRatio and decayConfig.decreaseOnSuccess) {
                fixRatio *= decayFactor(failRatio);
                if (solver.getOptions().verbosity >= Options::YACKING) {
                    std::cout << "-- decreasing fix ratio to " << fixRatio * 100
                              << "% after too many solver fails" << std::endl;
                }
            } else if (failRatio < decayConfig.minFailRatio) {
                fixRatio = std::min(1.0, fixRatio / decayConfig.decay);
                if (solver.getOptions().verbosity >= Options::YACKING) {
                    std::cout << "-- increasing fix ratio to " << fixRatio * 100
                              << "%" << std::endl;
                }
            }

            solverFailCount = fails;
            failCount = 0;
        }

        /**
         * Runs the GNN inference and fills the assumption cache
         */
        void runInference() {
            using namespace std::views;
            if (maxNumLiterals() == 0) {
                gnnCache.clear();
                return;
            }

            tempo::util::ScopeWatch sw(profiler, "gnn lns policy update");
            predictor.updateConfidence(solver);
            auto lits = predictor.getLiterals();
            auto selection = lits | take_while([m = decayConfig.minCertainty](auto &tpl) { return std::get<1>(tpl) > m; })
                             | common;
            gnnCache = std::vector(selection.begin(), selection.end());
            if (solver.getOptions().verbosity >= Options::YACKING) {
                std::cout << "-- Updating GNN confidence values" << std::endl;
            }
        }

        /**
         * Maximum number of literals the predictor will try to fix
         * @return upper bound on number of literals to fix
         */
        auto maxNumLiterals() const noexcept {
           return static_cast<std::size_t>(predictor.numLiterals() * fixRatio);
        }

        /**
         * Call this function to reset confidences and caches and rerun inference. This is useful when the gnn will be
         * applied to a new region in the search space
         * @param region The literals that were already fixed
         */
        template<concepts::typed_range<Literal<T>> Region>
        void notifyNewRegion(const Region &region) {
            failCount = 0;
            if (std::ranges::empty(region)) {
                fixRatio = 0;
                return;
            }

            predictor.reinitialize(solver);
            fixRatio = decayConfig.fixRatio;
            runInference();
            fixPolicy.reset();
        }

        /**
         * Indicates whether a new region should be explored
         * @return true if ratio of literals to fix is lower than exhaustion threshold, false otherwise
         */
        bool exhausted() const noexcept {
            auto ratio = std::min(numFixed, maxNumLiterals()) / static_cast<double>(predictor.numLiterals());
            bool res = ratio < decayConfig.exhaustionThreshold;
            if (res and solver.getOptions().verbosity >= Options::SOLVERINFO) {
                std::cout << "-- repair exhausted: num fixed = " << numFixed << ", max num literals = "
                          << maxNumLiterals() << ", ratio = " << ratio << "\n";
            }
            return res;
        }

        /**
         * Call this method when the last relaxation lead to an UNSAT problem. Updates the GNN prediction and decreases
         * the relaxation ratio
         */
        void notifyFailure(unsigned fails) {
            if (++failCount > decayConfig.retryLimit) {
                fixRatio *= decayFactor(calcFailRatio(fails));
                if (maxNumLiterals() > 0 and solver.getOptions().verbosity >= Options::YACKING) {
                    std::cout << std::setprecision(2) << "-- setting fix ratio = "
                              << fixRatio * 100 << "%" << std::endl;
                }
            } else {
                if (solver.getOptions().verbosity >= Options::YACKING) {
                    std::cout << "-- backbone prediction failed, retrying more carefully\n";
                }
            }

            solverFailCount = fails;
        }
    };
}

#endif //TEMPO_GNNREPAIR_HPP
