/**
* @author Tim Luchterhand
* @date 18.09.24
* @brief GNN based relaxation policy based on edge predictor
*/

#ifndef TEMPO_GNNRELAXATIONPOLICY_HPP
#define TEMPO_GNNRELAXATIONPOLICY_HPP

#include <filesystem>
#include <vector>
#include <ranges>
#include <iostream>

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

    /**
     * @brief GNN based relaxation policy. Uses precedence predictor to fix most probable precedences
     * @tparam T timing type
     * @tparam R resource type
     */
    template<concepts::scalar T, SchedulingResource R>
    class GNNRelaxationPolicy {
        GNNPrecedencePredictor<T, R> predictor;
        double relaxationRatio;
        double relaxationDecay;
        double minCertainty;
        const Solver<T> &solver;
        tempo::util::Profiler profiler;
        std::vector<Literal<T>> assumptionCache;

    public:
        GNNRelaxationPolicy(const GNNRelaxationPolicy &) = default;
        GNNRelaxationPolicy(GNNRelaxationPolicy &&) = default;

        /**
         * Dtor. Prints timing information if solver verbosity allows it
         */
        ~GNNRelaxationPolicy() {
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
        GNNRelaxationPolicy(const Solver<T> &solver, const fs::path &modelLocation,
                            const fs::path &featureExtractorConfigLocation,
                            const SchedulingProblemHelper<T, R> &problemInstance,
                            double relaxationRatio, double relaxationDecay, double minCertainty) :
                predictor(modelLocation, featureExtractorConfigLocation, problemInstance,
                          problemInstance.getSearchLiterals(solver)),
                relaxationRatio(relaxationRatio), relaxationDecay(relaxationDecay), minCertainty(minCertainty),
                solver(solver) {
            if (relaxationRatio < 0 or relaxationRatio > 1) {
                throw std::runtime_error("invalid relaxation ratio");
            }

            if (relaxationDecay < 0 or relaxationDecay > 1) {
                throw std::runtime_error("invalid relaxation ratio");
            }

            updateCache();
        }

        /**
         * Relaxation policy interface. Selects a set of literals according to the last GNN prediction
         * @param literals out param literals
         */
        void select(std::vector<Literal<T>> &literals) const {
            auto numLiterals = static_cast<std::size_t>(predictor.numLiterals() * relaxationRatio);
            if (numLiterals == 0) {
                literals.clear();
                return;
            }

            literals = assumptionCache;
            if (solver.getOptions().verbosity >= Options::NORMAL) {
                std::cout << "-- trying to fix " << literals.size() << " literals" << std::endl;
            }
        }

        /**
         * Call this method when the last relaxation was a success. Currently does nothing
         */
        void notifySuccess() {

        }

        void updateCache() {
            using namespace std::views;
            auto numLiterals = static_cast<std::size_t>(predictor.numLiterals() * relaxationRatio);
            if (numLiterals == 0) {
                return;
            }

            tempo::util::ScopeWatch sw(profiler, "gnn lns policy update");
            {
                tempo::util::ScopeWatch _(profiler, "gnn inference");
                predictor.updateConfidence(solver);
            }

            auto lits = predictor.getLiterals();
            auto selection = lits | filter([m = minCertainty](auto &tpl) { return std::get<1>(tpl) > m; }) |
                             take(numLiterals) | elements<0> | common;
            assumptionCache = std::vector(selection.begin(), selection.end());
        }

        /**
         * Call this method when the last relaxation lead to an UNSAT problem. Updates the GNN prediction and decreases
         * the relaxation ratio
         */
        void notifyFailure() {
            relaxationRatio *= relaxationDecay;
            auto numLiterals = static_cast<std::size_t>(predictor.numLiterals() * relaxationRatio);
            if (numLiterals > 0 and solver.getOptions().verbosity >= Options::NORMAL) {
                std::cout << std::setprecision(2) << "-- failed to relax, setting relaxation ratio = "
                          << relaxationRatio * 100 << "%" << std::endl;
            }

            updateCache();
        }
    };
}

#endif //TEMPO_GNNRELAXATIONPOLICY_HPP
