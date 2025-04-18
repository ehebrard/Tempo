/**
* @author Tim Luchterhand
* @date 25.10.24
* @brief
*/

#ifndef TEMPO_GNN_RELAXATION_HPP
#define TEMPO_GNN_RELAXATION_HPP

#include <filesystem>
#include <stdexcept>
#include <vector>


#include "Literal.hpp"
#include "Solver.hpp"
#include "util/traits.hpp"
#include "util/SchedulingProblemHelper.hpp"
#include "GNNRelaxationPredictor.hpp"
#include "heuristics/LNS/relaxation_interface.hpp"
#include "heuristics/LNS/fix_policies.hpp"
#include "heuristics/LNS/PolicyDecay.hpp"
#include "util/Profiler.hpp"
#include "heuristics/LNS/SolutionPool.hpp"
#include "util/SubscribableEvent.hpp"

namespace tempo::nn {

    /**
     * @brief GNN relaxation policy
     * @tparam T timing type
     * @tparam R resource type
     */
    template<concepts::scalar T, SchedulingResource R>
    class GNNRelax {
        // --- helpers
        using FixPolicy = lns::VariantFix<lns::BestN<lns::OrderType::Ascending>,
                lns::GreedyFix<T, lns::OrderType::Ascending>, lns::SampleFix<true>, lns::TaskFix<T, true>>;
        GNNRelaxationPredictor<T, R> predictor;
        mutable tempo::util::Profiler profiler;
        FixPolicy fixPolicy;
        lns::PolicyDecay policyDecay;
        lns::SolutionPool<T> solutions;
        SubscriberHandle handle;

        // --- state
        std::vector<std::pair<Literal<T>, DataType>> gnnCache;
        std::size_t numFixed = 0;
        double qualityFactor = 1;
        bool newInference = false;

        // --- config
        double exhaustionThreshold;
        double exhaustionProbability;
        int verbosity;
        bool reuseSolutions;

    public:
        GNNRelax(const GNNRelax &) = delete;
        GNNRelax(GNNRelax &&) = delete;
        GNNRelax &operator=(const GNNRelax &) = delete;
        GNNRelax &operator=(GNNRelax &&) = delete;

        /**
         * Ctor
         * @param solver instance of the solver
         * @param modelLocation path to GNN location
         * @param featureExtractorConfigLocation path to feature extractor config
         * @param problemInstance problem instance helper
         * @param decayConfig fix ratio decay policy config
         * @param assumptionMode how to make assumptions
         * @param exhaustionThreshold fix ratio threshold at which a new solution is explored
         * @param exhaustionProbability probability with which a new solution is explored even if not exhausted
         * @param reuseSolutions whether the same solution may be used twice
         * @param sampleSmoothingFactor smoothing factor for sample fix policy
         * @param resourceConstraints Resource expression needed for task fix policy, default empty
         * @param poolSize optional maximum pool size (default unlimited)
         */
        template<resource_range RR = std::vector<NoOverlapExpression<>>>
        GNNRelax(const Solver<T> &solver, const fs::path &modelLocation,
                 const fs::path &featureExtractorConfigLocation, const SchedulingProblemHelper<T, R> &problemInstance,
                 const lns::PolicyDecayConfig &decayConfig, lns::AssumptionMode assumptionMode,
                 double exhaustionThreshold, double exhaustionProbability, bool reuseSolutions = false,
                 double sampleSmoothingFactor = 0, const RR &resourceConstraints = {},
                 std::optional<std::size_t> poolSize = {}) :
                predictor(modelLocation, featureExtractorConfigLocation, problemInstance,
                          problemInstance.getSearchLiterals(solver)),
                policyDecay(decayConfig, predictor.numLiterals(), solver.getOptions().verbosity),
                solutions(problemInstance.schedule().duration, poolSize), handle(solver.SolutionFound.subscribe_handled(
                    [this](const auto &s) { solutions.addSolution(s); })),
                exhaustionThreshold(exhaustionThreshold), exhaustionProbability(exhaustionProbability),
                verbosity(solver.getOptions().verbosity), reuseSolutions(reuseSolutions) {
            using enum lns::AssumptionMode;
            using GF = lns::GreedyFix<T, lns::OrderType::Ascending>;
            switch (assumptionMode) {
                case BestN:
                    break;
                case GreedySkip:
                    fixPolicy.template emplace<GF>(false);
                    break;
                case GreedyInverse:
                    fixPolicy.template emplace<GF>(true);
                    break;
                case Sample:
                    fixPolicy.template emplace<lns::SampleFix<true>>(sampleSmoothingFactor);
                    break;
                case TaskFull:
                    fixPolicy.template emplace<lns::TaskFix<T, true>>(problemInstance.tasks(),
                                                                      resourceConstraints, true);
                    break;
                case TaskReduced:
                    fixPolicy.template emplace<lns::TaskFix<T, true>>(problemInstance.tasks(),
                                                                      resourceConstraints, false);
                    break;
                default:
                    throw std::runtime_error("unsupported assumption mode " + penum_to_string(assumptionMode));
            }

            if (verbosity >= Options::YACKING) {
                std::cout << "-- GNN relaxation policy config\n" << std::boolalpha
                          << "\t-- fix policy " << assumptionMode << "\n"
                          << "\t-- exhaustion threshold: " << exhaustionThreshold << "\n"
                          << "\t-- exhaust probability: " << exhaustionProbability << "\n"
                          << "\t-- reuse solutions: " << reuseSolutions << "\n";
                if (assumptionMode == Sample) {
                    std::cout << "\t-- sample smoothing factor: " << sampleSmoothingFactor << "\n";
                }
                std::cout << decayConfig << std::noboolalpha << std::endl;
            }
        }

        /**
         * Maximum number of literals the predictor will try to fix
         * @return upper bound on number of literals to fix
         */
        auto maxNumLiterals() const noexcept {
            return static_cast<std::size_t>(predictor.numLiterals() * policyDecay.getFixRatio() * qualityFactor);
        }

        /**
         * Relaxation policy interface
         * @tparam AI assumption interface type
         * @param proxy assumption interface
         */
        template<lns::assumption_interface AI>
        void relax(AI &proxy) {
            const auto &solver = proxy.getSolver();
            if (gnnCache.empty() and not runInference(solver)) {
                if (verbosity >= Options::YACKING) {
                    std::cout << "-- solution pool exhausted, doing root search\n";
                }

                return;
            }

            tempo::util::ScopeWatch sw(profiler, "relax");
            numFixed = fixPolicy.select(proxy, maxNumLiterals(), newInference, gnnCache);
            newInference = false;
            if (verbosity >= Options::YACKING) {
                if (proxy.getState() == lns::AssumptionState::Fail) {
                    std::cout << "-- failed to fix literals\n";
                } else {
                    std::cout << "-- fixing " << numFixed << " / " << predictor.numLiterals()
                              << " literals" << (policyDecay.getFailCount() > 0 ? " after fail\n" : "\n");
                }
            }
        }

        /**
         * Relaxation policy interface
         * @param numFailures current number of solver fails
         */
        void notifySuccess(unsigned numFailures) {
            policyDecay.notifySuccess(numFailures);
            exhaustIfNecessary();
        }

        /**
         * Relaxation policy interface
         * @param numFailures current number of solver fails
         */
        void notifyFailure(unsigned numFailures) {
            policyDecay.notifyFailure(numFailures);
            exhaustIfNecessary();
        }

        /**
         * Dtor. Prints profiling information if verbosity allows it
         */
        ~GNNRelax() {
            if (verbosity >= Options::YACKING) {
                profiler.printAll<std::chrono::milliseconds>(std::cout);
            }
        }

    private:
        void exhaustIfNecessary() {
            const auto fixRatio = std::min(numFixed, maxNumLiterals()) / static_cast<double>(predictor.numLiterals());
            bool randomExhaust = random_event_occurred(exhaustionProbability);
            if (fixRatio >= exhaustionThreshold and not randomExhaust) {
                return;
            }

            if (randomExhaust and verbosity >= Options::YACKING) {
                std::cout << "-- randomly trying out new solution\n";
            }

            gnnCache.clear();
            policyDecay.resetFixRatio();
            policyDecay.resetFailCount();
            fixPolicy.reset();
        }

        bool runInference(const Solver<T> &solver) {
            using namespace std::views;
            newInference = true;
            if (maxNumLiterals() == 0 or solutions.empty()) {
                gnnCache.clear();
                return false;
            }


            tempo::util::ScopeWatch sw(profiler, "gnn lns policy update");
            auto sol = reuseSolutions ? solutions.peekLast() : solutions.popLast();
            qualityFactor = solutions.getMakespan(sol) / static_cast<double>(solutions.bestMakespan());
            predictor.updateConfidence(solver, sol);
            gnnCache = predictor.getLiterals();
            if (verbosity >= Options::YACKING) {
                std::cout << "-- Updating GNN confidence values" << std::endl;
            }

            return true;
        }
    };
}

#endif //TEMPO_GNN_RELAXATION_HPP
