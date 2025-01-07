/**
* @author Tim Luchterhand
* @date 06.01.25
* @file RelaxationEvaluator.hpp
* @brief Relaxation policy evaluation wrapper
*/

#ifndef RELAXATIONEVALUATOR_HPP
#define RELAXATIONEVALUATOR_HPP


#include "util/traits.hpp"
#include "relaxation_interface.hpp"
#include "Solution.hpp"

namespace tempo::lns {
    /**
     * @brief Relaxation policy wrapper that tracks heuristic performance relative to a known solution
     * @tparam Policy Base policy type
     * @tparam T timing type
     */
    template<relaxation_policy Policy, concepts::scalar T = int>
    class RelaxationEvaluator {
        Policy policy;
        Solution<T> solution;
        std::vector<double> assumptionAccuracy{};
        std::vector<bool> successfulRuns{};
        unsigned long totalSet = 0;
        unsigned long totalErrors = 0;
        unsigned numFailedRuns = 0;

    public:
        /**
         * Ctor
         * @tparam Args arguments types for base heuristic ctor
         * @param solution reference solution
         * @param args arguments for base heuristic ctor
         */
        template<typename... Args>
        explicit RelaxationEvaluator(Solution<T> solution, Args &&... args) : policy(std::forward<Args>(args)...),
                                                                              solution(std::move(solution)) {}

        /**
         * Relaxation interface
         * @tparam AI
         * @param proxy
         */
        template<assumption_interface AI>
        void relax(AI &proxy) {
            AssumptionCollector<T, AI> ac(proxy);
            policy.relax(ac);
            const auto &assumptions = ac.getAssumptions();
            if (assumptions.empty()) {
                assumptionAccuracy.emplace_back(std::numeric_limits<double>::quiet_NaN());
                return;
            }

            totalSet += assumptions.size();
            unsigned numErrors = 0;
            for (auto lit : assumptions) {
                totalErrors += not solution.boolean.consistent(lit);
                numErrors += not solution.boolean.consistent(lit);
            }

            assumptionAccuracy.emplace_back(1.0 - static_cast<double>(numErrors) / assumptions.size());
        }

        void notifySuccess(unsigned numFails) {
            policy.notifySuccess(numFails);
            successfulRuns.emplace_back(true);
        }

        void notifyFailure(unsigned numFails) {
            policy.notifyFailure(numFails);
            successfulRuns.emplace_back(false);
            ++numFailedRuns;
        }

        /**
         * Relaxation accuracy (ratio between correctly fixed literals and all fixed literals) for each run
         * @return
         */
        [[nodiscard]] const std::vector<double> & assumptionAccuracyPerRun() const { return assumptionAccuracy; }

        /**
         * Success status for each run
         * @return
         */
        [[nodiscard]] const std::vector<bool> & runStatus() const { return successfulRuns; }

        /**
         * Overall relaxation accuracy, i.e. the ratio between the number of correctly fixed literals and all fixed
         * literals over all runs
         * @return
         */
        [[nodiscard]] double totalAssumptionAccuracy() const { return 1 - static_cast<double>(totalErrors) / totalSet; }

        /**
         * Number of successful runs over all runs
         * @return
         */
        [[nodiscard]] double runAccuracy() const noexcept {
            return 1 - static_cast<double>(numFailedRuns) / successfulRuns.size();
        }
    };

    /**
     * Helper function for creating an evaluating relaxation heuristic from an arbitrary relaxation policy
     * @tparam Policy policy type
     * @tparam T timing type
     * @param policy base relaxation policy
     * @param solution reference solution
     * @return
     */
    template<relaxation_policy Policy, concepts::scalar T = int>
    auto make_evaluator(Policy &&policy, Solution<T> solution) {
        return RelaxationEvaluator<Policy, T>(std::move(solution), std::forward<Policy>(policy));
    }
}

#endif //RELAXATIONEVALUATOR_HPP
