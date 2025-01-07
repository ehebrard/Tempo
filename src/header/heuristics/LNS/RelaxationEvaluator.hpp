/**
* @author Tim Luchterhand
* @date 06.01.25
* @file RelaxationEvaluator.hpp
* @brief Relaxation policy evaluation wrapper
*/

#ifndef RELAXATIONEVALUATOR_HPP
#define RELAXATIONEVALUATOR_HPP


#include <vector>
#include <ranges>

#include "util/traits.hpp"
#include "relaxation_interface.hpp"
#include "Solution.hpp"
#include "Model.hpp"

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
        std::vector<std::pair<std::size_t, std::size_t>> literalStats{};
        std::vector<unsigned> discrepancy{};
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
            unsigned discr = 0;
            for (auto [varId, pol] : iterators::const_enumerate(solution.boolean, var_t(0))) {
                discr += pol != proxy.getSolver().boolean.value(BooleanVar(varId));
            }

            discrepancy.emplace_back(discr);
            AssumptionCollector<T, AI> ac(proxy);
            policy.relax(ac);
            const auto &assumptions = ac.getAssumptions();
            if (assumptions.empty()) {
                literalStats.emplace_back(0, 0);
                return;
            }

            totalSet += assumptions.size();
            unsigned numErrors = 0;
            for (auto lit : assumptions) {
                totalErrors += not solution.boolean.consistent(lit);
                numErrors += not solution.boolean.consistent(lit);
            }

            literalStats.emplace_back(numErrors, assumptions.size());
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
        [[nodiscard]] auto assumptionAccuracyPerRun() const {
            return literalStats | std::views::transform([](const auto &pair) {
                return 1 - static_cast<double>(pair.first) / pair.second;
            });
        }

        [[nodiscard]] auto assumptionsPerRun() const noexcept -> const std::vector<std::pair<std::size_t, std::size_t>>& {
            return literalStats;
        }

        [[nodiscard]] auto solutionDiscrepancy() const noexcept -> const std::vector<unsigned>& {
            return discrepancy;
        }

        /**
         * Success status for each run
         * @return
         */
        [[nodiscard]] const std::vector<bool> & runStatus() const noexcept { return successfulRuns; }

        /**
         * Overall relaxation accuracy, i.e. the ratio between the number of correctly fixed literals and all fixed
         * literals over all runs
         * @return
         */
        [[nodiscard]] double totalAssumptionAccuracy() const noexcept { return 1 - static_cast<double>(totalErrors) / totalSet; }

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
