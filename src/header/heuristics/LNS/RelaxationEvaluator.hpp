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
#include <Iterators.hpp>

#include "util/traits.hpp"
#include "util/Profiler.hpp"
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
        std::optional<Solution<T>> solution;
        std::vector<std::size_t> assumptionSize{};
        std::vector<std::size_t> errors{};
        std::vector<unsigned> numFails{};
        std::vector<unsigned> discrepancy{};
        std::vector<unsigned> solverRunTimes{};
        std::vector<bool> successfulRuns{};
        util::StopWatch stopWatch;
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
        explicit RelaxationEvaluator(std::optional<Solution<T>> solution,
                                     Args &&... args) : policy(std::forward<Args>(args)...),
                                                        solution(std::move(solution)) {}

        /**
         * Relaxation interface
         * @tparam AI
         * @param proxy
         */
        template<assumption_interface AI>
        void relax(AI &proxy) {
            if (solution.has_value()) {
                unsigned discr = 0;
                for (auto [varId, pol] : iterators::const_enumerate(solution->boolean, var_t(0))) {
                    discr += pol != proxy.getSolver().boolean.value(BooleanVar(varId));
                }

                discrepancy.emplace_back(discr);
            }

            stopWatch.start();
            AssumptionCollector<T, AI> ac(proxy);
            policy.relax(ac);
            const auto &assumptions = ac.getAssumptions();
            totalSet += assumptions.size();
            assumptionSize.emplace_back(assumptions.size());
            if (solution.has_value()) {
                unsigned numErrors = 0;
                for (auto lit: assumptions) {
                    totalErrors += not solution->boolean.consistent(lit);
                    numErrors += not solution->boolean.consistent(lit);
                }

                errors.emplace_back(numErrors);
            }
        }

        void notifySuccess(unsigned numFails) {
            policy.notifySuccess(numFails);
            successfulRuns.emplace_back(true);
            auto prev = this->numFails.empty() ? 0u : this->numFails.back();
            this->numFails.emplace_back(numFails - prev);
            solverRunTimes.emplace_back(stopWatch.elapsed<std::chrono::milliseconds>());
        }

        void notifyFailure(unsigned numFails) {
            policy.notifyFailure(numFails);
            successfulRuns.emplace_back(false);
            ++numFailedRuns;
            auto prev = this->numFails.empty() ? 0u : this->numFails.back();
            this->numFails.emplace_back(numFails - prev);
            solverRunTimes.emplace_back(stopWatch.elapsed<std::chrono::milliseconds>());
        }

        /**
         * Relaxation accuracy (ratio between correctly fixed literals and all fixed literals) for each run
         * @return
         */
        [[nodiscard]] auto assumptionAccuracyPerRun() const {
            return iterators::zip(traits::as_mut(errors), traits::as_mut(assumptionSize)) | std::views::transform(
                       [](const auto &tuple) {
                           const auto [e, total] = tuple;
                           return 1 - static_cast<double>(e) / total;
                       });
        }

        [[nodiscard]] auto normalizedAssumptionAccuracyPerRun() const {
            return iterators::zip(traits::as_mut(errors), traits::as_mut(assumptionSize), traits::as_mut(discrepancy)) |
                   std::views::transform(
                       [numVars = solution.has_value() ? solution->boolean.size() : 0ul](const auto &tpl) {
                           const auto [err, total, dsc] = tpl;
                           auto expected = 0.5 - static_cast<double>(dsc) / numVars;
                           return 1 - static_cast<double>(err) / total - expected;
                       });
        }

        [[nodiscard]] auto assumptionsPerRun() const noexcept -> const std::vector<std::size_t>& {
            return assumptionSize;
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

        const auto &failsPerRun() const noexcept { return numFails; }

        const auto &searchTimesPerRun() const noexcept { return solverRunTimes; }
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
    auto make_evaluator(Policy &&policy, std::optional<Solution<T>> solution) {
        return RelaxationEvaluator<Policy, T>(std::move(solution), std::forward<Policy>(policy));
    }
}

#endif //RELAXATIONEVALUATOR_HPP
