/**
* @author Tim Luchterhand
* @date 06.01.25
* @file RelaxationEvaluator.hpp
* @brief Relaxation policy evaluation wrapper
*/

#ifndef RELAXATION_EVALUATORS_HPP
#define RELAXATION_EVALUATORS_HPP


#include <vector>
#include <ranges>
#include <Iterators.hpp>

#include "util/traits.hpp"
#include "util/Profiler.hpp"
#include "util/parsing/scheduling_collection.hpp"
#include "util/Options.hpp"
#include "util/ThreadPool.hpp"
#include "util/enum.hpp"
#include "util/KillHandler.hpp"
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
        std::vector<unsigned> assumptionTime{};
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

            AssumptionCollector<T, AI> ac(proxy);
            stopWatch.start();
            policy.relax(ac);
            assumptionTime.emplace_back(stopWatch.elapsed<std::chrono::milliseconds>());
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

            stopWatch.start();
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
        const auto &assumptionTimePerRun() const noexcept { return assumptionTime; }
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

    /**
     * @brief Local search result
     */
    PENUM(RegionStatus, Optimal, Unsat, Cancelled, Empty)

    template<concepts::scalar T>
    using RegionResult = std::pair<T, RegionStatus>;

    /**
     * @brief When to do local search
     */
    enum class ExecutionPolicy {
        Blocking, ///< local search pauses the actual search for its duration
        Lazy, ///< local search is deferred until the result is needed
        MultiThreaded ///< local search is performed during search in a separate thread pool
    };

    /**
     * @brief Relaxation policy wrapper that performs a deep search for each sub problem obtained by relaxation
     * @tparam T timing type
     * @tparam Policy base policy type
     */
    template<concepts::scalar T, relaxation_policy Policy>
    class RelaxationRegionEvaluator {
        Policy basePolicy;
        Options options;
        ThreadPool tp;
        std::vector<std::future<RegionResult<T>>> localOptima{};
        unsigned timeout;
        ExecutionPolicy policy;
    public:
        /**
         * Ctor
         * @tparam Args policy ctor argument types
         * @param options solver options
         * @param timeout timout for local search in ms (if 0, no local search is performed)
         * @param executionPolicy execution policy for running local search
         * @param numThreads number of threads to use for local searches (if 0, local search is performed at end of search)
         * @param args arguments for base policy ctor
         */
        template<typename... Args>
        RelaxationRegionEvaluator(Options options, unsigned timeout, ExecutionPolicy executionPolicy,
                                  unsigned numThreads, Args &&... args)
            : basePolicy(std::forward<Args>(args)...), options(std::move(options)),
              tp(executionPolicy != ExecutionPolicy::MultiThreaded ? 0 : numThreads),
              timeout(timeout), policy(executionPolicy) {}

        /**
         * Gets the underlying relaxation policy
         * @return reference to underlying policy
         */
        auto getBasePolicy() -> auto & {
            return basePolicy;
        }

        /**
         * @copydoc getBasePolicy
         */
        auto getBasePolicy() const -> const auto & {
            return basePolicy;
        }

        template<assumption_interface AI>
        void relax(AI &proxy) {
            using enum RegionStatus;
            if (timeout == 0) {
                basePolicy.relax(proxy);
                return;
            }

            AssumptionCollector<T, AI> ac(proxy);
            basePolicy.relax(ac);
            auto job = [this, state = ac.getState(),
                        assumptions = std::move(ac.getAssumptions())]() -> RegionResult<T> {
                if (state == AssumptionState::Fail) {
                    return {0, Unsat};
                } else if (state == AssumptionState::Empty) {
                    return {0, Empty};
                }

                auto opt = options;
                opt.verbosity = Options::SILENT;
                auto p = loadSchedulingProblem(opt);
                AssumptionProxy s(*p.solver);
                s.makeAssumptions(assumptions);
                assert(s.getState() == AssumptionState::Success);
                bool cancelled = false;
                p.solver->PropagationCompleted.subscribe_unhandled(
                    [this, sw = util::StopWatch{}, &cancelled, &solver = *p.solver](auto &) {
                        if (sw.elapsed<std::chrono::milliseconds>() > this->timeout) {
                            solver.cancelSearch();
                            cancelled = true;
                        }
                    });
                p.solver->minimize(p.instance.schedule().duration);
                cancelled |= KillHandler::instance().signalReceived();
                if (p.solver->numeric.hasSolution()) {
                    auto makespan = p.solver->numeric.lower(p.instance.schedule().duration);
                    return {makespan, cancelled ? Cancelled : Optimal};
                }

                return {0, cancelled ? Cancelled : Unsat};
            };

            switch (policy) {
                case ExecutionPolicy::Blocking:
                    localOptima.emplace_back(std::async(std::launch::deferred, std::move(job)));
                    localOptima.back().wait();
                    break;
                case ExecutionPolicy::Lazy:
                    localOptima.emplace_back(std::async(std::launch::deferred, std::move(job)));
                    break;
                case ExecutionPolicy::MultiThreaded:
                    localOptima.emplace_back(tp.submit(std::move(job)));
                    break;
            }
        }

        void notifySuccess(unsigned numFails) {
            basePolicy.notifySuccess(numFails);
        }

        void notifyFailure(unsigned numFails) {
            basePolicy.notifyFailure(numFails);
        }

        /**
         * Gets the results if the local search
         * @return reference to vector with local search results
         * @note Results are run asynchronously and need to be awaited
         */
        auto getResults() const -> const auto & {
            return localOptima;
        }

        /**
         * @copydoc getResults
         */
        auto getResults() -> auto & {
            return localOptima;
        }

    };

    /**
     * Helper function for creating a region evaluating relaxation heuristic from an arbitrary relaxation policy
     * @tparam T timing type
     * @tparam Policy base policy type
     * @param policy base relaxation policy
     * @param options solver options
     * @param timeout timout for local search
     * @param execPolicy local search execution policy
     * @param numThreads number of threads to use for local search
     * @return
     */
    template<concepts::scalar T, relaxation_policy Policy>
    auto make_region_evaluator(Policy &&policy, Options options, unsigned timeout, ExecutionPolicy execPolicy,
                               unsigned numThreads) {
        return RelaxationRegionEvaluator<T, Policy>(std::move(options), timeout, execPolicy, numThreads,
                                                    std::forward<Policy>(policy));
    }
}

#endif //RELAXATION_EVALUATORS_HPP
