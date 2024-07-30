/**
* @author Tim Luchterhand
* @date 24.07.24
* @brief
*/

#ifndef TEMPO_TRACEWATCHER_HPP
#define TEMPO_TRACEWATCHER_HPP

#include <vector>
#include <Iterators.hpp>
#include <cassert>
#include <ranges>
#include <algorithm>

#include "util/traits.hpp"
#include "util/SubscribableEvent.hpp"
#include "Model.hpp"

namespace tempo {

    template<typename T>
    class Solver;

    /**
     * @brief Class that memorizes choices taken on the path to a solution.
     * @details @copybrief Use it to determine if a partial solution is a subset of a previously found solution
     * @note only considers binary literals
     */
    class TraceWatcher {
        std::vector<bool> varPolarity;
        bool onTrack;
        var_t offset;
    public:
        using Conflicts = std::vector<std::pair<var_t, bool>>;

        /**
         * CTor
         * @param vars range containing ids of all search variables
         */
        template<concepts::ctyped_range<var_t> Vars>
        requires(std::ranges::sized_range<Vars>)
        explicit TraceWatcher(const Vars &vars): varPolarity(std::ranges::size(vars), false), onTrack(false),
                                                 offset(std::ranges::min(vars)) {
            if (std::ranges::max(vars) - offset + 1 != varPolarity.size()) {
                throw std::runtime_error("expected continuous range of search variables");
            }
        }

        /**
         * registers a complete solution
         * @tparam TF truth function type
         * @param truthFunction function object that maps variable ids to boolean values
         */
        template<concepts::callable_r<bool, var_t> TF>
        void registerSolution(TF &&truthFunction) {
            for (auto [var, val] : iterators::enumerate(varPolarity, offset)) {
                val = std::forward<TF>(truthFunction)(var);
            }

            onTrack = true;
        }

        /**
         * Checks if still on track to last solution and sets the internal on track flag to false if this is not the
         * case. Also returns a list with variables in conflict with the previous solution along with their
         * corresponding truth values
         * @tparam AF alignment function
         * @param isAligned function object that checks if the truth value of a variable aligns with a given boolean.
         * For example: undefined aligns with true and false but false only aligns with false.
         * @return Vector of pairs of conflicting variables along with their corresponding truth value from the last
         * solution
         */
        template<concepts::callable_r<bool, var_t, bool> AF>
        auto updateOnTrack(AF &&isAligned) -> Conflicts {
            Conflicts ret;
            for (auto [var, val] : iterators::const_enumerate(varPolarity, offset)) {
                if (not std::forward<AF>(isAligned)(var, val)) {
                    ret.emplace_back(var, val);
                }
            }

            onTrack = ret.empty();
            return ret;
        }

        /**
         * Registers a decision
         * @tparam T timing type
         * @param decision decision literal
         */
        template<concepts::scalar T>
        void step(Literal<T> decision) noexcept {
            assert(decision.isBoolean());
            onTrack &= varPolarity[decision.variable() - offset] == decision.sign();
        }

        /**
         * Manually sets the on track flag
         * @param truthVal target value
         */
        void setOnTrack(bool truthVal) noexcept;

        /**
         * Whether still on track to last solution
         * @return true if still on track to last solution, false otherwise
         */
        [[nodiscard]] bool isOnTrack() const noexcept;

        /**
         * Get the last registered solution
         * @return Vector with truth value for each variable (indexed by variable id)
         */
        [[nodiscard]] auto getLastSolution() const noexcept -> const std::vector<bool>&;

        [[nodiscard]] var_t getOffset() const noexcept;

    };

    template<concepts::scalar T>
    class Aligned {
        const Solver<T> &s;
    public:
        explicit constexpr Aligned(const Solver<T> &solver) noexcept: s(solver) {}
        constexpr bool operator()(var_t var, bool truthVal) const noexcept {
            return s.boolean.isUndefined(var) or s.boolean.isTrue(var) == truthVal;
        }
    };


    /**
     * @brief Type of deviation from path to last solution
     */
    enum class DeviationType {
        Propagation, ///< deviation after propagation
        Fail, ///< deviation due to a fail
        Decision ///< deviation occured after decision
    };


    /**
     * @brief Class that follows the search of the solver and triggers an event when the solver deviates from the path
     * to the last solution.
     * @details @copybrief
     */
    class Tracer {
        TraceWatcher watcher;
        SubscriberHandle solutionHandler;
        SubscriberHandle decisionHandler;
        SubscriberHandle conflictHandler;
        SubscriberHandle backtrackHandler;
        SubscriberHandle propagationHandler;
    public:
        SubscribableEvent<DeviationType, TraceWatcher::Conflicts> DeviationOccurred;
        ///< triggered when the solver deviates from path to the last solution. Arguments: deviation type,
        ///< conflicts after propagation

        /**
         * Ctor
         * @tparam T timing type
         * @param solver solver to trace
         */
        template<concepts::scalar T>
        explicit Tracer(const Solver<T> &solver) : watcher(solver.boolean_search_vars),
        solutionHandler(solver.SolutionFound.subscribe_handled([this](const auto &solver) {
            this->handleSolution(solver);
        })),
        decisionHandler(solver.ChoicePoint.subscribe_handled([this](auto lit) {
            this->handleDecision(lit);
        })),
        conflictHandler(solver.ConflictEncountered.subscribe_handled([this](const auto &) {
            this->handleConflict();
        })),
        backtrackHandler(solver.BackTrackCompleted.subscribe_handled([this, &solver]() {
            this->watcher.updateOnTrack(Aligned(solver));
        })),
        propagationHandler(solver.PropagationCompleted.subscribe_handled([this](const auto &solver) {
            this->handlePropagation(solver);
        })) {}

        /**
         * const access to trace watcher
         * @return const ref to internal trace watcher
         */
        [[nodiscard]] auto getWatcher() const noexcept -> const TraceWatcher &;

    private:
        template<concepts::scalar T>
        void handlePropagation(const Solver<T> &solver) {
            if (not watcher.isOnTrack()) {
                return;
            }

            auto conflicts = watcher.updateOnTrack(Aligned(solver));
            if (not watcher.isOnTrack()) {
                DeviationOccurred.trigger(DeviationType::Propagation, std::move(conflicts));
            }
        }

        template<concepts::scalar T>
        void handleSolution(const Solver<T> &solver) {
            watcher.registerSolution([&solver](auto var) {
                assert(not solver.boolean.isUndefined(var));
                return solver.boolean.isTrue(var);
            });
        }

        template<concepts::scalar T>
        void handleDecision(Literal<T> lit) {
            bool ot = watcher.isOnTrack();
            watcher.step(lit);
            if (ot != watcher.isOnTrack()) {
                DeviationOccurred.trigger(DeviationType::Decision, TraceWatcher::Conflicts{});
            }
        }

        void handleConflict();
    };
}

#endif //TEMPO_TRACEWATCHER_HPP
