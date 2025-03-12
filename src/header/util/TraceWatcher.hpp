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
#include "util/serialization.hpp"
#include "util/enum.hpp"
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
        serialization::Branch varsOnTrack{};
        bool onTrack;
        var_t offset;
    public:
        enum class TruthVal {
            False = 0,
            True = 1,
            Undefined = -1
        };

        using Conflicts = std::vector<std::pair<var_t, bool>>;

        /**
         * CTor
         * @param vars range containing ids of all search variables
         */
        template<concepts::ctyped_range<var_t> Vars>
        requires(std::ranges::sized_range<Vars>)
        explicit TraceWatcher(const Vars &vars): varPolarity(std::ranges::size(vars), false), onTrack(false), offset() {
            if (std::ranges::empty(vars)) {
                throw std::runtime_error("empty variable range");
            }

            auto [min, max] = std::ranges::minmax(vars);
            offset = min;
            if (max - offset + 1 != varPolarity.size()) {
                throw std::runtime_error("expected continuous range of search variables");
            }

            varsOnTrack.reserve(varPolarity.size());
        }

        /**
         * registers a complete solution
         * @tparam TF truth function type
         * @param truthFunction function object that maps variable ids to truth values
         */
        template<concepts::callable_r<TruthVal, var_t> TF>
        void registerSolution(TF &&truthFunction) {
            for (auto [var, val] : iterators::enumerate(varPolarity, offset)) {
                auto tv = std::forward<TF>(truthFunction)(var);
                assert(tv != TruthVal::Undefined);
                val = to_underlying(tv);
            }

            onTrack = true;
            varsOnTrack.clear();
        }

        /**
         * Checks if still on track to last solution and sets the internal on track flag to false if this is not the
         * case. Also returns a list with variables in conflict with the previous solution along with their
         * corresponding truth values
         * @tparam TF truth function
         * @param truthFunction function object that maps variable ids to truth values
         * For example: undefined aligns with true and false but false only aligns with false.
         * @return Vector of pairs of conflicting variables along with their corresponding truth value from the last
         * solution
         */
        template<concepts::callable_r<TruthVal, var_t> TF>
        auto updateOnTrack(TF &&truthFunction) -> Conflicts {
            Conflicts ret;
            varsOnTrack.clear();
            for (auto [var, val] : iterators::const_enumerate(varPolarity, offset)) {
                if (isEqual(std::forward<TF>(truthFunction)(var), val)) {
                    varsOnTrack.emplace_back(var, val);
                }

                if (not isAligned(std::forward<TF>(truthFunction)(var), val)) {
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
            if (onTrack) {
                varsOnTrack.emplace_back(decision.variable(), decision.sign());
            }
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

        /**
         * Gets the variable id offset. This value is necessary when iterating a solution because a variable x does
         * not correspond to a polarity value at index x but rather x - offset
         * @return variable id offset
         */
        [[nodiscard]] var_t getOffset() const noexcept;

        /**
         * Gets the decided variables up to this point that are contained in the last solution (that are on track).
         * @return All set variables that have the same polarity as in the last solutions
         * @note the contents are not updated when isOnTrack() evaluates to false. Also note that contained variables
         * are might have been fixed by propagation or pruning
         */
        [[nodiscard]] auto getVariablesOnTrack() const noexcept -> const serialization::Branch &;

    protected:
        static bool isAligned(TruthVal tv, bool polarity) noexcept {
            return tv == TruthVal::Undefined or to_underlying(tv) == polarity;
        }

        static bool isEqual(TruthVal tv, bool polarity) noexcept {
            return tv != TruthVal::Undefined and to_underlying(tv) == polarity;
        }
    };

    /**
     * @brief Functor object that gets the truth value of a binary variable from the solver
     * @tparam T timing type
     */
    template<concepts::scalar T>
    class TruthFunction {
        const Solver<T> &s;
    public:
        explicit constexpr TruthFunction(const Solver<T> &solver) noexcept: s(solver) {}

        constexpr auto operator()(var_t var) const noexcept -> TraceWatcher::TruthVal {
            if (s.boolean.isUndefined(var)) {
                return TraceWatcher::TruthVal::Undefined;
            }

            return static_cast<TraceWatcher::TruthVal>(s.boolean.isTrue(var));
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
        SubscriberHandle propCompletedHandler;
    public:
        SubscribableEvent<DeviationType, const TraceWatcher::Conflicts &,
                const serialization::Branch &> DeviationOccurred;
        ///< triggered when the solver deviates from path to the last solution. Arguments: deviation type,
        ///< conflicts after propagation, fixed variables on track

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
        decisionHandler(solver.ChoicePoint.subscribe_handled([this](auto &&, auto lit) {
            this->handleDecision(lit);
        })),
        conflictHandler(solver.ConflictEncountered.subscribe_handled([this](const auto &) {
            this->handleConflict();
        })),
        backtrackHandler(solver.BackTrackCompleted.subscribe_handled([this](const auto &solver) {
            this->watcher.updateOnTrack(TruthFunction(solver));
        })),
        propCompletedHandler(solver.PropagationCompleted.subscribe_handled([this](const auto &solver) {
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

            auto variablesOnTrack = watcher.getVariablesOnTrack();
            auto conflicts = watcher.updateOnTrack(TruthFunction(solver));
            if (not watcher.isOnTrack()) {
                DeviationOccurred.trigger(DeviationType::Propagation, conflicts, variablesOnTrack);
            }
        }

        template<concepts::scalar T>
        void handleSolution(const Solver<T> &solver) {
            watcher.registerSolution(TruthFunction(solver));
        }

        template<concepts::scalar T>
        void handleDecision(Literal<T> lit) {
            bool ot = watcher.isOnTrack();
            watcher.step(lit);
            if (ot != watcher.isOnTrack()) {
                DeviationOccurred.trigger(DeviationType::Decision, TraceWatcher::Conflicts{},
                                          watcher.getVariablesOnTrack());
            }
        }

        void handleConflict();
    };
}

#endif //TEMPO_TRACEWATCHER_HPP
