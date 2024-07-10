/**
* @author Tim Luchterhand
* @date 09.07.24
* @brief
*/

#ifndef TEMPO_HEURISTIC_INTERFACE_HPP
#define TEMPO_HEURISTIC_INTERFACE_HPP

#include <concepts>
#include <memory>
#include <variant>
#include <functional>

#include "Literal.hpp"
#include "util/traits.hpp"

namespace tempo {
    template<typename T>
    class Solver;
}

namespace tempo::heuristics {
    enum class VariableType {
        Boolean,
        Numeric
    };

    using VariableSelection = std::pair<var_t, VariableType>;

    template<typename H, typename T>
    concept heuristic = requires(H heuristic, const Solver<T> solver) {
        { heuristic.branch(solver) } -> std::same_as<Literal<T>>;
    };

    template<typename H, typename T>
    concept variable_heuristic = requires(H heuristic, const Solver<T> solver) {
        { heuristic.nextVariable(solver) } -> std::same_as<VariableSelection>;
    };

    /**
     * Interface for value selection heuristics
     * @tparam H heuristic type
     * @tparam Solver information provider (usually the scheduler)
     */
    template <typename H, typename Solver>
    concept value_heuristic = requires(H heuristic, VariableSelection x, const Solver &solver) {
        { heuristic.valueDecision(x, solver) } -> concepts::same_template<Literal>;
    };

    /**
     * Wrapper class for heuristics that do not support moving
     * @tparam H type of heuristic
     */
    template<typename H>
    class MovableHeuristic {
        std::unique_ptr<H> heuristic;
    public:
        /**
         * Ctor. Constructs the heuristic in place
         * @tparam Args argument types
         * @param args arguments to the Ctor of H
         */
        template<typename ...Args>
        explicit MovableHeuristic(Args &&...args) : heuristic(std::make_unique<H>(std::forward<Args>(args)...)) {}


        /**
         * Variable heuristic interface
         * @tparam T timing type
         * @param solver
         * @return variable selection
         */
        template<concepts::scalar T> requires(variable_heuristic<H, T>)
        auto nextVariable(const Solver<T> &solver) const -> VariableSelection {
            return heuristic->nextVariable(solver);
        }

        /**
         * Value heuristic interface
         * @tparam T timing type
         * @param x variable selection
         * @param solver
         * @return branching literal
         */
        template<concepts::scalar T> requires(value_heuristic<H, Solver<T>>)
        auto valueDecision(VariableSelection x, const Solver<T> &solver) const -> Literal<T> {
            return heuristic->valueDecision(x, solver);
        }

        /**
         * Heuristic interface
         * @tparam T timing type
         * @param solver
         * @return branching literal
         */
        template<concepts::scalar T> requires(tempo::heuristics::heuristic<H, T>)
        auto branch(const Solver<T> &solver) const -> Literal<T> {
            return heuristic->branch(solver);
        }
    };

    template<concepts::scalar T>
    class PolymorphicVariableHeuristic {
        std::function<Literal<T>(const Solver<T> &)> run;
    public:
        template<variable_heuristic<T> H>
        PolymorphicVariableHeuristic(H &&heuristic) : run(
                [h_ = std::forward<H>(heuristic)](auto &s) { return h_.nextVariable(s); }) {}

        auto variableDecision(const Solver<T> &solver) {
            return run(solver);
        }
    };

    template<typename ...Heuristics>
    struct PolymorphicHeuristic : std::variant<Heuristics...> {
        using std::variant<Heuristics...>::variant;
        using std::variant<Heuristics...>::emplace;

        template<concepts::scalar T> requires(variable_heuristic<Heuristics, T> && ...)
        auto nextVariable(const Solver<T> &solver) {
            return std::visit([&solver](auto &h) {return h.nextVariable(solver);}, *this);
        }

        template<concepts::scalar T> requires(value_heuristic<Heuristics, Solver<T>> && ...)
        auto valueDecision(VariableSelection x, const Solver<T> &solver) {
            return std::visit([x, &solver](auto &h) {return h.valueDecision(x, solver);}, *this);
        }

        template<concepts::scalar T> requires(heuristic<Heuristics, T> && ...)
        auto branch(const Solver<T> &solver) {
            return std::visit([&solver](auto &h) {return h.branch(solver);}, *this);
        }
    };


    template<typename VarH, typename ValH>
    class CompoundHeuristic {
        VarH variableSelector;
        ValH valueSelector;
    public:
        CompoundHeuristic(VarH variableHeuristic, ValH valueHeuristic) : variableSelector(std::move(variableHeuristic)),
                                                                         valueSelector(std::move(valueHeuristic)) {}

        template<concepts::scalar T>
        requires(variable_heuristic<VarH, T> and value_heuristic<ValH, Solver<T>>)
        auto branch(const Solver<T> &solver) -> Literal<T> {
            return valueSelector.valueDecision(variableSelector.nextVariable(solver), solver);
        }
    };

}

#endif //TEMPO_HEURISTIC_INTERFACE_HPP
