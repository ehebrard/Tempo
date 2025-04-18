/**
* @author Tim Luchterhand
* @date 09.07.24
* @brief
*/

#ifndef TEMPO_HEURISTIC_INTERFACE_HPP
#define TEMPO_HEURISTIC_INTERFACE_HPP

#include <concepts>
#include <memory>


#include "Literal.hpp"
#include "util/traits.hpp"
#include "util/PolyFactory.hpp"
#include "util/Options.hpp"

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

    /**
     * Interface concept for branching heuristics (combined variable and value selection)
     * @tparam H heuristic type
     * @tparam T timing type
     */
    template<typename H, typename T>
    concept heuristic = requires(H heuristic, const Solver<T> solver) {
        { heuristic.branch(solver) } -> std::same_as<Literal<T>>;
    };

    /**
     * Interface concept for variable selection heuristics
     * @tparam H heuristic type
     * @tparam T timing type
     */
    template<typename H, typename T>
    concept variable_heuristic = requires(H heuristic, const Solver<T> solver) {
        { heuristic.nextVariable(solver) } -> std::same_as<VariableSelection>;
    };

    /**
     * Interface concept for value selection heuristics
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
    class MovableHeuristic : public std::unique_ptr<H> {
    public:
        using std::unique_ptr<H>::unique_ptr;
        /**
         * Ctor. Constructs the heuristic in place
         * @tparam Args argument types
         * @param args arguments to the Ctor of H
         */
        template<typename... Args>
        explicit MovableHeuristic(Args &&... args) : std::unique_ptr<H>(
            std::make_unique<H>(std::forward<Args>(args)...)) {}


        /**
         * Variable heuristic interface
         * @tparam T timing type
         * @param solver
         * @return variable selection
         */
        template<concepts::scalar T> requires(variable_heuristic<H, T>)
        auto nextVariable(const Solver<T> &solver) const -> VariableSelection {
            return this->get()->nextVariable(solver);
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
            return this->get()->valueDecision(x, solver);
        }

        /**
         * Heuristic interface
         * @tparam T timing type
         * @param solver
         * @return branching literal
         */
        template<concepts::scalar T> requires(tempo::heuristics::heuristic<H, T>)
        auto branch(const Solver<T> &solver) const -> Literal<T> {
            return this->get()->branch(solver);
        }
    };


    /**
     * @brief Abstract variable selection heuristic
     * @details @copybrief
     * @tparam T timing type
     * @note You don't have to use this base class when implementing new variable selection heuristics. They only have
     * to model the concepts further above. Use this base class when you want to create your heuristic using a factory
     */
    template<concepts::scalar T>
    class BaseVariableHeuristic {
    public:
        BaseVariableHeuristic() = default;
        BaseVariableHeuristic(const BaseVariableHeuristic &) = default;
        BaseVariableHeuristic(BaseVariableHeuristic &&) = default;
        BaseVariableHeuristic &operator=(const BaseVariableHeuristic &) = default;
        BaseVariableHeuristic &operator=(BaseVariableHeuristic &&) = default;
        virtual ~BaseVariableHeuristic() = default;

        /**
         * variable selection interface
         * @param solver solver for which to select a new variable
         * @return variable selection
         */
        virtual auto nextVariable(const Solver<T> &solver) -> VariableSelection = 0;
    };

    /**
     * @brief Polymorphic variable selection heuristic type returned from factories
     * @tparam T timing type
     */
    template<concepts::scalar T>
    using VariableHeuristic = MovableHeuristic<BaseVariableHeuristic<T>>;

    /**
     * @brief These are all the scalar types supported by the heuristic factories.
     * If you need more, just extend the list.
     */
    using ScalarTypes = types_t<int, unsigned, long, unsigned long, float, double>;

    /**
     * @brief Helper CRTP type for creating variable selection heuristic factories (have your factory derive from this)
     * @tparam Impl Class implementing the CRTP
     */
    template<typename Impl>
    using MakeVariableHeuristicFactory = MakeFactory<Impl, VariableHeuristic, ScalarTypes, RefTag<Solver>>;

    /**
     * @brief Singleton variable selection heuristic factory type. Use this type to register your factories
     */
    using VariableHeuristicFactory = PolyFactory<VariableHeuristic, ScalarTypes, Options::ChoicePointHeuristics,
        RefTag<Solver>>;


    /**
     * @brief abstract value selection heuristic
     * @details @copybrief
     * @tparam T timing type
     * @note You don't have to use this base class when implementing new value selection heuristics. They only have
     * to model the concepts further above. Use this base class when you want to create your heuristic using a factory
     */
    template<concepts::scalar T>
    class BaseValueHeuristic {
    public:
        BaseValueHeuristic() = default;
        BaseValueHeuristic(const BaseValueHeuristic &) = default;
        BaseValueHeuristic(BaseValueHeuristic &&) = default;
        BaseValueHeuristic &operator=(const BaseValueHeuristic &) = default;
        BaseValueHeuristic &operator=(BaseValueHeuristic &&) = default;
        virtual ~BaseValueHeuristic() = default;

        /**
         * value selection interface
         * @param x selected variable
         * @param solver solver for which to select the variable value
         * @return branching literal
         */
        virtual auto valueDecision(const VariableSelection &x, const Solver<T> &solver) -> Literal<T> = 0;
    };

    /**
     * @brief Polymorphic value selection heuristic type returned from factories
     * @tparam T timing type
     */
    template<concepts::scalar T>
    using ValueHeuristic = MovableHeuristic<BaseValueHeuristic<T>>;

    /**
     * @brief Helper CRTP type for creating value selection heuristic factories (have your factory derive from this)
     * @tparam Impl Class implementing the CRTP
     */
    template<typename Impl>
    using MakeValueHeuristicFactory = MakeFactory<Impl, ValueHeuristic, ScalarTypes, RefTag<Solver>>;

    /**
     * @brief Singleton value selection heuristic factory type. Use this type to register your factories
     */
    using ValueHeuristicFactory = PolyFactory<ValueHeuristic, ScalarTypes, Options::PolarityHeuristic, RefTag<Solver>>;


    namespace detail {
        template<concepts::scalar T>
        struct PolyHeuristicCallableBase {
            PolyHeuristicCallableBase() = default;
            virtual ~PolyHeuristicCallableBase() = default;
            PolyHeuristicCallableBase(PolyHeuristicCallableBase&&) = default;
            PolyHeuristicCallableBase &operator=(PolyHeuristicCallableBase&&) = default;
            virtual auto invoke(const Solver<T> &) -> Literal<T> = 0;
        };

        template<concepts::scalar T, heuristic<T> H>
        struct PolyHeuristicCallable : PolyHeuristicCallableBase<T> {
            H impl;
            template<typename ...Args>
            explicit PolyHeuristicCallable(Args &&...args) : impl(std::forward<Args>(args)...) {}
            auto invoke(const Solver<T> &solver) -> Literal<T> override {
                return impl.branch(solver);
            }
        };
    }

    /**
     * @brief Polymorphic wrapper for arbitrary branching heuristics
     * @details @copybrief
     * (combined variable selection + value selection)
     * @tparam T timing type
     * @note use std::movable_function from c++23 instead of detail::PolyHeuristicCallable implementation in the future
     */
    template<concepts::scalar T>
    class PolymorphicHeuristic {
        std::unique_ptr<detail::PolyHeuristicCallableBase<T>> impl;
    public:
        /**
         * default ctor
         * Creates an invalid wrapper that must not be called
         */
        constexpr PolymorphicHeuristic() noexcept: impl(nullptr) {}

        /**
         * Ctor
         * @tparam H heuristic type
         */
        template<heuristic<T> H>
        PolymorphicHeuristic(H &&heuristic) :impl(
                std::make_unique<detail::PolyHeuristicCallable<T, std::decay_t<H>>>(std::forward<H>(heuristic))) {}

        /**
         * Branching heuristic interface
         * @param solver solver for which to select a branching literal
         * @return selected branching literal
         * @throws std::runtime_error if wrapper does not contain valid heuristic
         */
        auto branch(const Solver<T> &solver) -> Literal<T> {
            if (nullptr == impl) {
                throw std::runtime_error("no branching heuristic set");
            }

            return impl->invoke(solver);
        }

        /**
         * Whether the wrapper holds a valid heuristic
         */
        [[nodiscard]] bool isValid() const noexcept {
            return nullptr != impl;
        }
    };

    /**
     * @brief General branching heuristic wrapper that holds a separate variable and value selection heuristic
     * @details @copybrief
     * @tparam VariableHeuristic type of variable selection heuristic
     * @tparam ValueHeuristic type of value selection heuristic
     */
    template<typename VariableHeuristic, typename ValueHeuristic>
    class CompoundHeuristic {
        VariableHeuristic variableSelector;
        ValueHeuristic valueSelector;
    public:

        /**
         * Ctor
         * @tparam VarH type of variable heuristic
         * @tparam ValH type of value heuristic
         * @param variableHeuristic the variable heuristic
         * @param valueHeuristic the value heuristic
         */
        template<typename VarH, typename ValH>
        CompoundHeuristic(VarH &&variableHeuristic, ValH &&valueHeuristic) :
                variableSelector(std::forward<VarH>(variableHeuristic)),
                valueSelector(std::forward<ValH>(valueHeuristic)) {}

        template<concepts::scalar T>
        requires(variable_heuristic<VariableHeuristic, T> and value_heuristic<ValueHeuristic, Solver<T>>)
        auto branch(const Solver<T> &solver) -> Literal<T> {
            return valueSelector.valueDecision(variableSelector.nextVariable(solver), solver);
        }
    };

    /**
     * Helper function to create compound heuristics
     * @tparam VarH type of variable heuristic
     * @tparam ValH type of value heuristic
     * @param variableHeuristic the variable heuristic
     * @param valueHeuristic the value heuristic
     */
    template<typename VarH, typename ValH>
    auto make_compound_heuristic(VarH &&variableHeuristic, ValH &&valueHeuristic) {
        return CompoundHeuristic<VarH, ValH>(std::forward<VarH>(variableHeuristic), std::forward<ValH>(valueHeuristic));
    }
}

#endif //TEMPO_HEURISTIC_INTERFACE_HPP
