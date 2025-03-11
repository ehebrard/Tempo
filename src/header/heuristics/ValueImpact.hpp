/**
* @author Tim Luchterhand
* @date 11.03.25
* @file ValueImpact.hpp
* @brief generic value heuristic for binary variables based on impact measure
*/

#ifndef VALUEIMPACT_HPP
#define VALUEIMPACT_HPP

#include "util/traits.hpp"
#include "impl/ImpactMap.hpp"
#include "heuristic_interface.hpp"
#include "BaseBooleanHeuristic.hpp"

namespace tempo::heuristics {
    /**
     * @brief Generic value heuristic for binary variables based on a measure of impact
     * @tparam T timing type
     */
    template<concepts::scalar T>
    class ValueImpact : public BaseBooleanHeuristic<ValueImpact<T>, T> {
        impl::ImpactMap impact;
    public:
        /**
         * Ctor
         * @param solver target solver
         * @param epsilon epsilon value
         */
        explicit ValueImpact(const Solver<T> &solver, double epsilon = 0)
            : BaseBooleanHeuristic<ValueImpact, T>(epsilon), impact(solver) {}

        /**
         * Heuristic interface
         * @param x variable selection
         * @param solver target solver
         * @return literal with minimum impact
         */
        auto choose(var_t x, const Solver<T> &solver) -> Literal<T> {
            auto pos = solver.boolean.getLiteral(true, x);
            auto neg = solver.boolean.getLiteral(false, x);
            return impact.get(pos) > impact.get(neg) ? neg : pos;
        }
    };

    struct ValueImpactFactory : MakeValueHeuristicFactory<ValueImpactFactory> {
        ValueImpactFactory();

        template<concepts::scalar T>
        [[nodiscard]] auto build_impl(const Solver<T> &solver) const -> ValueHeuristic<T> {
            return std::make_unique<ValueImpact<T>>(solver, solver.getOptions().polarity_epsilon);
        }
    };
}

#endif //VALUEIMPACT_HPP
