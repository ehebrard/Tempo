/**
* @author Tim Luchterhand
* @date 11.03.25
* @file ValueImpact.hpp
* @brief generic value heuristic for binary variables based on impact measure
*/

#ifndef VALUEIMPACT_HPP
#define VALUEIMPACT_HPP

#include <stdexcept>

#include "util/traits.hpp"
#include "impl/ImpactMap.hpp"
#include "heuristic_interface.hpp"
#include "BaseBooleanHeuristic.hpp"
#include "TightestValue.hpp"

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

        /**
         * Get the percentage of literals for which impact information is available
         * @return percentage in [0, 1]
         */
        [[nodiscard]] double getInformationDensity() const {
            double numSet = 0;
            for (auto i : impact.getMap()) {
                numSet += i != 0;
            }

            return numSet / static_cast<double>(impact.getMap().size());
        }
    };

    /**
     * @brief Value heuristic that switches form tightest to impact based as soon as enough impact information is
     * available
     * @details @copybrief
     * @tparam T Timing type
     */
    template<concepts::scalar T>
    class TightestImpact: public BaseValueHeuristic<T> {
        TightestValue<T> tightest;
        ValueImpact<T> impact;
        double impactPercentage;
        bool warmedUp = false;

    public:
        /**
         * Ctor
         * @param solver target solver
         * @param impactPercentage percentage threshold at which the heuristic switches from tightest to impact based
         * @param epsilon epsilon value
         */
        TightestImpact(const Solver<T> &solver, double impactPercentage, double epsilon = 0)
            : tightest(epsilon), impact(solver, epsilon), impactPercentage(impactPercentage) {
            if (impactPercentage < 0 or impactPercentage > 1) {
                throw std::invalid_argument("impactPercentage must be between 0 and 1");
            }
        }

        /**
         * value heuristic interface
         * @param x varaible selection
         * @param solver target solver
         * @return literal chosen by tightest if not enough impact information, literal chosen by impact otherwise
         */
        auto valueDecision(const VariableSelection &x, const Solver<T> &solver) -> Literal<T> override {
            if (not warmedUp and impact.getInformationDensity() >= impactPercentage) {
                if (solver.getOptions().verbosity >= Options::YACKING) {
                    std::cout << "-- Tightest + Impact value heuristic now switching to impact only" << std::endl;
                }
                warmedUp = true;
            }

            if (warmedUp) {
                return impact.valueDecision(x, solver);
            }

            return tightest.valueDecision(x, solver);
        }
    };

    struct ValueImpactFactory : MakeValueHeuristicFactory<ValueImpactFactory> {
        ValueImpactFactory();

        template<concepts::scalar T>
        [[nodiscard]] auto build_impl(const Solver<T> &solver) const -> ValueHeuristic<T> {
            return std::make_unique<ValueImpact<T>>(solver, solver.getOptions().polarity_epsilon);
        }
    };

    struct TightestImpactFactory : MakeValueHeuristicFactory<TightestImpactFactory> {
        TightestImpactFactory();

        template<concepts::scalar T>
        [[nodiscard]] auto build_impl(const Solver<T> &solver) const -> ValueHeuristic<T> {
            return std::make_unique<TightestImpact<T>>(solver, 0.5, solver.getOptions().polarity_epsilon);
        }
    };
}

#endif //VALUEIMPACT_HPP
