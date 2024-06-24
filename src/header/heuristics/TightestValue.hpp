/**
 * @author Tim Luchterhand
 * @date 06.05.24
 * @brief Tightest value selection
 */

#ifndef TEMPO_TIGHTESTVALUE_HPP
#define TEMPO_TIGHTESTVALUE_HPP

#include "BaseBooleanHeuristic.hpp"
#include "Constant.hpp"
#include "DistanceConstraint.hpp"
#include "Global.hpp"
#include "Literal.hpp"
#include "util/factory_pattern.hpp"
#include "util/traits.hpp"


namespace tempo::heuristics {

    namespace detail {
        template<typename Solver, typename T>
        concept distance_provider = requires(const Solver s, Literal<T> l, event e) {
            { s.boolean.getEdge(l) } -> concepts::same_template<DistanceConstraint>;
            { s.numeric.upper(e) } -> concepts::scalar;
            { s.numeric.lower(e) } -> concepts::scalar;
        };
    }

/**
 * @brief Tightest value selection heuristic.
 * @details @copybrief
 * Chooses the polarity that would leave the most slack in the timing network
 */
    class TightestValue : public BaseBooleanHeuristic<TightestValue> {
    public:
        /**
         * Ctor
         * @param epsilon see tempo::heuristics::BaseValueHeuristic
         */
        explicit TightestValue(double epsilon)
                : BaseBooleanHeuristic<TightestValue>(epsilon) {}

        /**
         * heuristic interface
         * @tparam Solver class that provides distances between events and a mapping
         * from literals to edges
         * @param cp choice point
         * @param solver scheduler instance
         * @return either POS(cp) or NEG(cp)
         */
        template<concepts::scalar T, typename Solver>
        static bool choose(Literal<T> lit, const Solver &solver) {
            static_assert(detail::distance_provider<Solver, T>);
            // @TODO no gap info available -> what should I return?
            if (not lit.hasSemantic()) {
                return true;
            }

            auto edgePos = solver.boolean.getEdge(lit);
            auto edgeNeg = solver.boolean.getEdge(~lit);
            auto gapPos = solver.numeric.upper(edgePos.from) - solver.numeric.lower(edgePos.to);
            auto gapNeg = solver.numeric.upper(edgeNeg.from) - solver.numeric.lower(edgeNeg.to);
            return gapPos <= gapNeg;
        }
    };

    MAKE_FACTORY(TightestValue, const ValueHeuristicConfig &config) {
            return TightestValue(config.epsilon);
        }
    };
}

#endif // TEMPO_TIGHTESTVALUE_HPP
