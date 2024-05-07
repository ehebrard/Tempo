/**
* @author Tim Luchterhand
* @date 06.05.24
* @brief Tightest value selection
*/

#ifndef TEMPO_TIGHTESTVALUE_HPP
#define TEMPO_TIGHTESTVALUE_HPP

#include "BaseValueHeuristic.hpp"
#include "Global.hpp"
#include "Constant.hpp"
#include "util/factory_pattern.hpp"
#include "util/traits.hpp"

namespace tempo {
    template<typename T>
    class Scheduler;
}

namespace tempo::heuristics {

    /**
     * @brief Tightest value selection heuristic.
     * @details @copybrief
     * Chooses the polarity that would leave the most slack in the timing network
     */
    class TightestValue : public BaseValueHeuristic<TightestValue>{
    public:
        /**
         * Ctor
         * @param epsilon see tempo::heuristics::BaseValueHeuristic
         */
        explicit TightestValue(double epsilon): BaseValueHeuristic<TightestValue>(epsilon) {}

        /**
         * heuristic interface
         * @tparam T timing type
         * @param cp choice point
         * @param scheduler scheduler instance
         * @return either POS(cp) or NEG(cp)
         */
        template<concepts::scalar T>
        static lit choose(var cp, const Scheduler<T> &scheduler) {
            lit d = NoVar;
            auto prec_a{scheduler.getEdge(POS(cp))};
            auto prec_b{scheduler.getEdge(NEG(cp))};

            T gap_a{0};
            if (prec_a != Constant::NoEdge<T>) {
                gap_a = scheduler.upper(prec_a.from) - scheduler.lower(prec_a.to);
            }

            T gap_b{0};
            if (prec_b != Constant::NoEdge<T>) {
                gap_b = scheduler.upper(prec_b.from) - scheduler.lower(prec_b.to);
            }
            d = (gap_a < gap_b ? NEG(cp) : POS(cp));

#ifdef DBG_TRACE
            if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
                std::cout << scheduler << "\n-- new decision: " << prettyLiteral(EDGE(d))
                          << std::endl;
            }
#endif

            return d;
        }
    };

    MAKE_FACTORY(TightestValue, const ValueHeuristicConfig &config) {
            return TightestValue(config.epsilon);
    }};
}


#endif //TEMPO_TIGHTESTVALUE_HPP
