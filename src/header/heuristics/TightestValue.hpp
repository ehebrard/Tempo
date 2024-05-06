/**
* @author Tim Luchterhand
* @date 06.05.24
* @brief
*/

#ifndef TEMPO_TIGHTESTVALUE_HPP
#define TEMPO_TIGHTESTVALUE_HPP

#include "BaseValueHeuristic.hpp"
#include "Global.hpp"
#include "Constant.hpp"
#include "util/factory_pattern.hpp"
#include "util/traits.hpp"
#include "ValueHeuristicConfig.hpp"

namespace tempo {
    template<typename T>
    class Scheduler;
}

namespace tempo::heuristics {
    class TightestValue : public BaseValueHeuristic<TightestValue>{
    public:
        explicit TightestValue(const ValueHeuristicConfig &config): BaseValueHeuristic<TightestValue>(config.epsilon) {}

        template<concepts::scalar T>
        [[nodiscard]] lit choose(var cp, const Scheduler<T> &scheduler) const {
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

            //      double g{1};
            ////          std::cout << static_cast<double>(gap_a) << " <> " <<
            /// static_cast<double>(gap_b) << ": " << gap_ratio << " -> ";
            //          if(gap_a < gap_b) {
            ////              std::cout << "(" <<
            ///(static_cast<double>(gap_a)/static_cast<double>(gap_b)) << ") ";
            //
            //              g = (static_cast<double>(gap_a)/static_cast<double>(gap_b));
            //          } else {
            ////              std::cout << "(" <<
            ///(static_cast<double>(gap_b)/static_cast<double>(gap_a)) << ") ";
            //
            //              g = (gap_a > 0 ?
            //              static_cast<double>(gap_b)/static_cast<double>(gap_a) : 1);
            //          }
            ////          std::cout << gap_ratio << "\n";
            //
            //      gap_ratio += g;
            //      if(g == 1) {
            //          ++num_tight;
            //      }

            d = (gap_a < gap_b ? NEG(cp) : POS(cp));
#ifdef DBG_TRACE
            if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
    std::cout << *this << "\n-- new decision: " << prettyLiteral(EDGE(d))
              << std::endl;
  }
#endif

            return d;
        }
    };

    MAKE_FACTORY(TightestValue, const ValueHeuristicConfig &config) {
            return TightestValue(config);
    }};
}


#endif //TEMPO_TIGHTESTVALUE_HPP
