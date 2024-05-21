//
// Created by tim on 15.11.22.
//

#ifndef TEMPO_VSIDS_HPP
#define TEMPO_VSIDS_HPP
//#include <vector>
//#include <optional>
#include "BaseHeuristic.hpp"
//#include "util/traits.hpp"
#include "heuristics/impl/DecayingEventActivityMap.hpp"
#include "util/SubscribableEvent.hpp"

namespace tempo::heuristics {
    /**
     * @brief VSIDS (Variable State Independent Decaying Sum) heuristic
     * @details @copybrief
     * Selects variables which occur in the largest number of clauses (proportionally to highest activity). Activity for
     * a literal l is calculated the following way:
     * \f$ A_l = \sum_{i=1}^n \gamma^{n - i} w_i \f$ where \f$ \gamma \f$ is a constant decay factor and \f$ n \f$ is
     * the the number of calls to the function VSIDS::updateActivity
     * The final score of a literal is inversely proportional to its activity (@see getCost)
     */
    template<typename T>
    class VSIDS : public BaseHeuristic<VSIDS<T>> {
    public:

        /**
         * CTor. Infers parameters from given scheduler. Subscribes to events of interest
         * @tparam T
         * @param scheduler
         * @param options command line options
         */
        explicit VSIDS(Scheduler<T> &scheduler) :
//                sched(scheduler),
                activity(scheduler, scheduler.getOptions().vsids_decay),
                handlerToken(scheduler.ClauseAdded.subscribe_handled(
                        [this,&scheduler](const auto &arg) { this->activity.update(arg, scheduler); })) {}
        
        
        explicit VSIDS(Solver<T> &solver) :
//                solver(s),
                activity(solver, solver.getOptions().vsids_decay),
                handlerToken(solver.ClauseAdded.subscribe_handled(
                        [this,&solver](const auto &arg) { this->activity.update(arg, solver); })) {}

        // Copy and move are disabled. Otherwise, a call to the subscribed event handler will cause undefined behavior
        VSIDS(const VSIDS &) = delete;
        VSIDS(VSIDS &&) = delete;
        VSIDS &operator=(const VSIDS&) = delete;
        VSIDS &operator=(VSIDS &&) = delete;
        ~VSIDS() = default;

        /**
         * Cost for a choice point
         * @tparam T type of DistanceConstraint
         * @param choicePoint choice point to evaluate
         * @return maximum of the distance between from and to in both directions in the temporal network divided by
         * the choice points activity
         */
        [[nodiscard]] double getCost(const var x, const Scheduler<T>& sched) const {
          auto prec_a{sched.getEdge(POS(x))};
          auto prec_b{sched.getEdge(NEG(x))};
          auto gap_a = sched.upper(prec_a.from) - sched.lower(prec_a.to);
          auto gap_b = sched.upper(prec_b.from) - sched.lower(prec_b.to);
          return static_cast<double>(std::max(gap_a, gap_b)) / activity.get(x, sched);
        }
        
        [[nodiscard]] double getCost(const var_t x, const Solver<T>& solver) const {
            
//                std::cout << "var = " << x << "\n";
            
            double dom{1};
            
            if(solver.boolean.hasSemantic(x)) {
                auto p{solver.boolean.getLiteral(true,x)};
                auto n{solver.boolean.getLiteral(false,x)};
                
                //            std::cout << p << " <> " << n << std::endl;
                
                auto prec_a{solver.boolean.getEdge(p)};
                auto prec_b{solver.boolean.getEdge(n)};
                
                
                //            std::cout << prec_a << " <> " << prec_b << std::endl;
                
                auto gap_a = solver.numeric.upper(prec_a.from) - solver.numeric.lower(prec_a.to);
                auto gap_b = solver.numeric.upper(prec_b.from) - solver.numeric.lower(prec_b.to);
                
                dom = static_cast<double>(std::max(gap_a, gap_b));
            }
            
//            std::cout << "var = " << x << " act = " << activity.get(x, solver) << "\n";
            
          return dom / activity.get(x, solver);
        }

    private:
//      const Scheduler<T> &sched;
      impl::DecayingEventActivityMap<T> activity;
      SubscriberHandle handlerToken;
    };
}

#endif //TEMPO_VSIDS_HPP
