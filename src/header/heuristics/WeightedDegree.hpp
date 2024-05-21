//
// Created by tluchterha on 22/11/22.
//

#ifndef TEMPO_WEIGHTEDDEGREE_HPP
#define TEMPO_WEIGHTEDDEGREE_HPP
#include "util/Options.hpp"
#include "util/SubscribableEvent.hpp"
#include "heuristics/BaseHeuristic.hpp"
#include "heuristics/impl/DecayingEventActivityMap.hpp"
//#include <array>

namespace tempo::heuristics {
/**
 * @brief Weighted degree heuristic. Prioritized choice points according to
 * their recent activity
 * @details @copybrief
 * Registers the activity of literals. The activity is a decaying counter which
 * counts the number of times the literal has been involved in a conflict.
 * @note The difference to VSIDS is that the activity of a literal depends on
 * how often it is involved in a conflict and not how often it is contained in a
 * learned clause
 */
template <typename T>
class WeightedDegree : public BaseHeuristic<WeightedDegree<T>> {
public:
  /**
   * @copydoc VSIDS::VSIDS
   */
//  explicit WeightedDegree(Scheduler<T> &scheduler, const bool critpath)
//      : sched(scheduler), activity(scheduler, sched.getOptions().vsids_decay),
//        handlerToken((critpath
//                          ? scheduler.ConflictEncountered.subscribe_handled(
//                                [this](auto &expl) {
//                                  this->clause.clear();
//                                  if (expl.getType() == CYCLEEXPL)
//                                    expl.explain(NoLit, this->clause);
//                                  else
//                                    sched.getCriticalPath(this->clause);
//                                  this->activity.update(this->clause);
//                                })
//                          : scheduler.ConflictEncountered.subscribe_handled(
//                                [this](auto &expl) {
//                                  this->clause.clear();
//                                  expl.explain(NoLit, this->clause);
//                                  this->activity.update(this->clause);
//                                })))
//  //                                                                     :
//  //                                [this](auto &expl) {
//  //                this->clause.clear();
//  //                expl.explain(NoLit,this->clause);
//  //              this->activity.update(this->clause);
//  //            })
//  //                     )
//  //
//  //
//  //                     )
//  {}
    
    explicit WeightedDegree(Scheduler<T> &scheduler)
        : activity(scheduler, scheduler.getOptions().vsids_decay),
          handlerToken(scheduler.ConflictEncountered.subscribe_handled(
                                  [this,&scheduler](auto &expl) {
                                    this->clause.clear();
                                    expl.explain(NoLit, this->clause);
                                    this->activity.update(this->clause, scheduler);
                                  }))
    {}
    
    explicit WeightedDegree(Solver<T> &solver)
        : activity(solver, solver.getOptions().vsids_decay),
          handlerToken(solver.ConflictEncountered.subscribe_handled(
                                  [this,&solver](auto &expl) {
                                    this->nclause.clear();
                                    expl.explain(Solver<T>::Contradiction, this->nclause);
                                    this->activity.update(this->nclause, solver);
                                  }))
    {}

  /**
   * @copydoc VSIDS::getCost
   */
  [[nodiscard]] double getCost(const var x, const Scheduler<T>& sched) const {
    auto prec_a{sched.getEdge(POS(x))};
    auto prec_b{sched.getEdge(NEG(x))};
    auto gap_a = sched.upper(prec_a.from) - sched.lower(prec_a.to);
    auto gap_b = sched.upper(prec_b.from) - sched.lower(prec_b.to);
    return static_cast<double>(std::max(gap_a, gap_b)) / activity.get(x, sched);
  }
    
    [[nodiscard]] double getCost(const var x, const Solver<T>& solver) const {
        double dom{1};
        
        if(solver.boolean.hasSemantic(x)) {
            auto p{solver.boolean.getLiteral(true,x)};
            auto n{solver.boolean.getLiteral(false,x)};
            
            auto prec_a{solver.boolean.getEdge(p)};
            auto prec_b{solver.boolean.getEdge(n)};
            
            auto gap_a = solver.numeric.upper(prec_a.from) - solver.numeric.lower(prec_a.to);
            auto gap_b = solver.numeric.upper(prec_b.from) - solver.numeric.lower(prec_b.to);
            
            dom = static_cast<double>(std::max(gap_a, gap_b));
        }
        
        return dom / activity.get(x, solver);
    }

  WeightedDegree(const WeightedDegree &) = delete;
  WeightedDegree(WeightedDegree &&) = delete;
  WeightedDegree &operator=(const WeightedDegree &) = delete;
  WeightedDegree &operator=(WeightedDegree &&) = delete;
  ~WeightedDegree() = default;

private:
//  Scheduler<T> &sched;
  impl::DecayingEventActivityMap<T> activity;
  SubscriberHandle handlerToken;
  std::vector<genlit> clause;
    std::vector<Literal<T>> nclause;
  //  int polarity{-1};
};
}

#endif //SCHEDCL_WEIGHTEDDEGREE_HPP
