//
// Created by tluchterha on 22/11/22.
//

#ifndef TEMPO_WEIGHTEDCRITICALPATH_HPP
#define TEMPO_WEIGHTEDCRITICALPATH_HPP
#include "RankingHeuristic.hpp"
#include "Options.hpp"
#include "heuristics/impl/DecayingEventActivityMap.hpp"
#include "util/SubscribableEvent.hpp"
#include <array>

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
class WeightedCriticalPath : public BaseHeuristic<WeightedCriticalPath<T>> {
public:
  /**
   * @copydoc VSIDS::VSIDS
   */
  explicit WeightedCriticalPath(const Scheduler<T> &scheduler, const Options &options)
      : sched(scheduler), activity(scheduler, options.vsids_decay),
        handlerToken(scheduler.ConflictEncountered.subscribe_handled(
            [this](auto &expl) {
                this->clause.clear();
                if(expl.getType() == CYCLEEXPL)
                    expl.explain(NoLit,this->clause);
                else
                    sched.getCriticalPath(this->clause);
              this->activity.update(this->clause);
            })) {
  }

  /**
   * @copydoc VSIDS::getCost
   */
  double getCost(const var x) const {
      auto prec_a{sched.getEdge(POS(x))};
      auto prec_b{sched.getEdge(NEG(x))};
      auto gap_a = sched.upper(prec_a.from) - sched.lower(prec_a.to);
      auto gap_b = sched.upper(prec_b.from) - sched.lower(prec_b.to);
      return static_cast<double>(std::max(gap_a, gap_b)) / activity.get(x);
  }

  WeightedCriticalPath(const WeightedCriticalPath &) = delete;
  WeightedCriticalPath(WeightedCriticalPath &&) = delete;
  WeightedCriticalPath &operator=(const WeightedCriticalPath &) = delete;
  WeightedCriticalPath &operator=(WeightedCriticalPath &&) = delete;
  ~WeightedCriticalPath() = default;

private:
   const Scheduler<T>& sched;
  impl::DecayingEventActivityMap<T> activity;
  SubscriberHandle handlerToken;
    std::vector<genlit> clause;
};
}

#endif //TEMPO_WEIGHTEDCRITICALPATH_HPP
