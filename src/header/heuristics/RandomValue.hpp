/**
 * @author Tim Luchterhand
 * @date 07.05.24
 * @brief Random value selection
 */

#ifndef TEMPO_RANDOMVALUE_HPP
#define TEMPO_RANDOMVALUE_HPP
#include <cassert>

#include "BaseValueHeuristic.hpp"
#include "Global.hpp"
#include "util/factory_pattern.hpp"

namespace tempo::heuristics {

/**
 * @brief Random value selection heuristic.
 * @details @copybrief
 * Chooses choice point polarity always randomly
 */
class RandomValue {
public:
  /**
   * heuristic interface
   * @tparam Sched class that provides additional information for the actual
   * implementation
   * @param cp choice point
   * @param scheduler scheduler instance
   * @return either POS(cp) or NEG(cp)
   */
  template <typename Sched>
  static lit choosePolarity(var cp, const Sched &) noexcept {
    return tempo::random() % 2 == 0 ? POS(cp) : NEG(cp);
  }
};

MAKE_DEFAULT_FACTORY(RandomValue, const ValueHeuristicConfig &)

} // namespace tempo::heuristics

#endif // TEMPO_RANDOMVALUE_HPP
