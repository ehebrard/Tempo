/**
* @author Tim Luchterhand
* @date 07.05.24
* @brief
*/

#ifndef TEMPO_RANDOMVALUE_HPP
#define TEMPO_RANDOMVALUE_HPP
#include <cassert>

#include "Global.hpp"
#include "BaseValueHeuristic.hpp"
#include "util/factory_pattern.hpp"

namespace tempo::heuristics {

    class RandomValue {
    public:
        template<concepts::scalar T>
        static lit choosePolarity(var cp, const Scheduler <T> &) noexcept {
            return tempo::random() % 2 == 0 ? POS(cp) : NEG(cp);
        }
    };

    MAKE_DEFAULT_FACTORY(RandomValue, const ValueHeuristicConfig &)

}

#endif //TEMPO_RANDOMVALUE_HPP
