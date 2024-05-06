/**
* @author Tim Luchterhand
* @date 06.05.24
* @brief
*/

#ifndef TEMPO_BASEVALUEHEURISTIC_HPP
#define TEMPO_BASEVALUEHEURISTIC_HPP

#include <concepts>
#include <stdexcept>

#include "Global.hpp"
#include "util/crtp.hpp"
#include "util/traits.hpp"

namespace tempo {
    template<typename T>
    class Scheduler;
}

namespace tempo::heuristics {
    template<typename H, typename T>
    concept value_heuristic = requires(H heuristic, var x, const Scheduler<T> scheduler) {
        {heuristic.choose(x, scheduler)} -> std::same_as<lit>;
    };

    template<typename Impl>
    class BaseValueHeuristic : public crtp<Impl, BaseValueHeuristic> {
        static constexpr auto EpsScale = 10000ul;
        unsigned long epsilon;
    public:
        explicit constexpr BaseValueHeuristic(double epsilon): epsilon(static_cast<unsigned long>(epsilon * EpsScale)) {
            if (epsilon > 1 or epsilon < 0) {
                throw std::runtime_error("invalid epsilon value");
            }
        }

        template<concepts::scalar T> requires(value_heuristic<Impl, T>)
        constexpr lit choosePolarity(var cp, const Scheduler<T> &scheduler) {
            if (tempo::random() % EpsScale < epsilon) {
                return tempo::random() % 2 == 0 ? POS(cp) : NEG(cp);
            }

            return this->getImpl().choose(cp, scheduler);
        }
    };
}

#endif //TEMPO_BASEVALUEHEURISTIC_HPP
