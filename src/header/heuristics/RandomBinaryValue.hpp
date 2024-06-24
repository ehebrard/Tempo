/**
 * @author Tim Luchterhand
 * @date 07.05.24
 * @brief Random value selection
 */

#ifndef TEMPO_RANDOMBINARYVALUE_HPP
#define TEMPO_RANDOMBINARYVALUE_HPP

#include <cassert>

#include "BaseBooleanHeuristic.hpp"
#include "Global.hpp"
#include "util/factory_pattern.hpp"

namespace tempo::heuristics {

/**
 * @brief Random value selection heuristic.
 * @details @copybrief
 * Chooses choice point polarity always randomly
 */
    class RandomBinaryValue {
    public:
        /**
         * heuristic interface
         * @tparam Sched class that provides additional information for the actual
         * implementation
         * @param cp choice point
         * @param scheduler scheduler instance
         * @return either POS(cp) or NEG(cp)
         */
        template<boolean_info_provider Solver>
        static auto valueDecision(const VariableSelection &selection, const Solver &solver) noexcept {
            assert(selection.second == VariableType::Boolean);
            auto lit = solver.boolean.getLiteral(true, selection.first);
            return tempo::random() % 2 == 0 ? lit : ~lit;
        }
    };

    MAKE_DEFAULT_FACTORY(RandomBinaryValue, const ValueHeuristicConfig &)

} // namespace tempo::heuristics

#endif // TEMPO_RANDOMBINARYVALUE_HPP
