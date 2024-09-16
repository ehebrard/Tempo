//
// Created by tim on 15.11.22.
//

#ifndef TEMPO_VSIDS_HPP
#define TEMPO_VSIDS_HPP

#include "Global.hpp"
#include "RankingHeuristic.hpp"
#include "heuristics/impl/DecayingEventActivityMap.hpp"
#include "util/SubscribableEvent.hpp"
#include "GapOverActivity.hpp"

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
    class VSIDS : public GapOverActivity {
    public:


        /**
         * CTor. Infers parameters from given scheduler. Subscribes to events of interest
         * @param solver target solver
         */
        template<concepts::scalar T>
        explicit VSIDS(Solver<T> &solver) :
                GapOverActivity(solver, solver.ClauseAdded.subscribe_handled(
                        [this, &solver](const auto &arg) { this->activity.update(arg, solver); })) {}

        // Copy and move are disabled. Otherwise, a call to the subscribed event handler will cause undefined behavior
        VSIDS(const VSIDS &) = delete;

        VSIDS(VSIDS &&) = delete;

        VSIDS &operator=(const VSIDS &) = delete;

        VSIDS &operator=(VSIDS &&) = delete;

        ~VSIDS() = default;

    };
}

#endif //TEMPO_VSIDS_HPP
