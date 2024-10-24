/**
* @author Tim Luchterhand
* @date 16.09.24
* @brief Variable selection heuristics that orders variables by their respective gap in the temporal network divided by
 * variable activity
*/

#ifndef TEMPO_GAPOVERACTIVITY_HPP
#define TEMPO_GAPOVERACTIVITY_HPP

#include <algorithm>

#include "Global.hpp"
#include "Constant.hpp"
#include "RankingHeuristic.hpp"
#include "heuristics/impl/DecayingEventActivityMap.hpp"
#include "util/SubscribableEvent.hpp"
#include "util/distance.hpp"

namespace tempo {
    template<typename T>
    class Solver;
}

namespace tempo::heuristics {

    /**
     * @brief Variable selection heuristics that orders variables by their respective gap in the temporal network
     * divided by variable activity.
     * @details @copybrief
     * This is the base class for heuristics like VSIDS and WeightedDegree
     */
    class GapOverActivity : public RankingHeuristic<GapOverActivity> {
    public:

        /**
         * CTor. Initializes activity map.
         * @param solver target solver
         * @param handle subscriber handle of event handler that updates activity map
         */
        template<concepts::scalar T>
        explicit GapOverActivity(Solver<T> &solver, SubscriberHandle handle) :
                activity(solver, solver.getOptions().vsids_decay), handlerToken(std::move(handle)) {}

        // Copy and move are disabled. Otherwise, a call to the subscribed event handler will cause undefined behavior
        GapOverActivity(const GapOverActivity&) = delete;
        GapOverActivity(GapOverActivity &&) = delete;
        GapOverActivity &operator=(const GapOverActivity&) = delete;
        GapOverActivity &operator=(GapOverActivity &&) = delete;
        ~GapOverActivity() = default;

        template<concepts::scalar T>
        [[nodiscard]] double getCost(const var_t x, const Solver<T> &solver) const {
            //@TODO: there shoud be a normalization thingy and Boolean variables without semantic should get the highest value
            double dom{1};
            if (solver.boolean.hasSemantic(x)) {
                auto gapA = boundEstimation(true, x, solver).value_or(Constant::Infinity<T>);
                auto gapB = boundEstimation(false, x, solver).value_or(Constant::Infinity<T>);
                if(gapA == Constant::Infinity<T>) {
                    dom = static_cast<double>(gapA/2 + gapB/2);
                } else if(gapB == Constant::Infinity<T>) {
                    dom = static_cast<double>(gapA/2 + gapB/2);
                } else {
                    dom = static_cast<double>(std::max(gapA, gapB));
                }
            }
            auto act{activity.get(x, solver)};
            return dom / act;
        }

        /**
         * Ranking heuristic interface. Out of two variables selects the one for which gap_v / activity_v is lower.
         * Herby gap_v is the maximum distance on the associated distance constraint arcs and activity_v is the
         * variables activity.
         * @tparam T timing type
         * @param x variable x
         * @param y variable y
         * @param solver solver that provides distance information
         * @return variable according to criteria described above
         */
        template<concepts::scalar T>
        [[nodiscard]] var_t chooseBest(var_t x, var_t y, const Solver<T> &solver) const {
            return getCost(x, solver) <= getCost(y, solver) ? x : y;
        }


        /**
         * Heuristic interface.
         * @param solver
         * @todo currently only selects boolean variables
         */
        template<concepts::scalar T>
        auto nextVariable(const Solver<T> &solver) const -> VariableSelection {
            return {this->bestVariable(solver.getBranch(), solver), VariableType::Boolean};
        }

    protected:
        impl::DecayingEventActivityMap activity;
        SubscriberHandle handlerToken;
    };
}



#endif //TEMPO_GAPOVERACTIVITY_HPP
