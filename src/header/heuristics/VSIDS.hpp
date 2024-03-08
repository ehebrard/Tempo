//
// Created by tim on 15.11.22.
//

#ifndef SCHEDCL_VSIDS_HPP
#define SCHEDCL_VSIDS_HPP
#include <vector>
#include <optional>
#include "BaseHeuristic.hpp"
#include "util/traits.hpp"
#include "heuristics/impl/DecayingActivityMap.hpp"
#include "util/SubscribableEvent.hpp"

namespace schedcl::heuristics {
    /**
     * @brief VSIDS (Variable State Independent Decaying Sum) heuristic
     * @details @copybrief
     * Selects variables which occur in the largest number of clauses (proportionally to highest activity). Activity for
     * a literal l is calculated the following way:
     * \f$ A_l = \sum_{i=1}^n \gamma^{n - i} w_i \f$ where \f$ \gamma \f$ is a constant decay factor and \f$ n \f$ is
     * the the number of calls to the function VSIDS::updateActivity
     * The final score of a literal is inversely proportional to its activity (@see getCost)
     */
    template<typename Distance>
    class VSIDS : public BaseHeuristic<VSIDS<Distance>> {
    public:

        /**
         * CTor. Infers parameters from given scheduler. Subscribes to events of interest
         * @tparam T
         * @param scheduler
         * @param options command line options
         */
        template<typename T, typename DistanceFun>
        VSIDS(const Scheduler<T> &scheduler, const Options &options, DistanceFun &&distanceFun) :
                activity(scheduler, options.vsids_decay), distance(std::forward<DistanceFun>(distanceFun)),
                handlerToken(scheduler.ClauseAdded.subscribe_handled(
                        [this](const auto &arg) { this->activity.update(arg); })) {}

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
        template<typename T>
        double getCost(const DistanceConstraint<T> & choicePoint, std::size_t) const {
            return std::max(distance(choicePoint.from, choicePoint.to), distance(choicePoint.to, choicePoint.from)) /
                   activity.get(choicePoint);
        }

        template<typename T>
        constexpr void preEvaluation(const Scheduler<T> &) const noexcept {}

    private:
        impl::DecayingActivityMap activity;
        Distance distance;
        SubscriberHandle handlerToken;
    };
}

#endif //SCHEDCL_VSIDS_HPP
