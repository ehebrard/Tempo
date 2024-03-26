//
// Created by tluchterha on 22/11/22.
//

#ifndef SCHEDCL_EVENTACTIVITYMAP_HPP
#define SCHEDCL_EVENTACTIVITYMAP_HPP

#include <concepts>

#include "Scheduler.hpp"

namespace tempo {
    template<typename T>
    class Scheduler;
}

namespace tempo::heuristics::impl {
    /**
     * @brief Class that can be used to record the activity on distance constraints
     */
template<typename T>
    class EventActivityMap {
    public:
        /**
         * CTor. Initializes activity of all literals with 1.
         * @tparam T type of scheduler
         * @param scheduler scheduler for which to construct the ActivityMap
         */
        
        explicit EventActivityMap(Scheduler<T> &scheduler) : sched(scheduler) {
//            numNodes = scheduler.numEvent();
            activity.resize(sched.numEvent(), 1);
            sched.setActivityMap(this);
        }

        
        /**
         * Gets the activity for a given choice point
         * @tparam T
         * @param bound constraint
         * @return
         */
        constexpr double get(const BoundConstraint<T>& c) const noexcept {
            return activity[EVENT(c.l)];
        }
        
        /**
         * Gets the activity for a given choice point
         * @tparam T
         * @param edge constraint
         * @return
         */
        constexpr double get(const DistanceConstraint<T>& c) const noexcept {
            return activity[c.from] + activity[c.to];
        }
        
        /**
         * Gets the activity for a given choice point
         * @tparam T
         * @param variable
         * @return
         */
        constexpr double get(const var x) const noexcept {
            return get(sched.getEdge(POS(x))) + get(sched.getEdge(NEG(x)));
//            DistanceConstraint<T> left{sched.getEdge(POS(x))};
//            DistanceConstraint<T> right{sched.getEdge(NEG(x))};
//            return activity[left.from] + activity[left.to] + activity[right.from] + activity[right.to];
        }

//        /**
//         * Gets the activity for a given choice point
//         * @tparam T
//         * @param choicePoint
//         * @return
//         */
//        constexpr double &get(const var x) noexcept {
//            DistanceConstraint<T> left{sched.getEdge(POS(x))};
//            DistanceConstraint<T> right{sched.getEdge(NEG(x))};
//            return activity[left.from] + activity[left.to] + activity[right.from] + activity[right.to];
//        }

        /**
         * Applies the given functor to all entries in the activity map
         * @param functor
         */
        void for_each(const std::invocable<double &> auto &functor) {
            for (auto &val : activity) {
                functor(val);
            }
        }
        
        /**
         * Displays weights
         */
//        template<typename T>
        std::ostream& display(std::ostream& os) {
            for(auto i{0}; i<sched.numEvent(); ++i) {
                if(activity[i] > 1)
                    os << " " << prettyEvent(i) << ":" << activity[i];
            }
            return os;
        }

    protected:

        Scheduler<T>& sched;
        std::vector<double> activity{};
    };
}

#endif //SCHEDCL_EVENTACTIVITYMAP_HPP
