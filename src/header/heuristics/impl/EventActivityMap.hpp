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
        
        explicit EventActivityMap(const Scheduler<T> &scheduler) : sched(scheduler) {
//            numNodes = scheduler.numEvent();
            activity.resize(sched.numEvent(), 1);
        }

        /**
         * Gets the activity for a given choice point
         * @tparam T
         * @param choicePoint
         * @return
         */
//        template<typename T>
        constexpr double get(const var x) const noexcept {
            DistanceConstraint<T> left{sched.getEdge(POS(x))};
            DistanceConstraint<T> right{sched.getEdge(NEG(x))};
            
//            std::cout << "act[" << prettyEvent(left.from) << "]=" << activity[left.from] << ", act[" << prettyEvent(left.to) << "]=" << activity[left.to] << ", act[" << prettyEvent(right.from) << "]=" << activity[right.from] << ", act[" << prettyEvent(right.to) << "]=" << activity[right.to] << std::endl;
            
            return activity[left.from] + activity[left.to] + activity[right.from] + activity[right.to];
        }

        /**
         * Gets the activity for a given choice point
         * @tparam T
         * @param choicePoint
         * @return
         */
//        template<typename T>
        constexpr double &get(const var x) noexcept {
            DistanceConstraint<T> left{sched.getEdge(POS(x))};
            DistanceConstraint<T> right{sched.getEdge(NEG(x))};
            
//            std::cout << "act[" << prettyEvent(left.from) << "]=" << activity[left.from] << ", act[" << prettyEvent(left.to) << "]=" << activity[left.to] << ", act[" << prettyEvent(right.from) << "]=" << activity[right.from] << ", act[" << prettyEvent(right.to) << "]=" << activity[right.to] << std::endl;
            
            return activity[left.from] + activity[left.to] + activity[right.from] + activity[right.to];
//            return activity[choicePoint.from] + activity[choicePoint.to];
        }

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
//            for(auto i{0}; i<numNodes; ++i) {
//                os << std::setw(5) << i;
//            }
//            os << std::endl;
            for(auto i{0}; i<sched.numEvent(); ++i) {
                if(activity[i] > 1)
                    os << " " << prettyEvent(i) << ":" << activity[i];
            }
            return os;
        }

    protected:

        const Scheduler<T>& sched;
        std::vector<double> activity{};
//        std::size_t numNodes{};
    };
}

#endif //SCHEDCL_EVENTACTIVITYMAP_HPP
