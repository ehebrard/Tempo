/**
 * @author Tim Luchterhand
 * @date 10.11.23
 */

#ifndef TEMPO_TASK_TIMINGS_HPP
#define TEMPO_TASK_TIMINGS_HPP
#include "Global.hpp"
#include "traits.hpp"

namespace tempo {
    /**
     * Reads out the minimum task duration from an event network
     * @tparam EvtFun type of event distance function
     * @param t task id
     * @param eventDistanceFunction event distance function from which to read out the timing values
     * @return minimum task duration
     */
    template<concepts::arbitrary_event_dist_fun EvtFun>
    constexpr auto minDuration(task t, const EvtFun &eventDistanceFunction) {
        return -eventDistanceFunction(END(t), START(t));
    }

    /**
     * Reads out the maximum task duration from an event network
     * @tparam EvtFun type of event distance function
     * @param t task id
     * @param eventDistanceFunction event distance function from which to read out the timing values
     * @return maximum task duration
     */
    template<concepts::arbitrary_event_dist_fun EvtFun>
    constexpr auto maxDuration(task t, const EvtFun &eventDistanceFunction) {
        return eventDistanceFunction(START(t), END(t));
    }

    /**
     * Reads out the earliest starting time of a task from an event network
     * @tparam EvtFun type of event distance function
     * @param t task id
     * @param eventDistanceFunction event distance function from which to read out the timing values
     * @return earliest task starting time
     */
    template<concepts::arbitrary_event_dist_fun EvtFun>
    constexpr auto earliestStartTime(task t, const EvtFun &eventDistanceFunction) {
        return eventDistanceFunction(START(t), ORIGIN);
    }

    /**
     * Reads out the latest deadline of a task from an event network
     * @tparam EvtFun type of event distance function
     * @param t task id
     * @param eventDistanceFunction event distance function from which to read out the timing values
     * @return latest task completion time
     */
    template<concepts::arbitrary_event_dist_fun EvtFun>
    constexpr auto latestCompletion(task t, const EvtFun &eventDistanceFunction) {
        return eventDistanceFunction(HORIZON, END(t));
    }

    /**
     * Reads out the distance between two tasks from an event network, i.e. Let a and b be two tasks, then the distance
     * between a and b is the distance between a's start event and b's end event
     * @tparam EvtFun type of event distance function
     * @param from id of the first task
     * @param to id of the second task
     * @param eventDistanceFunction event distance function from which to read out the timing values
     * @return distance between START(a) and END(b)
     */
    template<concepts::arbitrary_event_dist_fun EvtFun>
    constexpr auto taskDistance(task from, task to, const EvtFun &eventDistanceFunction) {
        return eventDistanceFunction(START(from), END(to));
    }

    /**
     * Reads out the upper bound from an event network, i.e. the distance between the origin and the horizon events
     * @tparam EvtFun type of event distance function
     * @param eventDistanceFunction event distance function from which to read out the timing values
     * @return distance between ORIGIN and HORIZON
     */
    template<concepts::arbitrary_event_dist_fun EvtFun>
    constexpr auto upperBound(const EvtFun &eventDistanceFunction) {
        return eventDistanceFunction(ORIGIN, HORIZON);
    }

    /**
     * Reads out the lower bound from an event network, i.e. the negative distance between the horizon and origin events
     * @tparam EvtFun type of event distance function
     * @param eventDistanceFunction event distance function from which to read out the timing values
     * @return negative distance between HORIZON and ORIGIN
     */
    template<concepts::arbitrary_event_dist_fun EvtFun>
    constexpr auto lowerBound(const EvtFun &eventDistanceFunction) {
        return -eventDistanceFunction(HORIZON, ORIGIN);
    }
}

#endif //TEMPO_TASK_TIMINGS_HPP
