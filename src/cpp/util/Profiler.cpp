/**
* @author Tim Luchterhand
* @date 09.08.2024
* @brief
*/

#include <utility>

#include "util/Profiler.hpp"

namespace tempo::util {
    TimingEvent::TimingEvent(const detail::TP &start, const detail::TP &end) noexcept : start(start), end(end) {}

    void StopWatch::start() {
        startTp = std::chrono::high_resolution_clock::now();
    }

    TimingEvent StopWatch::getTiming() const {
        return {startTp, std::chrono::high_resolution_clock::now()};
    }

    StopWatch::StopWatch() {
        start();
    }

    ScopeWatch::ScopeWatch(Profiler &profiler, std::string eventName)
            : StopWatch(), profiler(profiler), name(std::move(eventName)) {
        start();
    }

    ScopeWatch::~ScopeWatch() {
        profiler.addEvent(getTiming(), name);
    }

    void Profiler::addEvent(const TimingEvent &event, const std::string &name) {
        events[name].emplace_back(event);
    }

    void Profiler::addEvent(detail::TP start, detail::TP end, const std::string &name) {
        events[name].emplace_back(start, end);
    }
}
