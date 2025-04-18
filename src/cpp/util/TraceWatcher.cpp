/**
* @author Tim Luchterhand
* @date 25.07.24
* @brief
*/

#include "util/TraceWatcher.hpp"

namespace tempo {

    auto TraceWatcher::getLastSolution() const noexcept -> const std::vector<bool> & {
        return varPolarity;
    }

    bool TraceWatcher::isOnTrack() const noexcept {
        return onTrack;
    }

    void TraceWatcher::setOnTrack(bool truthVal) noexcept {
        onTrack = truthVal;
    }

    var_t TraceWatcher::getOffset() const noexcept {
        return offset;
    }

    auto TraceWatcher::getVariablesOnTrack() const noexcept -> const serialization::Branch & {
        return varsOnTrack;
    }

    void Tracer::handleConflict() {
        if (watcher.isOnTrack()) {
            DeviationOccurred.trigger(DeviationType::Fail, TraceWatcher::Conflicts{}, watcher.getVariablesOnTrack());
        }

        watcher.setOnTrack(false);
    }

    auto Tracer::getWatcher() const noexcept -> const TraceWatcher & {
        return watcher;
    }
}
