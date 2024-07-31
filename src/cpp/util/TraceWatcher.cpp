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

    void Tracer::handleConflict() {
        if (watcher.isOnTrack()) {
            DeviationOccurred.trigger(DeviationType::Fail, TraceWatcher::Conflicts{}, Decisions{});
        }

        watcher.setOnTrack(false);
    }

    auto Tracer::getWatcher() const noexcept -> const TraceWatcher & {
        return watcher;
    }
}
