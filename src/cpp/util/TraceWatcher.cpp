/**
* @author Tim Luchterhand
* @date 25.07.24
* @brief
*/

#include "util/TraceWatcher.hpp"

namespace tempo {

    TraceWatcher::TraceWatcher(std::size_t numVariables): varPolarity(numVariables, false), onTrack(false) {}

    auto TraceWatcher::getLastSolution() const noexcept -> const std::vector<bool> & {
        return varPolarity;
    }

    bool TraceWatcher::isOnTrack() const noexcept {
        return onTrack;
    }

    void TraceWatcher::setOnTrack(bool truthVal) noexcept {
        onTrack = truthVal;
    }
}
