/**
* @author Tim Luchterhand
* @date 14.03.25
* @file GNNDispatcher.cpp
* @brief
*/

#include "nn/GNNDispatcher.hpp"

namespace tempo::nn {
    bool FullGuidance::runInference() {
        return true;
    }

    bool Never::runInference() {
        return false;
    }

    bool SingleShotDispatcher::runInference() {
        if (inferenceAllowed) {
            inferenceAllowed = false;
            return true;
        }

        return false;
    }
}
