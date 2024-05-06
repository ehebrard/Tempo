/**
* @author Tim Luchterhand
* @date 06.05.24
* @brief
*/

#include <stdexcept>

#include "heuristics/ValueHeuristicsManager.hpp"

namespace tempo::heuristics {

    auto valHeuristicTypeToString(Options::PolarityHeuristic type) -> std::string {
        switch (type) {
            case Options::PolarityHeuristic::Identity:
                return "Identity";
            case Options::PolarityHeuristic::LocalExploration:
                return "SolutionGuided";
            case Options::PolarityHeuristic::Tightest:
                return "TightestValue";
            default:
                throw std::runtime_error("unknown value heuristic type");
        }
    }
}