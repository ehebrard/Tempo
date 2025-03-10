/**
* @author Tim Luchterhand
* @date 07.03.25
* @file WeightedDegree.cpp
* @brief
*/

#include "heuristics/WeightedDegree.hpp"

namespace tempo::heuristics {
    WeightedDegreeFactory::WeightedDegreeFactory() {
        VariableHeuristicFactory::get().registerFactory(Options::ChoicePointHeuristics::WeightedDegree, this);
    }

    static inline WeightedDegreeFactory _weightedDegreeFactory{};
}
