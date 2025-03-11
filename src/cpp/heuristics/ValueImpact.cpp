/**
* @author Tim Luchterhand
* @date 11.03.25
* @file ValueImpact.cpp
* @brief
*/

#include "util/Options.hpp"
#include "heuristics/ValueImpact.hpp"

namespace tempo::heuristics {
    ValueImpactFactory::ValueImpactFactory() {
        ValueHeuristicFactory::get().registerFactory(Options::PolarityHeuristic::Impact, this);
    }

    static ValueImpactFactory _valueImpactFactory;
}
