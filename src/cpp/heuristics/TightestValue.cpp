/**
* @author Tim Luchterhand
* @date 07.03.25
* @file TightestValue.cpp
* @brief
*/

#include "heuristics/TightestValue.hpp"
#include "util/Options.hpp"

namespace tempo::heuristics {
    TightestValueFactory::TightestValueFactory() {
        ValueHeuristicFactory::get().registerFactory(Options::PolarityHeuristic::Tightest, this);
    }

    static inline TightestValueFactory _tightestValueFactory;
}
