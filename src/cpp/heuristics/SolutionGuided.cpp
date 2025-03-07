/**
* @author Tim Luchterhand
* @date 07.03.25
* @file SolutionGuided.cpp
* @brief
*/

#include "heuristics/SolutionGuided.hpp"
#include "util/Options.hpp"

namespace tempo::heuristics {
    TSGFactory::TSGFactory() {
        ValueHeuristicFactory::get().registerFactory(Options::PolarityHeuristic::TSG, this);
    }

    RSGFactory::RSGFactory() {
        ValueHeuristicFactory::get().registerFactory(Options::PolarityHeuristic::RSG, this);
    }

    static inline TSGFactory _tsgFactory;
    static inline RSGFactory _rsgFactory;
}
