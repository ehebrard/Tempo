/**
* @author Tim Luchterhand
* @date 07.03.25
* @file RandomBinaryValue.cpp
* @brief
*/

#include "heuristics/RandomBinaryValue.hpp"
#include "util/Options.hpp"

namespace tempo::heuristics {
    RandomBinaryValueFactory::RandomBinaryValueFactory() {
        ValueHeuristicFactory::get().registerFactory(Options::PolarityHeuristic::Random, this);
    }

    static inline RandomBinaryValueFactory _randomBinaryValueFactory;
}
