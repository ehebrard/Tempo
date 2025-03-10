/**
* @author Tim Luchterhand
* @date 07.03.25
* @file RandomVariableSelection.cpp
* @brief
*/

#include "heuristics/RandomVariableSelection.hpp"

namespace tempo::heuristics {
    RandomVariableSelectionFactory::RandomVariableSelectionFactory() {
        VariableHeuristicFactory::get().registerFactory(Options::ChoicePointHeuristics::Random, this);
    }

    static inline RandomVariableSelectionFactory _randomVariableSelectionFactory{};
}
