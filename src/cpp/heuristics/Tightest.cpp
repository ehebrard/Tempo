/**
* @author Tim Luchterhand
* @date 07.03.25
* @file Tightest.cpp
* @brief
*/

#include "heuristics/Tightest.hpp"

namespace tempo::heuristics {
    TightestFactory::TightestFactory() {
        VariableHeuristicFactory::get().registerFactory(Options::ChoicePointHeuristics::Tightest, this);
    }

    static inline TightestFactory _tightestFactory{};
}
