/**
* @author Tim Luchterhand
* @date 07.03.25
* @file VSIDS.cpp
* @brief
*/

#include "heuristics/VSIDS.hpp"

namespace tempo::heuristics {
    VSIDSFactory::VSIDSFactory() {
        VariableHeuristicFactory::get().registerFactory(Options::ChoicePointHeuristics::VSIDS, this);
    }

    static inline VSIDSFactory _vsidsFactory{};
}
