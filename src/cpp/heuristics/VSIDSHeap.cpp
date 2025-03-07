/**
* @author Tim Luchterhand
* @date 07.03.25
* @file VSIDSHeap.cpp
* @brief
*/

#include "heuristics/VSIDSHeap.hpp"

namespace tempo::heuristics {
    VSIDSHeapFactory::VSIDSHeapFactory() {
        VariableHeuristicFactory::get().registerFactory(Options::ChoicePointHeuristics::VSIDSHeap, this);
    }

    static inline VSIDSHeapFactory _vsidsHeapFactory;
}
