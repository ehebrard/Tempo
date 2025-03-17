/**
* @author Tim Luchterhand
* @date 07.03.25
* @file LearningRateHeap.cpp
* @brief
*/

#include "heuristics/LearningRateHeap.hpp"

namespace tempo::heuristics {
LearningRateHeapFactory::LearningRateHeapFactory() {
        VariableHeuristicFactory::get().registerFactory(Options::ChoicePointHeuristics::LRB, this);
    }

    static inline LearningRateHeapFactory _lrbHeapFactory;
}
