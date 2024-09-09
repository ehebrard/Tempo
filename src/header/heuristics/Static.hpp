//
// Created by tluchterha on 21/11/22.
//

#ifndef TEMPO_STATIC_HPP
#define TEMPO_STATIC_HPP

#include "ReversibleObject.hpp"

namespace tempo::heuristics {

    class Static {
    public:
        template<concepts::scalar T>
        Static(Solver<T> &solver) : next(0, &(solver.getEnv())) {
            for(auto x : solver.getBranch()) {
                ordering.push_back(x);
            }
        }
        
        template<concepts::scalar T>
        Static(Solver<T> &solver, std::vector<var_t>& o) : ordering(o), next(0, &(solver.getEnv())) {}

        /**
         * @tparam T
         * @param solver
         * @todo currently only selects boolean variables
         */
        template<concepts::scalar T>
        [[nodiscard]] auto nextVariable(const Solver<T> &solver) -> VariableSelection {
            unsigned i = static_cast<unsigned>(next);
            while(i < ordering.size()) {
                if(solver.getBranch().has(ordering[i])) {
                    next = i;
                    return {ordering[i], VariableType::Boolean};
                }
                ++i;
            }
            next = i;
            return {solver.getBranch().front(), VariableType::Boolean};
        }

    private:
        std::vector<var_t> ordering;
        Reversible<unsigned> next;
    };
}


#endif //TEMPO_TIGHTEST_HPP
