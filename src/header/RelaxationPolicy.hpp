/************************************************
 * Tempo RelaxationPolicy.hpp
 *
 * Copyright 2024 Emmanuel Hebrard
 *
 * Tempo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 * Tempo is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Tempo.  If not, see <http://www.gnu.org/licenses/>.
 *
 ***********************************************/

#ifndef _TEMPO_RELAXATIONPOLICY_HPP
#define _TEMPO_RELAXATIONPOLICY_HPP

#include <vector>

#include "Model.hpp"

namespace tempo {

//! Relaxation Policies for LNS
template<typename T>
class RelaxationPolicy {
public:
    virtual void select(std::vector<Literal<T>>& fixed) = 0;
    virtual void notifySuccess() {}
    virtual void notifyFailure() {}
};



template<typename T>
class RelaxRandomDisjunctiveResource : public RelaxationPolicy<T> {
public:
    RelaxRandomDisjunctiveResource(Solver<T>& solver, std::vector<NoOverlapExpression<T>>& resources) : solver(solver), resources(resources) {}
     void select(std::vector<Literal<T>>& fixed) override;
    
private:
    Solver<T>& solver;
    std::vector<NoOverlapExpression<T>>& resources;
};


template<typename T>
void RelaxRandomDisjunctiveResource<T>::select(std::vector<Literal<T>>& fixed) {
    fixed.clear();
    int r{static_cast<int>(random() % resources.size())};
    for(auto &resource : resources) if(--r != 0) {
        for(auto bi{resource.begDisjunct()}; bi!=resource.endDisjunct(); ++bi) {
            fixed.push_back(*bi == solver.boolean.value(*bi));
        }
    }
}


template<typename T>
class FixRandomDisjunctiveResource : public RelaxationPolicy<T> {
public:
    FixRandomDisjunctiveResource(Solver<T>& solver, std::vector<NoOverlapExpression<T>>& resources) : solver(solver), resources(resources) {}
     void select(std::vector<Literal<T>>& fixed) override;
    
private:
    Solver<T>& solver;
    std::vector<NoOverlapExpression<T>>& resources;
};


template<typename T>
void FixRandomDisjunctiveResource<T>::select(std::vector<Literal<T>>& fixed) {
    fixed.clear();
    int r{static_cast<int>(random() % resources.size())};
        for(auto bi{resources[r].begDisjunct()}; bi!=resources[r].endDisjunct(); ++bi) {
            fixed.push_back(*bi == solver.boolean.value(*bi));
        }
}


} // namespace tempo

#endif
