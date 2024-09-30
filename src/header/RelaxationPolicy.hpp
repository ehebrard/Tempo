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
#include <ranges>

#include "Model.hpp"
#include "heuristics/RelaxationInterface.hpp"

namespace tempo {

//! Relaxation Policies for LNS
template<typename T>
class BaseRelaxationPolicy {
public:
    virtual void relax(heuristics::AssumptionInterface<T> &solver) = 0;
    virtual void notifySuccess(unsigned) {}
    virtual void notifyFailure() {}

    BaseRelaxationPolicy() = default;
    virtual ~BaseRelaxationPolicy() = default;
    BaseRelaxationPolicy(BaseRelaxationPolicy &&) = default;
    BaseRelaxationPolicy &operator=(BaseRelaxationPolicy &&) = default;
};



template<typename T>
class RelaxRandomDisjunctiveResource : public BaseRelaxationPolicy<T> {
public:
    RelaxRandomDisjunctiveResource(Solver<T>& solver, std::vector<NoOverlapExpression<T>>& resources) : solver(solver), resources(resources) {}
    void relax(heuristics::AssumptionInterface<T> &s) override;

private:
    Solver<T>& solver;
    std::vector<NoOverlapExpression<T>>& resources;
};


template<typename T>
void RelaxRandomDisjunctiveResource<T>::relax(heuristics::AssumptionInterface<T> &s) {
    int r{static_cast<int>(random() % resources.size())};

    std::cout << "relax resource " << r << "/" << resources.size() << std::endl;
    std::vector<Literal<T>> fixed;
    for (auto &resource: resources) {
        if (--r != 0) {
            for (auto bi{resource.begDisjunct()}; bi != resource.endDisjunct(); ++bi) {
                fixed.push_back(*bi == solver.boolean.value(*bi));
            }
        }
    }

    s.makeAssumptions(fixed);
}


template<typename T>
class FixRandomDisjunctiveResource : public BaseRelaxationPolicy<T> {
public:
    FixRandomDisjunctiveResource(Solver<T>& solver, std::vector<NoOverlapExpression<T>>& resources) : solver(solver), resources(resources) {}
    void relax(heuristics::AssumptionInterface<T> &s) override;
    
private:
    Solver<T>& solver;
    std::vector<NoOverlapExpression<T>>& resources;
};


template<typename T>
void FixRandomDisjunctiveResource<T>::relax(heuristics::AssumptionInterface<T> &s) {
    int r{static_cast<int>(random() % resources.size())};
    s.makeAssumptions(std::ranges::subrange(resources[r].begDisjunct(), resources[r].endDisjunct()));
}

template <typename T> class RandomSubset : public BaseRelaxationPolicy<T> {
public:
  RandomSubset(Solver<T> &solver, std::vector<BooleanVar<T>> &vars,
            const double ratio)
      : solver(solver), vars(vars), ratio(1.0 - ratio) {}
  void relax(heuristics::AssumptionInterface<T> &s) override;

private:
  Solver<T> &solver;
  std::vector<BooleanVar<T>> vars;
  double ratio;
};

template <typename T>
void RandomSubset<T>::relax(heuristics::AssumptionInterface<T> &s) {
  std::vector<Literal<T>> fixed;
  fixed.reserve(vars.size());
  for (auto x : vars) {
    fixed.push_back(x == solver.boolean.value(x));
  }

  size_t n{static_cast<size_t>(ratio * static_cast<double>(vars.size()))};
  for (size_t i{0}; i < n; ++i) {
    size_t r{i + static_cast<size_t>(random() % (n - i))};
    std::swap(fixed[i], fixed[r]);
  }
  fixed.resize(n);

  s.makeAssumptions(fixed);
}

} // namespace tempo

#endif
