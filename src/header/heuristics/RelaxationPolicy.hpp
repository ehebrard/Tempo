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
#include <filesystem>

#include "Model.hpp"
#include "RelaxationInterface.hpp"
#include "util/serialization.hpp"

namespace tempo::heuristics {

//! Relaxation Policies for LNS
template<typename T>
class BaseRelaxationPolicy {
public:
    virtual void relax(heuristics::AssumptionProxy<T> &solver) = 0;
    virtual void notifySuccess(unsigned) {}
    virtual void notifyFailure(unsigned) {}

    BaseRelaxationPolicy() = default;
    virtual ~BaseRelaxationPolicy() = default;
    BaseRelaxationPolicy(BaseRelaxationPolicy &&) = default;
    BaseRelaxationPolicy &operator=(BaseRelaxationPolicy &&) = default;
};



template<typename T>
class RelaxRandomDisjunctiveResource : public BaseRelaxationPolicy<T> {
public:
    RelaxRandomDisjunctiveResource(Solver<T>& solver, std::vector<NoOverlapExpression<T>>& resources) : solver(solver), resources(resources) {}
    void relax(heuristics::AssumptionProxy<T> &s) override;

private:
    Solver<T>& solver;
    std::vector<NoOverlapExpression<T>>& resources;
};


template<typename T>
void RelaxRandomDisjunctiveResource<T>::relax(heuristics::AssumptionProxy<T> &s) {
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
    void relax(heuristics::AssumptionProxy<T> &s) override;
    
private:
    Solver<T>& solver;
    std::vector<NoOverlapExpression<T>>& resources;
};


template<typename T>
void FixRandomDisjunctiveResource<T>::relax(heuristics::AssumptionProxy<T> &s) {
    int r{static_cast<int>(random() % resources.size())};
    s.makeAssumptions(std::ranges::subrange(resources[r].begDisjunct(), resources[r].endDisjunct()));
}

template <typename T> class RandomSubset : public BaseRelaxationPolicy<T> {
public:
  RandomSubset(Solver<T> &solver, std::vector<BooleanVar<T>> &vars,
            double ratio, double decay)
      : solver(solver), vars(vars), ratio(1.0 - ratio), decay(decay) {}
  void relax(heuristics::AssumptionProxy<T> &s) override;
  void notifyFailure(unsigned) override;

private:
  Solver<T> &solver;
  std::vector<BooleanVar<T>> vars;
  double ratio, decay;

  // helper
  std::vector<Literal<T>> fixed;
};

template<typename T>
void RandomSubset<T>::notifyFailure(unsigned) {
    ratio *= decay;
}

template <typename T>
void RandomSubset<T>::relax(heuristics::AssumptionProxy<T> &s) {
  fixed.clear();
  for (auto x : vars) {
    fixed.push_back(x == solver.boolean.value(x));
  }

  size_t n{static_cast<size_t>(ratio * static_cast<double>(vars.size()))};
  for (size_t i{0}; i < n; ++i) {
    size_t r{i + static_cast<size_t>(random() % (fixed.size() - i))};
    std::swap(fixed[i], fixed[r]);
  }
  fixed.resize(n);

  s.makeAssumptions(fixed);
}

/**
 * @brief Policy that exactly replays the assumptions made by another policy
 * @tparam T timing type
 */
template<typename T>
class PolicyReplay {
    std::vector<typename heuristics::LoggingAssumptionProxy<T>::PolicyTrace> trace;
    std::size_t idx = 0;
public:
    /**
     * Ctor
     * @param traceFile replay file
     */
    explicit PolicyReplay(const std::filesystem::path &traceFile) : trace(
            serialization::deserializeFromFile<decltype(trace)>(traceFile)) {}

    /**
     * Relaxation interface. Replays the assumptions previously made by the recorded policy
     * @param s assumption proxy to
     */
    void relax(heuristics::AssumptionProxy<T> &s) {
        using enum heuristics::PolicyAction;
        if (idx >= trace.size()) {
            throw std::runtime_error("Assumption list exhausted. Cannot make any more assumptions.");
        }

        for (const auto &[action, literals] : trace[idx]) {
            switch (action) {
                case Reset:
                    s.reset();
                    break;
                case TrySet:
                    assert(literals.size() == 1);
                    s.tryMakeAssumption(literals.front());
                    break;
                case Set:
                    s.makeAssumptions(literals);
                    break;
                default:
                    throw std::runtime_error("enum out of bounds");
            }
        }

        ++idx;
    }

    void notifyFailure(unsigned) const noexcept {}
    void notifySuccess(unsigned) const noexcept {}
};

} // namespace tempo

#endif
