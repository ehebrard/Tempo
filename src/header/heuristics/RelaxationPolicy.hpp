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

#include "Global.hpp"
#include "Model.hpp"
#include "RelaxationInterface.hpp"
#include "util/serialization.hpp"

namespace tempo::heuristics {


template<typename T>
class RelaxRandomDisjunctiveResource {
public:
    RelaxRandomDisjunctiveResource(Solver<T>& solver, std::vector<NoOverlapExpression<T>>& resources) : solver(solver), resources(resources) {}
    template<assumption_interface AI>
    void relax(AI &s);
    void notifyFailure(unsigned ) noexcept {}
    void notifySuccess(unsigned ) noexcept {}

private:
    Solver<T>& solver;
    std::vector<NoOverlapExpression<T>>& resources;
};


template<typename T>
template<assumption_interface AI>
void RelaxRandomDisjunctiveResource<T>::relax(AI &s) {
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
class FixRandomDisjunctiveResource {
public:
    FixRandomDisjunctiveResource(const Solver<T> &solver, std::vector<NoOverlapExpression<T>> resources) :
            solver(solver), resources(std::move(resources)) {}
    void notifyFailure(unsigned ) noexcept {}
    void notifySuccess(unsigned ) noexcept {}

    template<assumption_interface AI>
    void relax(AI &s);
    
private:
    const Solver<T>& solver;
    std::vector<NoOverlapExpression<T>>& resources;
};


template<typename T>
template<assumption_interface AI>
void FixRandomDisjunctiveResource<T>::relax(AI &s) {
    int r{static_cast<int>(random() % resources.size())};
    auto assumptions = std::ranges::subrange(resources[r].begDisjunct(), resources[r].endDisjunct()) |
                       std::views::transform([this](const auto &var) { return var == solver.boolean.value(var); });
    s.makeAssumptions(assumptions);
}

template <typename T> class RandomSubset {
public:
  RandomSubset(Solver<T> &solver, std::vector<BooleanVar<T>> vars,
            double ratio, double decay)
      : solver(solver), vars(std::move(vars)), ratio(1.0 - ratio), decay(decay) {}

  template<assumption_interface AI>
  void relax(AI &s);
  void notifyFailure(unsigned ) noexcept;
  void notifySuccess(unsigned ) noexcept {}

protected:
  Solver<T> &solver;
  std::vector<BooleanVar<T>> vars;
  double ratio, decay;

  // helper
  std::vector<Literal<T>> fixed;
};

template<typename T>
void RandomSubset<T>::notifyFailure(unsigned) noexcept {
    ratio *= decay;
}

template <typename T>
template<assumption_interface AI>
void RandomSubset<T>::relax(AI &s) {
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
 * @brief Relaxation policy wrapper that randomly triggers a search without relaxation
 * @details @copybrief
 * The root search is triggered with a probability that increases on failed relaxation runs
 * @tparam P base relaxation policy type
 */
template<relaxation_policy P>
class SporadicRootSearch {
    static constexpr auto Resolution = 100000;
    P basePolicy;
    double rootSearchProbability = 0;
    double probabilityIncrement;
    int verbosity;
public:
    /**
     * Ctor
     * @tparam Pol base relaxation policy type
     * @param probabilityIncrement increment to apply to the relaxation probability on failed relaxation
     * @param policy base relaxation policy
     * @param verbosity logging verbosity
     */
    template<relaxation_policy Pol>
    SporadicRootSearch(double probabilityIncrement, Pol &&policy, int verbosity = Options::NORMAL) :
            basePolicy(std::forward<Pol>(policy)), probabilityIncrement(probabilityIncrement), verbosity(verbosity) {
        if (probabilityIncrement > 1 or probabilityIncrement < 0) {
            throw std::runtime_error("invalid probability increment");
        }
    }

    /**
     * relaxation policy interface
     * @tparam AI assumption interface type
     * @param s solver proxy
     */
    template<assumption_interface AI>
    void relax(AI &s) {
        if (verbosity >= Options::YACKING) {
            std::cout << "-- sporadic probability " << std::setprecision(2) << rootSearchProbability * 100 << "%"
                      << std::endl;
        }
        const bool root = random() % Resolution < static_cast<unsigned long>(rootSearchProbability * Resolution);
        if (root) {
            if (verbosity >= Options::YACKING) {
                std::cout << "-- root search" << std::endl;
            }
            return;
        }

        basePolicy.relax(s);
    }

    /**
     * Notifies the underlying policy of a failure and increments the root search probability
     * @param numFailures current number of failures of the solver
     */
    void notifyFailure(unsigned numFailures) {
        rootSearchProbability += probabilityIncrement;
        basePolicy.notifyFailure(numFailures);
    }

    /**
     * Notifies the underlying policy of a successful relaxation and resets the root search probability to 0
     * @param numFailures current number of failures of the solver
     */
    void notifySuccess(unsigned numFailures) {
        rootSearchProbability = 0;
        basePolicy.notifySuccess(numFailures);
    }
};

/**
 * Helper method for constructing sporadic root search policies facilitating template argument deduction
 * @tparam P relaxation policy type
 * @param rootProbabilityIncrement increment to apply to the relaxation probability on failed relaxation
 * @param policy base relaxation policy
 * @param verbosity logging verbosity
 * @return constructed SporadicRootSearch policy
 */
template<relaxation_policy P>
auto make_sporadic_root_search(double rootProbabilityIncrement, P &&policy, int verbosity = Options::NORMAL) {
    return SporadicRootSearch<P>(rootProbabilityIncrement, std::forward<P>(policy), verbosity);
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
