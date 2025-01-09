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
#include <variant>
#include <filesystem>

#include "Global.hpp"
#include "Model.hpp"
#include "relaxation_interface.hpp"
#include "util/Options.hpp"
#include "util/serialization.hpp"
#include "util/factory_pattern.hpp"
#include "util/random.hpp"
#include "heuristics/LNS/PolicyDecay.hpp"
#include "util/Lookup.hpp"

namespace tempo::lns {


template<resource_expression R>
class RelaxRandomDisjunctiveResource {
public:
    explicit RelaxRandomDisjunctiveResource(std::vector<R> resources) noexcept: resources(std::move(resources)) {}

    template<resource_range RR>
    explicit RelaxRandomDisjunctiveResource(const RR& resources): resources(resources.begin(), resources.end()) {}

    template<assumption_interface AI>
    void relax(AI &s) const;
    void notifyFailure(unsigned ) noexcept {}
    void notifySuccess(unsigned ) noexcept {}

private:
    std::vector<R> resources;
};


template<resource_expression R>
template<assumption_interface AI>
void RelaxRandomDisjunctiveResource<R>::relax(AI &s) const {
    using namespace std::views;
    int r{static_cast<int>(random() % resources.size())};
    if (s.getSolver().getOptions().verbosity >= Options::YACKING) {
        std::cout << "relax resource " << r << "/" << resources.size() << std::endl;
    }

    auto skipView = resources | filter([r](const auto &) mutable { return r-- != 0; });
    auto vars = booleanVarsFromResources(skipView);
    s.makeAssumptions(vars | transform([&s](const auto &var) { return var == s.getSolver().boolean.value(var); }));
}


template<resource_expression R>
class FixRandomDisjunctiveResource {
public:
    explicit FixRandomDisjunctiveResource(std::vector<R> resources) : resources(std::move(resources)) {}

    template<resource_range RR>
    explicit FixRandomDisjunctiveResource(const RR& resources): resources(resources.begin(), resources.end()) {}

    void notifyFailure(unsigned ) noexcept {}
    void notifySuccess(unsigned ) noexcept {}

    template<assumption_interface AI>
    void relax(AI &s) const;
    
private:
    std::vector<R> resources;
};


template<resource_expression R>
template<assumption_interface AI>
void FixRandomDisjunctiveResource<R>::relax(AI &s) const {
    int r{static_cast<int>(random() % resources.size())};
    const auto variables = booleanVarsFromResources(resources[r]);
    auto assumptions = variables |
                       std::views::transform([&s](const auto &var) { return var == s.getSolver().boolean.value(var); });
    s.makeAssumptions(assumptions);
}

template <typename T> class RandomSubset {
public:
    RandomSubset(std::vector<BooleanVar<T>> vars, double relaxRatio, double decay) : vars(std::move(vars)),
        ratio(1.0 - relaxRatio), decay(decay) {}

  template<assumption_interface AI>
  void relax(AI &s) const;
  void notifyFailure(unsigned ) noexcept;
  void notifySuccess(unsigned ) noexcept {}

protected:
  std::vector<BooleanVar<T>> vars;
  double ratio, decay;

  // helper
};

template<typename T>
void RandomSubset<T>::notifyFailure(unsigned) noexcept {
    ratio *= decay;
}

template <typename T>
template<assumption_interface AI>
void RandomSubset<T>::relax(AI &s) const {
  std::vector<Literal<T>> fixed;
  fixed.reserve(vars.size());
  for (auto x : vars) {
    fixed.push_back(x == s.getSolver().boolean.value(x));
  }

  size_t n{static_cast<size_t>(ratio * static_cast<double>(vars.size()))};
  for (size_t i{0}; i < n; ++i) {
    size_t r{i + static_cast<size_t>(random() % (fixed.size() - i))};
    std::swap(fixed[i], fixed[r]);
  }

  fixed.resize(n);
  s.makeAssumptions(fixed);
  if (s.getSolver().getOptions().verbosity >= Options::YACKING) {
      std::cout << "-- fixing " << fixed.size() << " literals\n";
  }
}

namespace detail {
    template<concepts::scalar T>
    class TaskVarMap {
        Lookup<Interval<T>, std::vector<std::pair<int, BooleanVar<T>>>, IdProjection> map{}; // target task id, corresponding variable
        mutable Lookup<Interval<T>, bool, IdProjection> relaxed;
    public:
        template<concepts::typed_range<Interval<T>> Tasks, resource_range RR>
        TaskVarMap(const Tasks &tasks, const RR &resources) {
            using namespace std::views;
            map = decltype(map)(tasks);
            relaxed = decltype(relaxed)(tasks, false, {}, std::pair{map.keyOffset(), map.maxKey()});
            for (const auto &t : tasks) {
                auto &tVars = map[t];
                for (const auto &r : resources) {
                    auto res = std::ranges::find(r, t);
                    if (res == std::ranges::end(r)) {
                        continue;
                    }

                    const auto idx = std::ranges::distance(std::ranges::begin(r), res);
                    for (auto [targetTask, lit]: iterators::enumerate(r.getDisjunctiveLiterals().row_unsafe(idx))) {
                        if (lit == Contradiction<T>) {
                            continue;
                        }

                        tVars.emplace_back(targetTask, lit);
                    }
                }

                auto comp = [](const auto &pair) { return pair.second.id(); };
                std::ranges::sort(tVars, {}, comp);
                auto res = std::ranges::unique(tVars, {}, comp);
                tVars.erase(res.begin(), res.end());
                tVars.shrink_to_fit();
            }
        }

        template<concepts::same_template<Interval> Task>
        auto operator()(const Task &t) const noexcept -> const auto & {
            return map[t];
        }

        template<concepts::typed_range<Interval<T>> Tasks>
        auto getTaskLiterals(const Tasks &tasks, bool allEdges) const -> std::vector<BooleanVar<T>> {
            using namespace std::views;
            if (not allEdges) {
                relaxed.data().assign(relaxed.size(), false);
                for (const auto &t : tasks) {
                    relaxed[t] = true;
                }
            }

            auto varsView = tasks | transform([this](const auto &t) -> decltype(auto) { return (*this)(t); }) | join
                            | filter([this, allEdges](const auto &tpl) {
                                return allEdges or relaxed.data()[std::get<0>(tpl)];
                            }) | elements<1>;
            std::vector<BooleanVar<T>> vars;
            std::ranges::copy(varsView, std::back_inserter(vars));
            std::ranges::sort(vars);
            auto res = std::ranges::unique(vars);
            vars.erase(res.begin(), res.end());
            return vars;
        }
    };

}

/**
 * @brief Relaxation policy that relaxes a number of tasks
 * @tparam T timing type
 */
template<concepts::scalar T>
class RelaxTasks {
    std::vector<Interval<T>> tasks;
    detail::TaskVarMap<T> map;
    PolicyDecay decayHandler;
    bool allTaskEdges;
public:
    /**
     * Ctor
     * @tparam RR resource range type
     * @param tasks vector of all tasks to consider
     * @param resources resource expressions in the problem
     * @param decayConfig dynamic relaxation ratio decay config
     * @param allTaskEdges whether to fix all edges of a task or only these between other fixed tasks
     * @param verbosity logging verbosity
     */
    template<resource_range RR>
    RelaxTasks(std::vector<Interval<T>> tasks, const RR &resources, const PolicyDecayConfig &decayConfig,
               bool allTaskEdges, int verbosity = Options::NORMAL):
            tasks(std::move(tasks)), map(this->tasks, resources),
            decayHandler(decayConfig,map.getTaskLiterals(this->tasks, allTaskEdges).size(),verbosity),
            allTaskEdges(allTaskEdges) {
        if (verbosity >= Options::YACKING) {
            std::cout << std::boolalpha << "-- fix all task edges: " << allTaskEdges << std::noboolalpha << std::endl;
            std::cout << decayConfig << std::endl;
        }
    }

    void notifyFailure(unsigned numFails) noexcept {
        decayHandler.notifyFailure(numFails);
    }

    void notifySuccess(unsigned numFails) noexcept {
        decayHandler.notifySuccess(numFails);
    }

    template<assumption_interface AI>
    void relax(AI &proxy) {
        using namespace std::views;
        const auto numFix = static_cast<std::size_t>(decayHandler.getFixRatio() * tasks.size());
        if (numFix == 0) {
            return;
        }

        std::ranges::shuffle(tasks, RNG{});
        auto vars = map.getTaskLiterals(counted(tasks.begin(), numFix), allTaskEdges);
        if (proxy.getSolver().getOptions().verbosity >= Options::YACKING) {
            std::cout << "-- fixing " << numFix << " / " << tasks.size()
                      << " tasks (" << vars.size() << ") variables" << std::endl;
        }

        proxy.makeAssumptions(
                vars | transform([&b = proxy.getSolver().boolean](const auto &var) { return var == b.value(var); }));
    }
};

/**
 * @brief Relaxation policy that chronologically relaxes slices of the best known schedule.
 * @details @copybrief
 * On a fail, the policy tries to relax the next slice in the schedule.
 * @tparam T timing type
 * @note using this policy as is makes the LNS incomplete. Use it together with other policies (e.g. SporadicRootSearch)
 */
template<concepts::scalar T>
class RelaxChronologically {
    std::vector<Interval<T>> tasks;
    detail::TaskVarMap<T> map;
    unsigned numberOfSlices;
    unsigned sliceWidth;
    unsigned sliceIdx = 0;
    bool allTaskEdges;

public:

    /**
     * Ctor
     * @tparam RR resource range type
     * @param tasks vector of all tasks to consider
     * @param resources resource expressions in the problem
     * @param numberOfSlices number of slices to split the schedule
     * @param allTaskEdges whether to fix all edges of a task or only these between other fixed tasks
     */
    template<resource_range RR>
    RelaxChronologically(std::vector<Interval<T>> tasks, const RR &resources, unsigned numberOfSlices,
                         bool allTaskEdges): tasks(std::move(tasks)), map(this->tasks, resources),
                                             numberOfSlices(
                                                 std::max(numberOfSlices, static_cast<unsigned>(tasks.size()))),
                                             sliceWidth(ceil_division(this->tasks.size(),
                                                                      static_cast<std::size_t>(this->numberOfSlices))),
                                             allTaskEdges(allTaskEdges) {
        if (numberOfSlices == 0) {
            throw std::runtime_error("number of slices cannot be 0");
        }
    }

    void notifyFailure(unsigned ) {
        sliceIdx = (sliceIdx + 1) % numberOfSlices;
    }

    void notifySuccess(unsigned ) {
        sliceIdx = 0;
    }

    template<assumption_interface AI>
    void relax(AI &proxy) {
        if (proxy.getSolver().getOptions().verbosity >= Options::YACKING) {
            std::cout << "-- relaxing slice " << sliceIdx << " of last schedule (~" << sliceWidth << " / "
                      << tasks.size() << " tasks)" << std::endl;
        }

        std::ranges::sort(tasks, {}, [&s = proxy.getSolver()](const auto &t) { return t.getLatestStart(s); });
        for (std::size_t idx = 0; idx < numberOfSlices; ++idx) {
            if (idx == sliceIdx) {
                continue;
            }

            std::ranges::subrange tRange(tasks.cbegin() + idx * sliceWidth,
                                         tasks.cbegin() + std::min((idx + 1) * sliceWidth, tasks.size()));
            auto vars = map.getTaskLiterals(tRange, allTaskEdges);
            bool success = proxy.makeAssumptions(vars | std::views::transform(
                    [&b = proxy.getSolver().boolean](const auto &var) { return var == b.value(var); }));
            if (not success) {
                return;
            }
        }
    }
};

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
public:
    /**
     * Ctor
     * @tparam Pol base relaxation policy type
     * @param probabilityIncrement increment to apply to the relaxation probability on failed relaxation
     * @param policy base relaxation policy
     */
    template<relaxation_policy Pol>
    SporadicRootSearch(double probabilityIncrement, Pol &&policy) :
            basePolicy(std::forward<Pol>(policy)), probabilityIncrement(probabilityIncrement) {
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
        if (s.getSolver().getOptions().verbosity >= Options::YACKING) {
            std::cout << "-- sporadic probability " << std::setprecision(2) << rootSearchProbability * 100 << "%"
                      << std::endl;
        }

        if (random_event_occurred<Resolution>(rootSearchProbability)) {
            if (s.getSolver().getOptions().verbosity >= Options::YACKING) {
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
 * @return constructed SporadicRootSearch policy
 */
template<relaxation_policy P>
auto make_sporadic_root_search(double rootProbabilityIncrement, P &&policy) {
    return SporadicRootSearch<P>(rootProbabilityIncrement, std::forward<P>(policy));
}

/**
 * @brief Policy that exactly replays the assumptions made by another policy
 * @tparam T timing type
 */
template<typename T>
class PolicyReplay {
    std::vector<typename LoggingAssumptionProxy<T>::PolicyTrace> trace;
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
    void relax(AssumptionProxy<T> &s) {
        using enum PolicyAction;
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


template<relaxation_policy ...R>
struct VariantPolicy : public std::variant<R...> {
    using std::variant<R...>::variant;

    template<assumption_interface AI>
    DYNAMIC_DISPATCH(relax, AI& proxy, &, proxy, EMPTY)

    DYNAMIC_DISPATCH(notifyFailure, unsigned numFail, =, numFail, EMPTY)

    DYNAMIC_DISPATCH(notifySuccess, unsigned numFail, =, numFail, EMPTY)

};

} // namespace tempo

#endif
