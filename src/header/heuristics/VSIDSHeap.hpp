/**
* @author Emmanuel Hebrard
* @date 12.09.24
* @brief
*/

#ifndef TEMPO_VSIDSHEAP_HPP
#define TEMPO_VSIDSHEAP_HPP


//#include <algorithm>

#include "heuristics/impl/ActivityMap.hpp"
#include "util/Heap.hpp"
#include "util/SparseSet.hpp"
#include "util/traits.hpp"

#include "heuristic_interface.hpp"
#include "ReversibleObject.hpp"

namespace tempo {
template <typename T> class Solver;
}

namespace tempo::heuristics {
/**
 * @brief Random variable selection strategy
 * @details @copybrief Randomly chooses a variable from the remenaing search
 * variables
 * @note Right now, only binary variables are selected
 */
struct VSIDSHeap {

  template <concepts::scalar T>
  VSIDSHeap(Solver<T> &solver)
      : handlerToken(solver.ClauseAdded.subscribe_handled(
            [this](const auto &arg) { this->updateActivity(arg); })),
        activity(solver.getOptions().vsids_decay) {

    var_heap.resize(solver.boolean.size());
    std::iota(var_heap.begin(), var_heap.end(), 0);
    index.resize(solver.boolean.size(), 0);
    trail.push_back(var_heap.size());

    //@TODO: better initialisation
    for (size_t i{0}; i < var_heap.size(); ++i) {
      std::swap(var_heap[i], var_heap[i + random() % (var_heap.size() - i)]);
      index[var_heap[i]] = static_cast<index_t>(i);
    }

    activity.resize(solver.boolean.size(), impl::ActivityMap::baseIncrement);
  }

  VSIDSHeap(const VSIDSHeap &) = delete;
  VSIDSHeap(VSIDSHeap &&) = delete;
  VSIDSHeap &operator=(const VSIDSHeap &) = delete;
  VSIDSHeap &operator=(VSIDSHeap &&) = delete;
  ~VSIDSHeap() = default;

  /**
   * Heuristic interface
   * @tparam T timing type
   * @param solver solver for which to select the variable
   * @return randomly selected variable
   */
  template <concepts::scalar T>
  auto nextVariable(const Solver<T> &solver) noexcept -> VariableSelection {
    const concepts::same_template<SparseSet> auto &variables =
        solver.getBranch();

    assert(not variables.empty());

    auto n{trail.back()};
    while (static_cast<int>(trail.size()) > solver.level()) {
      trail.pop_back();
    }

    while (n < trail.back()) {
      heap::percolate_up(var_heap.begin(), n, index,
                         [&](const var_t x, const var_t y) {
                           return activity[x] > activity[y];
                         });
      ++n;
    }

    while (static_cast<int>(trail.size()) < solver.level())
      trail.push_back(trail.back());

    assert(checkVars(solver));
    assert(checkHeap(0));

    auto last{trail.back()};
    var_t x;
    do {
      x = pickBest(last);
    } while (not variables.has(x));

    trail.push_back(last);

    assert(checkHeap(0));

    return {x, VariableType::Boolean};
  }

  template <concepts::scalar T> void updateActivity(const Solver<T> &solver) {

    var_buffer.clear();
    for (auto l : solver.lastLearnt()) {
      if (not l.isNumeric()) {
        var_buffer.push_back(l.variable());
      }
    }

    for (auto i : solver.cut.cached_) {
      auto l{solver.getLiteral(i)};
      if (not l.isNumeric())
        var_buffer.push_back(l.variable());
    }

    activity.update(var_buffer);

    for (auto x : var_buffer) {
      if (index[x] < trail.back()) {
        heap::percolate_up(var_heap.begin(), index[x], index,
                           [&](const var_t x, const var_t y) {
                             return activity[x] > activity[y];
                           });
      }
    }

    assert(checkHeap(0));
  }

  template <concepts::scalar T> bool checkVars(const Solver<T> &solver) {
    auto last{trail.back()};
    for (auto x : solver.getBranch()) {
      bool notin{true};
      for (auto v{var_heap.begin()}; notin and v != (var_heap.begin() + last);
           ++v) {
        if (*v == x) {
          notin = false;
        }
      }
      if (notin) {
        std::cout << x << " is not in the heap!\n";
        return false;
      }
    }
    return true;
  }

  bool checkHeap(const int i) {
    auto lc{heap::left(i)};
    auto rc{heap::right(i)};

    auto n{static_cast<int>(trail.back())};

    bool ok_left =
        ((lc >= n) or
         ((activity[var_heap[i]] >= activity[var_heap[lc]]) and checkHeap(lc)));

    if (ok_left) {
      auto ok_right =
          ((rc >= n) or ((activity[var_heap[i]] >= activity[var_heap[rc]]) and
                         checkHeap(rc)));
      return ok_right;
    }

    return false;
  }

  var_t pickBest(size_t &_end_) {
    var_t x{var_heap[0]};
    heap::remove_min(var_heap.begin(), var_heap.begin() + _end_, index,
                     [&](const var_t x, const var_t y) {
                       return activity[x] > activity[y];
                     });
    --_end_;
    return x;
  }

  // the successive sizes of the heap -> decreasing at each call to
  // nextVariable() |trail| should therefore be equal to solver.level(), if not,
  // there has been a backtrack then the values in trail can be used to
  // re-integrate the necessary variables
  std::vector<size_t> trail;

  // position of each variable in the heap (so that we can percolate)
  std::vector<index_t> index;

  std::vector<var_t> var_heap;

  SubscriberHandle handlerToken;

  std::vector<var_t> var_buffer;

  impl::ActivityMap activity;
};
}

#endif // TEMPO_VSIDSHEAP_HPP
