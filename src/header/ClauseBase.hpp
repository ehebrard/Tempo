/************************************************
 * Tempo ClauseBase.hpp
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

#ifndef _TEMPO_CLAUSEBASE_HPP
#define _TEMPO_CLAUSEBASE_HPP

#include <vector>
#include <assert.h>

#include "Clause.hpp"
#include "Failure.hpp"
#include "Global.hpp"
#include "constraints/Constraint.hpp"
#include "util/Heap.hpp"
#include "util/Options.hpp"
#include "util/SparseSet.hpp"
#include "util/SubscribableEvent.hpp"

//#define DBG_WATCHERS

namespace tempo {

template <typename T> class Solver;

//! ClauseBase
/*!
 CNF formula with unit-propagation method (two-watched literals)
*/
template <class T> class ClauseBase : public Constraint<T> {

public:
  /**
   * @name constructors
   */
  //@{
  ClauseBase(Solver<T> &);
  ~ClauseBase() = default;

  // notify that the solver has Boolean variable x
  void newBooleanVar(const var_t x);

  // notify that the solver has numeric variable x
  void newNumericVar(const var_t x);
  //@}

  /**
   * @name accessors
   */
  //@{
  // number of clauses
  size_t size() const;
  // current number of clauses that have been learnt
  size_t numLearnt() const;
  // total number of literals in the clauses
  size_t volume() const;
  // returns the clause with id i
  Clause<T> *operator[](const index_t i);
  // returns the last clause to have been added
  Clause<T> *back();

  // create and add a clause from a list of literals (and returns it) 'learnt'
  // flag distinguished constraints from cuts
  template <typename iter>
  Clause<T> *add(const iter first, const iter last, const bool learnt = false);
  // helpers to handle any type of literal
  bool satisfied(const Literal<T>) const;
  bool falsified(const Literal<T>) const;
  //@}

  /**
   * @name unit propagation
   */
  //@{
  // post the clause base as a constraint (to use only if 'propagate()' is to be
  // used)
  virtual void post(const int idx);

  // notify a change (with the literal and it's variable rank in the scope)
  virtual bool notify(const Literal<T>, const int);

  // set the 'r'-th watcher of clause 'cl' to be its 'i'-th literal
  void set_watcher(const int r, const index_t i, Clause<T> *cl);

  //  // set the 'r'-th watcher of clause 'cl' to be its 'i'-th literal (and
  //  order
  //  // the watch-list of this literal)
  //  void set_watcher_numeric(const int r, const index_t i, Clause<T> *cl);

  // works for both numeric and Boolean, does not order numeric watch-lists
  void unit_propagate(const Literal<T> l);

  // works only for Boolean literals
  void unit_propagate_boolean(const Literal<T> l);

  // works only for numeric literals, orders numeric watch-lists
  void unit_propagate_numeric(const Literal<T> l);

  //  // works for both numeric and Boolean, orders numeric watch-lists
  //  void unit_propagate_beta(const Literal<T> l);

  // unit propagates all the numeric triggers in 'triggered_bounds'
  virtual void propagate();

  // clears the set of numeric triggers ('triggered_bounds')
  void clearTriggers();

  // assigns literal 'l' with explanation 'e'
  void assign(const Literal<T> l, const Explanation<T> &e);
  //@}

  /**
   * @name printing
   */
  //@{
  std::ostream &display(std::ostream &os) const;
  std::ostream &displayWatchStruct(std::ostream &os) const;
  //@}

  /**
   * @name explanation
   */
  //@{
  void xplain(const Literal<T> l, const hint h, std::vector<Literal<T>> &Cl);
  std::ostream &print_reason(std::ostream &, const hint) const;
  int getType() const;
  //@}

  /**
   * @name clause forgetting
   */
  //@{
  // forgets a learnt clause
  void forget(Clause<T> *cl);
  // forgets the last learnt clause (worst if they have been sorted)
  void forget_worst();
  // forgets according to the forgetting policy in solver.getOptions()
  void forget();
  // forgets all learnt clauses
  void forgetAll();
  // literal score based on its semantic
  double looseness(const Literal<T> l);
  // literal score based on its activity
  double inverseActivity(const Literal<T> l);
  // literal score based on both its semantic and its activity
  double loosenessOverActivity(const Literal<T> l);
  //@}

  /**
   * @name debug
   */
  //@{
  Clause<T> *consistent();
  //@}

private:
  Solver<T> &solver;

  SubscriberHandle handlerToken;

  /**
   * @name clause memory store
   */
  //@{
  // a vector containing all the clauses [the id of a clause is its rank in this
  // vector]
  std::vector<Clause<T> *> base;
  // a score vector [same indexing as 'base']
  std::vector<double> score;
  // free indices in the vector base (so that we can remove clauses)
  // - contain the indices in 'base' that are available (no current clause has
  // this id)
  // - at the back of the sparse set are the indices used by learnt clauses
  // - at the front of the sparse set are the indices used by other clauses
  SparseSet<int> free_cl_indices;
  // the watch lists for every literel (watch[BOOLEAN] and watch[NUMERIC])
  std::vector<std::vector<Clause<T> *>> watch[2];
  //@}

  // the bounds literals that need unit-propagation
  SparseSet<var_t> triggered_bounds;
  // helper to unit-propagate ordered watch lists
  std::vector<int> search_stack;

  //    unsigned lazyness{3};
  //    std::vector<unsigned> count;

public:
  /**
   * @name statistics
   */
  //@{
  // number of calls to bound literals propagation
  unsigned long num_prop{0};
  // number of calls to unit-progation
  unsigned long num_up{0};
  // number of misses when unit-propagating bound literals
  unsigned long num_miss{0};
  // total number of literals in the clauses
  size_t total_size{0};
  //@}

  // indexing helper for the watch lists
  static const bool BOOLEAN{false};
  static const bool NUMERIC{true};

  // type helper for the literals -> return l.isNumeric() == NUMERIC is l is
  // numeric
  static bool litType(Literal<T> l) { return l.isNumeric(); }

#ifdef DBG_WATCHERS
  void verifyWatchers(const char *msg) const;
#endif
};


////// NEW CLAUSES
///

template <typename T>
ClauseBase<T>::ClauseBase(Solver<T> &c)
    : solver(c), handlerToken(solver.SearchRestarted.subscribe_handled(
                     [this]() { this->forget(); })) {

  Constraint<T>::priority = Priority::Low;
    }


template <typename T>
void ClauseBase<T>::post(const int idx) {
    Constraint<T>::cons_id = idx;
    if(solver.getOptions().full_up) {
        //    for(var_t x{0}; x<solver.boolean.size(); ++x) {
        //        solver.wake_me_on(solver.boolean.getLiteral(true, x), Constraint<T>::cons_id);
        //        solver.wake_me_on(solver.boolean.getLiteral(false, x), Constraint<T>::cons_id);
        //    }
        for(var_t x{0}; x<solver.numeric.size(); ++x) {
            solver.wake_me_on(lb<T>(x), Constraint<T>::cons_id);
            solver.wake_me_on(ub<T>(x), Constraint<T>::cons_id);
        }
        triggered_bounds.reserve(2 * solver.numeric.size());
    }
}

// propagate the constraint
template <typename T>
void ClauseBase<T>::propagate() {
    
    ++num_prop;
    
    for(auto b : triggered_bounds) {
      auto p{solver.numeric.getLiteral(Literal<T>::sgn(b), Literal<T>::var(b))};
      //      std::cout << " --> unitprop " << solver.pretty(p) << std::endl;
      if (solver.getOptions().order_bound_watch) {
        unit_propagate_numeric(p);
      } else {
        unit_propagate(p);
      }
//        unit_propagate_numeric(p);
    }
    clearTriggers();
    //    std::cout << triggered_bounds << std::endl;
}

template <typename T> void ClauseBase<T>::clearTriggers() {
  triggered_bounds.clear();
}

template <typename T>
bool ClauseBase<T>::notify(const Literal<T> l, const int) {

    assert(l.isNumeric());

    //    if(solver.num_choicepoints - count[l] > lazyness) { //and not
    //    triggered_bounds.has(l)) {
    //
    //        if(triggered_bounds.has(l))
    //        {
    //            std::cout << "weird\n";
    //            exit(1);
    //        }
    //        count[l] = solver.num_choicepoints;
    //
    //        triggered_bounds.add(l);
    //        return true;
    //    }

    if (not triggered_bounds.has(l)) {
      triggered_bounds.add(l);
      return true;
    }

    return false;
}

template<typename T>
size_t ClauseBase<T>::numLearnt() const {
    return free_cl_indices.backsize();
}

template <typename T> size_t ClauseBase<T>::size() const {

  if ((base.size() - free_cl_indices.size()) != (free_cl_indices.backsize() + free_cl_indices.frontsize())) {
    std::cout << "what ?\n";
    exit(1);
  }

  return free_cl_indices.backsize() + free_cl_indices.frontsize();
}

template <typename T> size_t ClauseBase<T>::volume() const {
  return total_size;
}

// template <typename T>
// void ClauseBase<T>::resize(const size_t n, const size_t m) {
//   watch[BOOLEAN].resize(2 * n);
//   watch[NUMERIC].resize(2 * m);
// }

template <typename T> Clause<T> *ClauseBase<T>::operator[](const index_t i) {
  if (i >= free_cl_indices.capacity() or free_cl_indices.has(i))
    return NULL;
  return base[i];
}

template <typename T> Clause<T> *ClauseBase<T>::back() {
  return base[*(free_cl_indices.bbegin())];
}

template <typename T> void ClauseBase<T>::newBooleanVar(const var_t x) {
  watch[BOOLEAN].resize(static_cast<size_t>(2 * x + 2));
}

template <typename T> void ClauseBase<T>::newNumericVar(const var_t x) {
    watch[NUMERIC].resize(static_cast<size_t>(2 * x + 2));
    //    count.resize(static_cast<size_t>(2 * x + 2), 0);
}

template <typename T>
bool ClauseBase<T>::satisfied(const Literal<T> l) const {
  if (l.isNumeric())
    return solver.numeric.satisfied(l);
  else
    return solver.boolean.satisfied(l);
}

template <typename T>
bool ClauseBase<T>::falsified(const Literal<T> l) const {
  if (l.isNumeric())
    return solver.numeric.falsified(l);
  else
    return solver.boolean.falsified(l);
}

// template <typename T>
// void ClauseBase<T>::assign(const Literal<T> l, const Explanation<T> &e)
// {
//   solver.set(l, e);
// }

template <typename T> Clause<T> *ClauseBase<T>::consistent() {
  //    for(var_t x{0}; x<solver.boolean.size(); ++x) {
  //        for(auto cl : watch[BOOLEAN][Literal<T>(true,x)]) {
  //            assert
  //        }
  //    }

  for (auto i{free_cl_indices.fbegin()}; i != free_cl_indices.fend(); ++i) {
    auto cl{base[*i]};
    auto w0{cl->watched(0)};
    auto w1{cl->watched(1)};
    if (satisfied(w0))
      continue;
    if (satisfied(w1))
      continue;
    if (falsified(w0)) {
      if (not satisfied(w1))
        return cl;
      for (auto p : *cl) {
        if (p != w0 and p != w1 and not falsified(p)) {
          return cl;
        }
      }
    } else if (falsified(w1)) {
      if (not satisfied(w0))
        return cl;
      for (auto p : *cl) {
        if (p != w0 and p != w1 and not falsified(p)) {
          return cl;
        }
      }
    }
  }

  for (auto i{free_cl_indices.bbegin()}; i != free_cl_indices.bend(); ++i) {
    auto cl{base[*i]};
    auto w0{cl->watched(0)};
    auto w1{cl->watched(1)};
    if (satisfied(w0))
      continue;
    if (satisfied(w1))
      continue;
    if (falsified(w0)) {
      if (not satisfied(w1))
        return cl;
      for (auto p : *cl) {
        if (p != w0 and p != w1 and not falsified(p)) {
          return cl;
        }
      }
    } else if (falsified(w1)) {
      if (not satisfied(w0))
        return cl;
      for (auto p : *cl) {
        if (p != w0 and p != w1 and not falsified(p)) {
          return cl;
        }
      }
    }
  }

  //    std::cout << "ok\n";
  return NULL;
}

// template <typename T>
// void ClauseBase<T>::unit_propagate_beta(const Literal<T> l) {
//
//   ++num_up;
//
//#ifdef DBG_WATCHERS
//   verifyWatchers("before UP");
//#endif
//
//#ifdef DBG_TRACE
//   if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
//     std::cout << "unit propagate true lit " << l << "\n";
//   }
//#endif
//
//   if (l.isNumeric()) {
//     unit_propagate_numeric(l);
//   } else {
//     unit_propagate_boolean(l);
//   }
////    unit_propagate_generic(l);
//
//#ifdef DBG_WATCHERS
//  verifyWatchers("after UP");
//#endif
//}

template <typename T>
void ClauseBase<T>::unit_propagate(const Literal<T> l) {
    
    ++num_up;

  #ifdef DBG_WATCHERS
    verifyWatchers("before UP");
  #endif

  auto lt{litType(l)};

  // watch[lt][l] contains all clauses watched by ~l.
  // unless l is a numeric literal x <= k, in which case watch[lt][l] contains
  // clauses watched by some literal -x <= v therefore the trigger is "real"
  // only if k+v < 0
  for (auto c{watch[lt][l].rbegin()}; c != watch[lt][l].rend(); ++c) {

    auto cl{*c};

#ifdef DBG_TRACE
    if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
      std::cout << " watched by " << *cl << std::endl;
    }
#endif

    bool watch_rank{cl->watch_rank(l)};
    index_t idx{cl->watched_index[watch_rank]};
    Literal<T> other{cl->watched(1 - watch_rank)};
    Literal<T> c_lit{(*cl)[idx]};

    assert(c_lit.sameVariable(l));
    assert(c_lit.sign() != l.sign());

    if (lt == NUMERIC) {
      if (c_lit.value() + l.value() >= 0) {

#ifdef DBG_TRACE
        if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                    std::cout << " false trigger (" << c_lit << " is not falsified)\n";
        }
#endif

        ++num_miss;
        continue;
      } else {

#ifdef DBG_TRACE
        if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
          std::cout << " true trigger (" << c_lit << " and " << l
                    << " are contradictory)\n";
        }
#endif
      }
    }

    if (satisfied(other)) {

#ifdef DBG_TRACE
      if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
        std::cout << ": satisfied by " << other << std::endl;
      }
#endif

      continue;
    }

    index_t i{idx};

#ifdef DBG_TRACE
    if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
      std::cout << "search another literal to watch";
    }
#endif
    while (true) {

      if (++i == cl->size())
        i = 0;
      if (i == idx)
        break;

      // look for a replacement
      auto p = (*cl)[i];

#ifdef DBG_TRACE
      if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
        std::cout << "  " << p;
      }
#endif

      if (p != other) {
        if (not falsified(p)) {

#ifdef DBG_TRACE
          if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
            std::cout << ": replace by " << p << " (" << i << ") as "
                      << watch_rank << "-th watcher ";
          }
#endif

          set_watcher(watch_rank, i, cl);

#ifdef DBG_TRACE
          if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
            std::cout << " and rm from " << l << "'s watches\n";
          }
#endif

          // remove clause from l's watch list
          swap(*c, *watch[lt][l].rbegin());
          watch[lt][l].pop_back();

          break;
        }
      }
#ifdef DBG_TRACE
      else if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
        std::cout << "*";
      }
#endif
    }

    if (i == idx) {

#ifdef DBG_TRACE
      if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
        std::cout << ": new unit " << other << std::endl;
      }
#endif

#ifdef DBG_WATCHERS
      verifyWatchers("at assign");
#endif

      solver.set(other, {this, cl->id});

      assert(cl == base[cl->id]);
    }
  }
    
    
#ifdef DBG_WATCHERS
  verifyWatchers("after UP");
#endif
}

template <typename T>
void ClauseBase<T>::unit_propagate_numeric(const Literal<T> l) {

  // watch[lt][l] contains all clauses watched by ~l.
  // unless l is a numeric literal x <= k, in which case watch[lt][l] contains
  // clauses watched by some literal -x <= v therefore the trigger is "real"
  // only if k+v < 0

  //  std::vector<int> search_stack;

  if (watch[NUMERIC][l].empty())
    return;

  search_stack.clear();

  search_stack.push_back(0);
  while (not search_stack.empty()) {

#ifdef DBG_TRACE
    if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
      for (index_t i{0}; i < watch[NUMERIC][l].size(); ++i) {
        std::cout << std::setw(3) << i << " "
                  << watch[NUMERIC][l][i]->watched(
                         watch[NUMERIC][l][i]->watch_rank(l))
                  << " " << *(watch[NUMERIC][l][i]) << std::endl;
      }
      std::cout << "stack:";
      for (auto i : search_stack) {
        std::cout << " " << i;
      }
      std::cout << std::endl;
    }
#endif

    auto k{search_stack.back()};

    if (static_cast<size_t>(k) >= watch[NUMERIC][l].size()) {
      search_stack.pop_back();
      continue;
    }

    auto cl{watch[NUMERIC][l][k]};

    bool watch_rank{cl->watch_rank(l)};
    index_t idx{cl->watched_index[watch_rank]};
    Literal<T> other{cl->watched(1 - watch_rank)};
    Literal<T> c_lit{(*cl)[idx]};

#ifdef DBG_TRACE
    if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
      std::cout << "explore clause " << *cl << " (" << c_lit << ")\n";
    }
#endif

    if (c_lit.value() + l.value() >= 0) {

#ifdef DBG_TRACE
      if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
        std::cout << " false trigger (" << c_lit << " is not falsified)\n";
      }
#endif

      ++num_miss;

      // false trigger, all descendant in the heap are also false trigger
      search_stack.pop_back();
      continue;
    }

    if (satisfied(other)) {

#ifdef DBG_TRACE
      if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
        std::cout << ": satisfied by " << other << std::endl;
      }
#endif

      // no unit from this clause, but descendant are not necessarily false
      // triggers
      search_stack.pop_back();
      auto child{heap::left(k)};
      int end_list{static_cast<int>(watch[NUMERIC][l].size())};
      if (child < end_list) {
        search_stack.push_back(child);
        child = heap::right(k);
        if (child < end_list)
          search_stack.push_back(child);
      }
      continue;
    }

    index_t i{idx};

#ifdef DBG_TRACE
    if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
      std::cout << "search another literal to watch";
    }
#endif
    while (true) {

      if (++i == cl->size())
        i = 0;
      if (i == idx)
        break;

      // look for a replacement
      auto p = (*cl)[i];

#ifdef DBG_TRACE
      if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
        std::cout << "  " << p;
      }
#endif

      if (p != other) {
        if (not falsified(p)) {

#ifdef DBG_TRACE
          if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
            std::cout << ": replace by " << p << " (" << i << ") as "
                      << watch_rank << "-th watcher ";
          }
#endif

          //            if(p.isNumeric())
          //                set_watcher_numeric(watch_rank, i, cl);
          //            else
          set_watcher(watch_rank, i, cl);

#ifdef DBG_TRACE
          if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
            std::cout << " and rm from " << l << "'s watches\n";
          }
#endif

          // remove clause from l's watch list
          swap(*(watch[NUMERIC][l].begin() + k), *watch[NUMERIC][l].rbegin());
          watch[NUMERIC][l].pop_back();

          if (static_cast<size_t>(k) < watch[NUMERIC][l].size()) {
            //                std::cout << "precolated down from cl[" << k << "]
            //                = " << *watch[NUMERIC][l][k] << std::endl;

            heap::percolate_down(
                watch[NUMERIC][l].begin(), watch[NUMERIC][l].end(), k,
                [l](const Clause<T> *c1, const Clause<T> *c2) {
                  return c1->watched(c1->watch_rank(l)).value() <
                         c2->watched(c2->watch_rank(l)).value();
                });
          }

          // resume search from the same point in the heap, no need to change
          // anything
          break;
        }
      }
#ifdef DBG_TRACE
      else if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
        std::cout << "*";
      }
#endif
    }

    if (i == idx) {

#ifdef DBG_TRACE
      if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
        std::cout << ": new unit " << other << std::endl;
      }
#endif

#ifdef DBG_WATCHERS
      verifyWatchers("at assign");
#endif

      solver.set(other, {this, cl->id});
      // there was pruning, descendant might trigger
      search_stack.pop_back();
      auto child{heap::left(k)};
      int end_list{static_cast<int>(watch[NUMERIC][l].size())};
      if (child < end_list) {
        search_stack.push_back(child);
        child = heap::right(k);
        if (child < end_list)
          search_stack.push_back(child);
      }

      assert(cl == base[cl->id]);
    }
  }
}

// template <typename T>
// void ClauseBase<T>::unit_propagate_numeric(const Literal<T> l) {
//
//
//   // watch[lt][l] contains all clauses watched by ~l.
//   // unless l is a numeric literal x <= k, in which case watch[lt][l]
//   contains
//   // clauses watched by some literal -x <= v therefore the trigger is "real"
//   // only if k+v < 0
//
////    std::vector<index_t> search_stack;
////
////    search_stack.
//
//  for (auto c{watch[NUMERIC][l].rbegin()}; c != watch[NUMERIC][l].rend(); ++c)
//  {
//
//    auto cl{*c};
//
//#ifdef DBG_TRACE
//    if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
//      std::cout << " watched by " << *cl << std::endl;
//    }
//#endif
//
//    bool watch_rank{cl->watch_rank(l)};
//    index_t idx{cl->watched_index[watch_rank]};
//    Literal<T> other{cl->watched(1 - watch_rank)};
//    Literal<T> c_lit{(*cl)[idx]};
//
////      std::cout << " self = " << c_lit << std::endl;
////      std::cout << "other = " << other << std::endl;
//
//      assert(c_lit.sameVariable(l));
//      assert(c_lit.sign() != l.sign());
//
////      if(l.sign() == c_lit.sign()) {
////          std::cout << l << " is watched by " << *cl << " @" << idx <<
/// std::endl; /          exit(1); /      }
//
////      std::cout << "==> " << c_lit << std::endl;
//
//
//          if(c_lit.value() + l.value() >= 0) {
//
//#ifdef DBG_TRACE
//              if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
//                  std::cout << " false trigger (" << c_lit << " is not
//                  falsified)\n";
//              }
//#endif
//
//              continue;
//          } else {
//
//#ifdef DBG_TRACE
//              if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
//                  std::cout << " true trigger (" << c_lit << " and " << l << "
//                  are contradictory)\n";
//              }
//#endif
//
//          }
//
//    if (satisfied(other)) {
//
//#ifdef DBG_TRACE
//      if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
//        std::cout << ": satisfied by " << other << std::endl;
//      }
//#endif
//
//      continue;
//    }
//
//    index_t i{idx};
//
//#ifdef DBG_TRACE
//      if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
//        std::cout << "search another literal to watch" ;
//      }
//#endif
//    while (true) {
//
//      if (++i == cl->size())
//        i = 0;
//      if (i == idx)
//        break;
//
//      // look for a replacement
//      auto p = (*cl)[i];
//
//#ifdef DBG_TRACE
//      if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
//        std::cout << "  " << p;
//      }
//#endif
//
//      if (p != other) {
//        if (not falsified(p)) {
//
//#ifdef DBG_TRACE
//          if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
//            std::cout << ": replace by " << p << " (" << i << ") as "
//                      << watch_rank << "-th watcher ";
//          }
//#endif
//
//          set_watcher(watch_rank, i, cl);
//
//#ifdef DBG_TRACE
//          if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
//            std::cout << " and rm from " << l << "'s watches\n";
//          }
//#endif
//
//          // remove clause from l's watch list
//          swap(*c, *watch[NUMERIC][l].rbegin());
//          watch[NUMERIC][l].pop_back();
//            heap::percolate_down(watch[NUMERIC][l].begin(),
//            watch[NUMERIC][l].end(), 0, [l](const Clause<T>* c1, const
//            Clause<T>* c2) {return c1->watched(c1->watch_rank(l)).value() >
//            c2->watched(c2->watch_rank(l)).value(); });
//
//          break;
//        }
//      }
//#ifdef DBG_TRACE
//      else if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
//        std::cout << "*";
//      }
//#endif
//    }
//
//    if (i == idx) {
//
//#ifdef DBG_TRACE
//      if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
//        std::cout << ": new unit " << other << std::endl;
//      }
//#endif
//
//#ifdef DBG_WATCHERS
//      verifyWatchers("at assign");
//#endif
//
//      solver.set(other, {this, cl->id});
//
//      assert(cl == base[cl->id]);
//    }
//  }
//
//}

template <typename T>
void ClauseBase<T>::unit_propagate_boolean(const Literal<T> l) {

  // watch[lt][l] contains all clauses watched by ~l.
  // unless l is a numeric literal x <= k, in which case watch[lt][l] contains
  // clauses watched by some literal -x <= v therefore the trigger is "real"
  // only if k+v < 0
  for (auto c{watch[BOOLEAN][l].rbegin()}; c != watch[BOOLEAN][l].rend(); ++c) {

    auto cl{*c};

#ifdef DBG_TRACE
    if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
      std::cout << " watched by " << *cl << std::endl;
    }
#endif

    bool watch_rank{cl->watch_rank(l)};
    index_t idx{cl->watched_index[watch_rank]};
    Literal<T> other{cl->watched(1 - watch_rank)};
    //    Literal<T> c_lit{(*cl)[idx]};

    if (satisfied(other)) {

#ifdef DBG_TRACE
      if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
        std::cout << ": satisfied by " << other << std::endl;
      }
#endif

      continue;
    }

    index_t i{idx};

#ifdef DBG_TRACE
    if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
      std::cout << "search another literal to watch";
    }
#endif
    while (true) {

      if (++i == cl->size())
        i = 0;
      if (i == idx)
        break;

      // look for a replacement
      auto p = (*cl)[i];

#ifdef DBG_TRACE
      if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
        std::cout << "  " << p;
      }
#endif

      if (p != other) {
        if (not falsified(p)) {

#ifdef DBG_TRACE
          if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
            std::cout << ": replace by " << p << " (" << i << ") as "
                      << watch_rank << "-th watcher ";
          }
#endif

          set_watcher(watch_rank, i, cl);

#ifdef DBG_TRACE
          if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
            std::cout << " and rm from " << l << "'s watches\n";
          }
#endif

          // remove clause from l's watch list
          swap(*c, *watch[BOOLEAN][l].rbegin());
          watch[BOOLEAN][l].pop_back();

          break;
        }
      }
#ifdef DBG_TRACE
      else if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
        std::cout << "*";
      }
#endif
    }

    if (i == idx) {

#ifdef DBG_TRACE
      if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
        std::cout << ": new unit " << other << std::endl;
      }
#endif

#ifdef DBG_WATCHERS
      verifyWatchers("at assign");
#endif

      solver.set(other, {this, cl->id});

      assert(cl == base[cl->id]);
    }
  }
}

template <typename T>
void ClauseBase<T>::set_watcher(const int r, const index_t i, Clause<T> *cl) {
  cl->watched_index[r] = i;
  Literal<T> l{~((*cl)[i])};
  if (l.isNumeric()) {
    watch[NUMERIC][l].push_back(cl);
    if (solver.getOptions().order_bound_watch) {
      heap::percolate_up(watch[NUMERIC][l].begin(),
                         watch[NUMERIC][l].size() - 1,
                         [l](const Clause<T> *c1, const Clause<T> *c2) {
                           return c1->watched(c1->watch_rank(l)).value() <
                                  c2->watched(c2->watch_rank(l)).value();
                         });
    }
  } else {
    watch[BOOLEAN][l].push_back(cl);
  }
}

// template <typename T>
// void ClauseBase<T>::set_watcher_numeric(const int r, const index_t i,
//                                    Clause<T> *cl) {
//
//   cl->watched_index[r] = i;
//   Literal<T> l{~((*cl)[i])};
//
//     assert(l.isNumeric());
//
//     auto j{watch[NUMERIC][l].size()};
//     watch[NUMERIC][l].push_back(cl);
//     heap::percolate_up(watch[NUMERIC][l].begin(), j,
//                        [l](const Clause<T> *c1, const Clause<T> *c2) {
//                          return c1->watched(c1->watch_rank(l)).value() <
//                                 c2->watched(c2->watch_rank(l)).value();
//                        });
// }

template <typename T> void ClauseBase<T>::forget(Clause<T> *cl) {
  for (auto r{0}; r < 2; ++r) {
    auto wl{~(cl->watched(r))};
    auto wt{litType(wl)};

    for (auto c{watch[wt][wl].begin()}; c != watch[wt][wl].end(); ++c)
      if (*c == cl) {
        swap(*c, watch[wt][wl].back());
        watch[wt][wl].pop_back();
        break;
      }
  }

  free_cl_indices.add(cl->id);

  total_size -= cl->size();

  cl->clear();
}

/// | |  1 2  6 5 3 4

template <typename T> void ClauseBase<T>::forget_worst() {
    
//    std::cout << free_cl_indices << std::endl;
//    
//    std::cout << (*(free_cl_indices.bbegin())) << std::endl;
//    
//    std::cout << base.size() << std::endl;
//    
//    std::cout << *(base[*(free_cl_indices.bbegin())]) << std::endl;
////    exit(1);
    
  forget(base[*(free_cl_indices.bbegin())]);
}

template <typename T>
double ClauseBase<T>::loosenessOverActivity(const Literal<T> l) {
  return solver.looseness(l) / solver.getActivityMap()->get(l, solver);
}

template <typename T>
double ClauseBase<T>::inverseActivity(const Literal<T> l) {
  return 1.0 / solver.getActivityMap()->get(l, solver);
}

template <typename T> double ClauseBase<T>::looseness(const Literal<T> l) {
  return solver.looseness(l);
}

template <typename T> void ClauseBase<T>::forgetAll() {
  while (free_cl_indices.backsize() > 0) {
    forget_worst();
  }
}

template <typename T> void ClauseBase<T>::forget() {

  if (size() < 1000)
    return;

#ifdef DBG_WATCHERS
  verifyWatchers("before forget");
#endif

  //    std::cout << "option: " << solver.getOptions().forget_strategy <<
  //    std::endl;

  if (solver.getOptions().forget_strategy == Options::LiteralScore::Looseness or
      (solver.getActivityMap() == NULL and
       solver.getOptions().forget_strategy != Options::LiteralScore::Size)) {
    for (auto idx{free_cl_indices.bbegin()}; idx != free_cl_indices.bend();
         ++idx) {
      score[*idx] = 0;
        
//        std::cout << "compute (1) score for " << *base[*idx] << std::endl;
        
      for (auto l : *base[*idx]) {
        score[*idx] += looseness(l);
      }
    }
  } else if (solver.getOptions().forget_strategy ==
             Options::LiteralScore::LoosenessOverActivity) {

    //  std::cout << "forget with l/a\n";

    for (auto idx{free_cl_indices.bbegin()}; idx != free_cl_indices.bend();
         ++idx) {
      score[*idx] = 0;
        
//        std::cout << "compute (2) score for " << *base[*idx] << std::endl;
        
      for (auto l : *base[*idx]) {
        score[*idx] += loosenessOverActivity(l);
      }
    }
  } else if (solver.getOptions().forget_strategy ==
             Options::LiteralScore::Activity) {
    for (auto idx{free_cl_indices.bbegin()}; idx != free_cl_indices.bend();
         ++idx) {
      score[*idx] = 0;
        
//        std::cout << "compute (2) score for " << *base[*idx] << std::endl;
        
      for (auto l : *base[*idx]) {
        score[*idx] += inverseActivity(l);
      }
    }
  }

  if (solver.getOptions().forget_strategy == Options::LiteralScore::Size) {
    std::sort(free_cl_indices.bbegin(), free_cl_indices.bend(),
              [&](const int i, const int j) {
                return (base[i]->size() > base[j]->size());
              });
  } else {
    std::sort(free_cl_indices.bbegin(), free_cl_indices.bend(),
              [&](const int i, const int j) {
                return (score[base[i]->id] > score[base[j]->id]);
              });
  }

  free_cl_indices.re_index(free_cl_indices.bbegin(), free_cl_indices.bend());

  //
  //    std::cout << "ordering, from worst to best:\n";
  //    for(auto i{free_cl_indices.bbegin()}; i!=free_cl_indices.bend(); ++i) {
  //        std::cout << *base[*i] << std::endl;
  //    }
  //
  auto target_size = static_cast<size_t>(
      static_cast<double>(numLearnt()) * (1.0 - solver.getOptions().forgetfulness));

  while (numLearnt() > target_size) {
    forget_worst();
  }
  //
  //    std::cout << "\n===>\n";
  //    for(auto i{free_cl_indices.bbegin()}; i!=free_cl_indices.bend(); ++i) {
  //        std::cout << *base[*i] << std::endl;
  //    }

#ifdef DBG_WATCHERS
  verifyWatchers("after forget");
#endif
}

template <typename T>
template <typename iter>
Clause<T> *ClauseBase<T>::add(const iter first, const iter last,
                                    const bool learnt) {

#ifdef DBG_WATCHERS
  verifyWatchers("before add");
#endif

  //  std::cout << "add clause (" << (learnt ? "learnt" : "base" ) << ")";
  //  for (auto l{first}; l != last; ++l) {
  //    std::cout << " " << *l;
  //  }
  //  std::cout << std::endl;

  Clause<T> *c{NULL};

  assert(first != last);

  if (first + 1 == last) {
    //    assign(*first, Constant::NewNoReason<T>);
    solver.set(*first, Constant::GroundFact<T>);
  } else {
    if (not free_cl_indices.empty()) {
      int id{free_cl_indices.back()};
      if (learnt)
        free_cl_indices.remove_back(id);
      else
        free_cl_indices.remove_front(id);
      c = base[id];
      assert(c->empty());
    } else {
      int id{static_cast<int>(base.size())};
      c = new Clause<T>(id);
      base.push_back(c);
      score.push_back(0);
      free_cl_indices.reserve(base.size());
        
//    std::cout << free_cl_indices << std::endl;
        
      if (not learnt) {
        free_cl_indices.add(id);
        free_cl_indices.remove_front(id);
      } 
        
//    std::cout << free_cl_indices << std::endl;
        
//      else {
//          free_cl_indices.add(id);
//          free_cl_indices.remove_back(id);
//      }
    }

    for (auto l{first}; l != last; ++l) {
      c->push_back(*l);
    }
      
    set_watcher(0, 0, c);
    set_watcher(1, 1, c);

    total_size += c->size();

      if(learnt) {
          Literal<T> l{(*c)[0]};
          //          assign(l, {this, c->id});
          solver.set(l, {this, c->id});
      }
  }

#ifdef DBG_WATCHERS
  verifyWatchers("after add");
#endif

  return c;
}

//template <typename T>
//std::ostream &ClauseBase<T>::displayClause(std::ostream &os,
//                                              const Clause<T> *cl) const {
//  //    os << "[" << cl->watch_index(0) << "|" << cl->watch_index(1) << "]";
//
//  //    assert(not cl->empty());
//
//  //    if (cl->empty()) {
//  //        os << "(" << cl->id <<":)";
//  //    } else {
//  //
//  //        os << "(";
//  //        os << cl->id <<":";
//  ////        os << "[" << prettyLiteral((*cl)[0]) << "]";
//  ////        if (0 == cl->watch_index(0) or 0 == cl->watch_index(1))
//  ////            os << "*";
//  ////        for (size_t i{1}; i < cl->size(); ++i) {
//  ////            os << " or [" << prettyLiteral((*cl)[i]) << "]";
//  ////            if (i == cl->watch_index(0) or i == cl->watch_index(1))
//  ////                os << "*";
//  ////        }
//  //        if(LTYPE((*cl)[0]) == BOUND_LIT)
//  //            os << FROM_GEN((*cl)[0]);
//  //        else
//  //            os << ".";
//  //        for (size_t i{1}; i < cl->size(); ++i) {
//  //            if(LTYPE((*cl)[i]) == BOUND_LIT)
//  //                os << " " << FROM_GEN((*cl)[i]);
//  //            else
//  //                os << " .";
//  //        }
//  //        os << ")";
//  //    }
//
//  if (cl->empty()) {
//    os << "()";
//  } else {
//    os << "(";
//    os << "[" << (*cl)[0] << "]";
//    if (0 == cl->watch_index(0) or 0 == cl->watch_index(1))
//      os << "*";
//    for (size_t i{1}; i < cl->size(); ++i) {
//      os << " or [" << (*cl)[i] << "]";
//      if (i == cl->watch_index(0) or i == cl->watch_index(1))
//        os << "*";
//    }
//    os << ")";
//  }
//  return os;
//}

template <typename T>
void ClauseBase<T>::xplain(const Literal<T> l, const hint h,
                              std::vector<Literal<T>> &Cl) {

  //    int i{0};
  //    for(auto cl : base) {
  //        assert(cl->id == i);
  //        ++i;
  //    }

  Clause<T> &reason(*(base[h]));

  //    assert(reason.id == h);

  if (l == Solver<T>::Contradiction) {

    //      std::cout << "explain contradiction\n";

    for (auto p : reason) {
      //          std::cout << " > " << p << std::endl;
      Cl.push_back(~p);
    }
  } else {

    //      std::cout << "explain " << l << std::endl;

    for (auto p : reason)
      if (not l.sameVariable(p) or l.sign() != p.sign()) {
        Cl.push_back(~p);
      }
  }
}

template <typename T>
std::ostream &ClauseBase<T>::print_reason(std::ostream &os,
                                             const hint h) const {
//  return displayClause(os, base[h]);
    os << *(base[h]);
    return os;
}

template <typename T> int ClauseBase<T>::getType() const {
  return CLAUSEEXPL;
}

template <typename T>
std::ostream &ClauseBase<T>::display(std::ostream &os) const {
  os << "clause base";
  return os;
}

template <typename T>
std::ostream &ClauseBase<T>::displayWatchStruct(
    std::ostream &os
    //                                            , const std::function<int>&
    //                                            f=TODIMACS
) const {
  //    os << "base:";
  //    for(auto cl : base)
  //        if(not free_cl_indices.has(cl->id))
  //            os << " " << *cl << std::endl;
  //        else
  //            os << "(deleted clause " << cl->id << ")\n";

  for (size_t x{0}; x < solver.boolean.size(); ++x) {
    auto l{solver.boolean.getLiteral(true, x)};
    if (not watch[BOOLEAN][l].empty()) {
      os << l << " is watched in";
      for (auto cl : watch[BOOLEAN][l]) {
        os << " " << *cl;
      }
      os << std::endl;
    }
    if (not watch[BOOLEAN][~l].empty()) {
      os << ~l << " is watched in";
      for (auto cl : watch[BOOLEAN][~l]) {
        os << " " << *cl;
      }
      os << std::endl;
    }
  }
  for (size_t x{0}; x < solver.numeric.size(); ++x) {
    auto l = lb<T>(x);
    if (not watch[NUMERIC][l].empty()) {
      os << l << " is watched in";
      for (auto cl : watch[NUMERIC][l]) {
        os << " " << *cl;
      }
      os << std::endl;
    }
    if (not watch[NUMERIC][~l].empty()) {
      os << ~l << " is watched in";
      for (auto cl : watch[NUMERIC][~l]) {
        os << " " << *cl;
      }
      os << std::endl;
    }
  }

  return os;
}

#ifdef DBG_WATCHERS
template <typename T>
void ClauseBase<T>::verifyWatchers(const char *msg) const {
    
    //    std::vector<>
    
    int i{0};
    for (auto cl : base) {
        if (cl->id != i) {
            std::cout << msg << "indexing error" << std::endl;
            exit(1);
        }
        ++i;
    }
    
    size_t num_watchers{0};
    for (size_t x{0}; x < solver.boolean.size(); ++x) {
        Literal<T> l{solver.boolean.getLiteral(true,x)};
        num_watchers += watch[BOOLEAN][l].size();
        if (not watch[BOOLEAN][l].empty()) {
            for (auto cl : watch[BOOLEAN][l]) {
                
                if (cl->size() < 2) {
                    std::cout << msg << ": error empty clause watching "
                    << l << std::endl;
                    exit(1);
                }
                
                if (free_cl_indices.has(cl->id)) {
                    std::cout << *cl << "'s id is free " << std::endl;
                    exit(1);
                }

                if (cl->watched(0) != ~l and cl->watched(1) != ~l) {
                  std::cout << msg << ": error on clause " << cl->id << " -- "
                            << *cl << " watching " << ~l << std::endl;
                  exit(1);
                }
            }
        }
        if (not watch[BOOLEAN][~l].empty()) {
            num_watchers += watch[BOOLEAN][~l].size();
            for (auto cl : watch[BOOLEAN][~l]) {
                if (cl->size() < 2) {
                    std::cout << msg << ": error empty clause watching "
                    << l << std::endl;
                    exit(1);
                }
                
                if (free_cl_indices.has(cl->id)) {
                    std::cout << *cl << "'s id is free " << std::endl;
                    exit(1);
                }

                if (cl->watched(0) != l and cl->watched(1) != l) {
                  std::cout << msg << ": error on clause " << cl->id << " -- "
                            << *cl << " watching " << l << std::endl;
                  exit(1);
                }
            }
        }
    }
    for (size_t x{0}; x < solver.numeric.size(); ++x) {
      Literal<T> l{lb<T>(x)};
      num_watchers += watch[NUMERIC][l].size();
      if (not watch[NUMERIC][l].empty()) {
        for (auto cl : watch[NUMERIC][l]) {
          if (cl->size() < 2) {
            std::cout << msg << ": error empty clause watching " << l
                      << std::endl;
            exit(1);
          }

          if (free_cl_indices.has(cl->id)) {
            std::cout << *cl << "'s id is free " << std::endl;
            exit(1);
          }

          if (not((cl->watched(0).sameVariable(l) and cl->watched(0).sign() != l.sign())  or
                  (cl->watched(1).sameVariable(l) and cl->watched(1).sign() != l.sign()))) {
            std::cout << msg << ": error on clause " << cl->id << " -- " << *cl
                      << " watching " << l << " (" << info_t(l) << ")" << std::endl;
            exit(1);
          }
        }
      }
        auto p{ub<T>(x)};
      if (not watch[NUMERIC][p].empty()) {
        num_watchers += watch[NUMERIC][p].size();
        for (auto cl : watch[NUMERIC][p]) {
          if (cl->size() < 2) {
            std::cout << msg << ": error empty clause watching " << p
                      << std::endl;
            exit(1);
          }

          if (free_cl_indices.has(cl->id)) {
            std::cout << *cl << "'s id is free " << std::endl;
            exit(1);
          }

          if (not((cl->watched(0).sameVariable(p) and cl->watched(0).sign() != p.sign())  or
                  (cl->watched(1).sameVariable(p) and cl->watched(1).sign() != p.sign()))) {
            std::cout << msg << ": error on clause " << cl->id << " -- " << *cl
                      << " watching " << p << " (" << info_t(p) << ")" << std::endl;
            exit(1);
          }
        }
      }
    }

    if (num_watchers != 2 * size()) {
        std::cout << msg << ": wrong number of watchers !\n";
        exit(1);
    }
    //    std::cout << msg << ": ok\n";
}
#endif

template <typename T>
std::ostream &operator<<(std::ostream &os, const ClauseBase<T> &x) {
  return x.display(os);
}

}

#endif // _TEMPO_CLAUSEBASE_HPP

