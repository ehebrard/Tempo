

#ifndef _TEMPO_CLAUSEBASE_HPP
#define _TEMPO_CLAUSEBASE_HPP

#include <vector>
#include <assert.h>

#include "Clause.hpp"
#include "Failure.hpp"
#include "Global.hpp"
#include "constraints/Constraint.hpp"
#include "util/Options.hpp"
#include "util/SparseSet.hpp"
#include "util/SubscribableEvent.hpp"

//#define DBG_WATCHERS

namespace tempo {

template <typename T> class Solver;

template <class T> class ClauseBase : public Constraint<T> {

public:
  //    std::vector<Clause *> base;
  //    std::vector<Clause *> learnt;

  ClauseBase(Solver<T> &);
  ~ClauseBase() = default;

  //  void resize(const size_t n, const size_t m);
  void newBooleanVar(const var_t x);
  void newNumericVar(const var_t x);
  size_t size() const;
    size_t numLearnt() const;

  size_t volume() const;

  //    void unit_propagate(const var l);
  void set_watcher(const int r, const index_t i, Clause<T> *cl);

  // "learn" makes several assumptions about the added clause, too lazy to write
  // the robust version yet
  //    template <typename iter>
  //    void add(const iter first, const iter last);

  template <typename iter>
  Clause<T> *add(const iter first, const iter last,
                    const bool learnt = false);

  void unit_propagate(const Literal<T> l);

  void clearTriggers();

  void assign(const Literal<T> l, const Explanation<T> &e);

  std::ostream &display(std::ostream &os) const;
  std::ostream &displayWatchStruct(std::ostream &os) const;

  //  std::ostream &displayClause(std::ostream &os, const Clause<T> *c)
  //  const;

  bool satisfied(const Literal<T>) const;
  bool falsified(const Literal<T>) const;

  //  lit newNegLiteral(const lit l);
  //  lit getReasonLit(const lit l) const;

  //  bool sameVar(const Literal<T>, const Literal<T>) const;
  //    lit getVarInClause(const lit l) const;
  //    lit getVarInTrigger(const lit l) const;

  // explanation stuff
  void xplain(const Literal<T> l, const hint h, std::vector<Literal<T>> &Cl);
  std::ostream &print_reason(std::ostream &, const hint) const;
  int getType() const;

  //  void notifyRemoval(const Literal<T> l);
  void forget(Clause<T> *cl);
  void forget_worst();
  void forget();
  void forgetAll();
  double looseness(const Literal<T> l);
  double inverseActivity(const Literal<T> l);
  double loosenessOverActivity(const Literal<T> l);

  Clause<T> *operator[](const index_t i) {

    //        std::cout << i << std::endl << free_cl_indices << std::endl;

    if (i >= free_cl_indices.capacity() or free_cl_indices.has(i))
      return NULL;
    return base[i];
  }

  Clause<T> *back() {
    //
    //    std::cout << base.size() << " / "
    //              << static_cast<size_t>(free_cl_indices.bbegin() -
    //                                     free_cl_indices.fbegin())
    //              << std::endl;

    return base[*(free_cl_indices.bbegin())];
  }

  //    std::string prettyLiteral(const genlit l) const;

    //
    virtual void post(const int idx);
    // propagate the constraint
    virtual void propagate();
    // notify a change (with the literal and it's variable rank in the scope)
    virtual bool notify(const Literal<T>, const int);

    Clause<T> *consistent();

  private:
    Solver<T> &solver;

    SubscriberHandle handlerToken;

    /// Clause memory store ///
    /// The set of clauses
    std::vector<Clause<T> *> base;
    std::vector<double> score;
    /// free indices in the vector base (so that we can remove clauses)
    SparseSet<int> free_cl_indices;
    // clauses, watch struct watch[BOUND_LIT] for bound literals,
    // watch[EDGE_LIT] for edges
    std::vector<std::vector<Clause<T> *>> watch[2];
    ////////////////////////////////////////

    //    /// Literal memory store ///
    //    // store the set of bound constraints appearing in a literal
    //    std::vector<BoundConstraint<T>> constraint;
    //    // the number of clauses in which this bound-constraint appears
    //    std::vector<int> cardinality;
    //    // free indices in the vector constraints (so that we can remove
    //    constraints when they are not used) SparseSet<int> free_lit_indices;
    //    // for every event-lit, the set of (pointers to) constraints ordered
    //    by bound std::vector<std::vector<int>> cons_list;
    //    ////////////////////////////////////////

    SparseSet<var_t> triggered_bounds;

  // statistics
  size_t total_size{0};
  int num_units{0};

  static const bool BOOLEAN{false};
  static const bool NUMERIC{true};

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
//        Constraint<T>::idempotent = true;
    }


template <typename T>
void ClauseBase<T>::post(const int idx) {
    Constraint<T>::cons_id = idx;
//    for(var_t x{0}; x<solver.boolean.size(); ++x) {
//        solver.wake_me_on(solver.boolean.getLiteral(true, x), Constraint<T>::cons_id);
//        solver.wake_me_on(solver.boolean.getLiteral(false, x), Constraint<T>::cons_id);
//    }
    
//        for(var_t x{0}; x<solver.numeric.size(); ++x) {
//            solver.wake_me_on(lb<T>(x), Constraint<T>::cons_id);
//            solver.wake_me_on(ub<T>(x), Constraint<T>::cons_id);
//        }
    triggered_bounds.reserve(2 * solver.numeric.size());
}

// propagate the constraint
template <typename T>
void ClauseBase<T>::propagate() {
    for(auto b : triggered_bounds) {
      auto p{solver.numeric.strongestLiteral(Literal<T>::sgn(b),
                                             Literal<T>::var(b))};
//      std::cout << " --> unitprop " << solver.pretty(p) << std::endl;
      unit_propagate(p);
    }
    clearTriggers();
    //    std::cout << triggered_bounds << std::endl;
}

template <typename T> void ClauseBase<T>::clearTriggers() {
  triggered_bounds.clear();
}

template <typename T>
bool ClauseBase<T>::notify(const Literal<T> l, const int) {
//    if(not l.isNumeric())
//        unit_propagate(l);
    
    assert(l.isNumeric());
    if(not triggered_bounds.has(l)) {
        triggered_bounds.add(l);
        return true;
//        std::cout << triggered_bounds << std::endl;
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

template <typename T> void ClauseBase<T>::newBooleanVar(const var_t x) {
  watch[BOOLEAN].resize(static_cast<size_t>(2 * x + 2));
}

template <typename T> void ClauseBase<T>::newNumericVar(const var_t x) {
  watch[NUMERIC].resize(static_cast<size_t>(2 * x + 2));
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

template <typename T>
void ClauseBase<T>::unit_propagate(const Literal<T> l) {

#ifdef DBG_WATCHERS
  verifyWatchers("before UP");
#endif

#ifdef DBG_TRACE
  if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
    std::cout << "unit propagate true lit " << l << "\n";
  }
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
      
//      std::cout << " self = " << c_lit << std::endl;
//      std::cout << "other = " << other << std::endl;
      
      assert(c_lit.sameVariable(l));
      assert(c_lit.sign() != l.sign());

//      if(l.sign() == c_lit.sign()) {
//          std::cout << l << " is watched by " << *cl << " @" << idx << std::endl;
//          exit(1);
//      }
      
//      std::cout << "==> " << c_lit << std::endl;

      if (lt == NUMERIC) {
          if(c_lit.value() + l.value() >= 0) {
              
#ifdef DBG_TRACE
              if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                  std::cout << " false trigger (" << c_lit << " is not falsified)\n";
              }
#endif
              
              continue;
          } else {
              
#ifdef DBG_TRACE
              if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                  std::cout << " true trigger (" << c_lit << " and " << l << " are contradictory)\n";
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
        std::cout << "search another literal to watch" ;
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


//template <typename T>
//void ClauseBase<T>::unit_propagate_num(const Literal<T> l) {
//
//#ifdef DBG_WATCHERS
//  verifyWatchers("before UP");
//#endif
//
//#ifdef DBG_TRACE
//  if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
//    std::cout << "unit propagate true lit " << l << "\n";
//  }
//#endif
//
//
//  // watch[lt][l] contains all clauses watched by ~l.
//  // unless l is a numeric literal x <= k, in which case watch[lt][l] contains
//  // clauses watched by some literal -x <= v therefore the trigger is "real"
//  // only if k+v < 0
//  for (auto c{watch[NUMERIC][l].rbegin()}; c != watch[NUMERIC][l].rend(); ++c) {
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
////          std::cout << l << " is watched by " << *cl << " @" << idx << std::endl;
////          exit(1);
////      }
//      
////      std::cout << "==> " << c_lit << std::endl;
//
//      if (lt == NUMERIC) {
//          if(c_lit.value() + l.value() >= 0) {
//              
//#ifdef DBG_TRACE
//              if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
//                  std::cout << " false trigger (" << c_lit << " is not falsified)\n";
//              }
//#endif
//              
//              continue;
//          } else {
//              
//#ifdef DBG_TRACE
//              if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
//                  std::cout << " true trigger (" << c_lit << " and " << l << " are contradictory)\n";
//              }
//#endif
//              
//          }
//      }
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
//          swap(*c, *watch[lt][l].rbegin());
//          watch[lt][l].pop_back();
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
//#ifdef DBG_WATCHERS
//  verifyWatchers("after UP");
//#endif
//}

//template <typename T>
//void ClauseBase<T>::unit_propagate_bool(const Literal<T> l) {
//
//#ifdef DBG_WATCHERS
//  verifyWatchers("before UP");
//#endif
//
//#ifdef DBG_TRACE
//  if (DBG_CBOUND and (DBG_TRACE & UNITPROPAGATION)) {
//    std::cout << "unit propagate true lit " << l << "\n";
//  }
//#endif
//
//  // watch[lt][l] contains all clauses watched by ~l.
//  // unless l is a numeric literal x <= k, in which case watch[lt][l] contains
//  // clauses watched by some literal -x <= v therefore the trigger is "real"
//  // only if k+v < 0
//  for (auto c{watch[BOOLEAN][l].rbegin()}; c != watch[BOOLEAN][l].rend(); ++c) {
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
//          swap(*c, *watch[BOOLEAN][l].rbegin());
//          watch[BOOLEAN][l].pop_back();
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
//#ifdef DBG_WATCHERS
//  verifyWatchers("after UP");
//#endif
//}

template <typename T>
void ClauseBase<T>::set_watcher(const int r, const index_t i,
                                   Clause<T> *cl) {
    
    
//    std::cout << *cl << std::endl << "current watchers: " << cl->watched(0) << " & " << cl->watched(1) << std::endl;

  cl->watched_index[r] = i;

  Literal<T> l{~((*cl)[i])};

//      std::cout << "make " << *cl << " (" << cl->id << ") watch " << l << " (" << info_t(l) << ") because it contains " <<
//    (*cl)[i] << " (" << info_t((*cl)[i]) << ")"
//    << std::endl;
//    
//    std::cout << "new watchers: " << cl->watched(0) << " & " << cl->watched(1) << std::endl;

  watch[l.isNumeric()][l].push_back(cl);
}

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

