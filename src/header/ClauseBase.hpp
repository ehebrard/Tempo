

#ifndef _TEMPO_CLAUSEBASE_HPP
#define _TEMPO_CLAUSEBASE_HPP

#include <vector>
#include <assert.h>

#include "Clause.hpp"
#include "Failure.hpp"
#include "Global.hpp"
#include "util/Options.hpp"
#include "util/SparseSet.hpp"
#include "util/SubscribableEvent.hpp"

//#define DBG_WATCHERS

namespace tempo {

template <typename T> class Solver;
template<typename T> class Scheduler;
template<typename T> class BoundConstraint;

template<class T>
class ClauseBase : Explainer {
    
public:
  //    std::vector<Clause *> base;
  //    std::vector<Clause *> learnt;

  ClauseBase(Scheduler<T> &);
  ~ClauseBase() = default;

  void resize(const size_t n, const size_t m);
  size_t size() const;
  size_t volume() const;

  //    void unit_propagate(const var l);
  void set_watcher(const int r, const index_t i, Clause *cl);

  // "learn" makes several assumptions about the added clause, too lazy to write
  // the robust version yet
  //    template <typename iter>
  //    void add(const iter first, const iter last);

  template <typename iter>
  Clause *add(const iter first, const iter last, const bool learnt = false);

  void unit_propagate(const lit l);

  void assign(const lit l, Explanation e);

  std::ostream &display(std::ostream &os) const;

  std::ostream &displayClause(std::ostream &os, const Clause *c) const;

  bool satisfied(const genlit) const;
  bool falsified(const genlit) const;

  lit newNegLiteral(const lit l);
  lit getReasonLit(const lit l) const;

  bool sameLit(const lit, const lit) const;
  //    lit getVarInClause(const lit l) const;
  //    lit getVarInTrigger(const lit l) const;

  // explanation stuff
  void xplain(const lit l, const hint h, std::vector<lit> &Cl);
  std::ostream &print_reason(std::ostream &, const hint) const;
  int getType() const;

  void notifyRemoval(const lit l);
  void forget(Clause *cl);
  void forget_worst();
  void forget();
  void forgetAll();
  double looseness(const lit l);
  double inverseActivity(const lit l);
  double loosenessOverActivity(const lit l);

  Clause *back() {
//
//    std::cout << base.size() << " / "
//              << static_cast<size_t>(free_cl_indices.bbegin() -
//                                     free_cl_indices.fbegin())
//              << std::endl;

    return base[*(free_cl_indices.bbegin())];
  }

    std::string prettyLiteral(const genlit l) const;
    
private:
    
    Scheduler<T>& caller;

    SubscriberHandle handlerToken;

    /// Clause memory store ///
    /// The set of clauses
    std::vector<Clause *> base;
    std::vector<double> score;
    /// free indices in the vector base (so that we can remove clauses)
    SparseSet<int> free_cl_indices;
    // clauses, watch struct watch[BOUND_LIT] for bound literals,
    // watch[EDGE_LIT] for edges
    std::vector<std::vector<Clause*>> watch[2];
    ////////////////////////////////////////

    /// Literal memory store ///
    // store the set of bound constraints appearing in a literal
    std::vector<BoundConstraint<T>> constraint;
    // the number of clauses in which this bound-constraint appears
    std::vector<int> cardinality;
    // free indices in the vector constraints (so that we can remove constraints when they are not used)
    SparseSet<int> free_lit_indices;
    // for every event-lit, the set of (pointers to) constraints ordered by bound
    std::vector<std::vector<int>> cons_list;
    ////////////////////////////////////////

    // statistics
    size_t total_size{0};
    int num_units{0};

#ifdef DBG_WATCHERS
    void verifyWatchers(const char *msg) const;
#endif
};

template <class T> class NewClauseBase : NewExplainer<T> {

public:
  //    std::vector<Clause *> base;
  //    std::vector<Clause *> learnt;

  NewClauseBase(Solver<T> &);
  ~NewClauseBase() = default;

  //  void resize(const size_t n, const size_t m);
  void newBooleanVar(const var_t x);
  void newNumericVar(const var_t x);
  size_t size() const;
    size_t numLearnt() const;

  size_t volume() const;

  //    void unit_propagate(const var l);
  void set_watcher(const int r, const index_t i, NewClause<T> *cl);

  // "learn" makes several assumptions about the added clause, too lazy to write
  // the robust version yet
  //    template <typename iter>
  //    void add(const iter first, const iter last);

  template <typename iter>
  NewClause<T> *add(const iter first, const iter last,
                    const bool learnt = false);

  void unit_propagate(const Literal<T> l);

  void assign(const Literal<T> l, const NewExplanation<T> &e);

  std::ostream &display(std::ostream &os) const;

//  std::ostream &displayClause(std::ostream &os, const NewClause<T> *c) const;

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
  void forget(NewClause<T> *cl);
  void forget_worst();
  void forget();
  void forgetAll();
  double looseness(const Literal<T> l);
  double inverseActivity(const Literal<T> l);
  double loosenessOverActivity(const Literal<T> l);

  NewClause<T> *back() {
    //
    //    std::cout << base.size() << " / "
    //              << static_cast<size_t>(free_cl_indices.bbegin() -
    //                                     free_cl_indices.fbegin())
    //              << std::endl;

    return base[*(free_cl_indices.bbegin())];
  }

  //    std::string prettyLiteral(const genlit l) const;

private:
  Solver<T> &caller;

  SubscriberHandle handlerToken;

  /// Clause memory store ///
  /// The set of clauses
  std::vector<NewClause<T> *> base;
  std::vector<double> score;
  /// free indices in the vector base (so that we can remove clauses)
  SparseSet<int> free_cl_indices;
  // clauses, watch struct watch[BOUND_LIT] for bound literals,
  // watch[EDGE_LIT] for edges
  std::vector<std::vector<NewClause<T> *>> watch[2];
  ////////////////////////////////////////

  //    /// Literal memory store ///
  //    // store the set of bound constraints appearing in a literal
  //    std::vector<BoundConstraint<T>> constraint;
  //    // the number of clauses in which this bound-constraint appears
  //    std::vector<int> cardinality;
  //    // free indices in the vector constraints (so that we can remove
  //    constraints when they are not used) SparseSet<int> free_lit_indices;
  //    // for every event-lit, the set of (pointers to) constraints ordered by
  //    bound std::vector<std::vector<int>> cons_list;
  //    ////////////////////////////////////////

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

template<typename T>
std::string ClauseBase<T>::prettyLiteral(const genlit l) const {
    if(LTYPE(l) == EDGE_LIT) {
        
        //        std::cout <<
        
        assert(FROM_GEN(l) >= 0 and FROM_GEN(l) < static_cast<lit>(2*caller.numVariable()));
        
        std::stringstream ss;
        ss << caller.getEdge(FROM_GEN(l));
        return ss.str();
    } else {
        
        
        //        std::cout << FROM_GEN(l) << " / " << constraint.size() << std::endl;
        
        assert(FROM_GEN(l) >= 0 and FROM_GEN(l) < static_cast<lit>(constraint.size()));
        
        
        
        std::stringstream ss;
        ss << constraint[FROM_GEN(l)];
        return ss.str();
    }
}

template <typename T>
ClauseBase<T>::ClauseBase(Scheduler<T> &c)
    : caller(c), handlerToken(caller.SearchRestarted.subscribe_handled(
                     [this]() { this->forget(); })) {}



template<typename T>
size_t ClauseBase<T>::size() const {

  if ((base.size() - free_cl_indices.size()) != free_cl_indices.backsize()) {
    std::cout << "what ?\n";
    exit(1);
  }

  return free_cl_indices.backsize() +
         free_cl_indices.frontsize(); // base.size() - free_cl_indices.size();
                                      // // + learnt.size();
  //    return var_level.size();
}

template<typename T>
size_t ClauseBase<T>::volume() const {
    return total_size; // + learnt.size();
    //    return var_level.size();
}

template <typename T>
void ClauseBase<T>::resize(const size_t n, const size_t m) {
  watch[EDGE_LIT].resize(2 * m);
  watch[BOUND_LIT].resize(2 * n);
  cons_list.resize(2 * n);
}

template<typename T>
bool ClauseBase<T>::satisfied(const genlit gl) const {
    lit l{FROM_GEN(gl)};
    auto lt{LTYPE(gl)};
    if(lt == EDGE_LIT)
        return caller.satisfied(l);
    return caller.satisfied(constraint[l]);
}

template<typename T>
bool ClauseBase<T>::falsified(const genlit gl) const {
    lit l{FROM_GEN(gl)};
    auto lt{LTYPE(gl)};
    if(lt == EDGE_LIT)
        return caller.falsified(l);
    return caller.falsified(constraint[l]);
}

template<typename T>
void ClauseBase<T>::assign(const lit l, Explanation e) {
    
    if(LTYPE(l) == EDGE_LIT) {
        
        //        std::cout << " + true egde literal " << prettyLiteral(l) << std::endl;
        
        caller.set(FROM_GEN(l), e);
    } else {
        
//                std::cout << " + true bound literal " << prettyLiteral(l) << std::endl;
        
        caller.set(constraint[FROM_GEN(l)], e);
    }
}

//template<typename T>
//index_t ClauseBase<T>::getWatchRank(const lit gl, const Clause* cl) {
//    if(LTYPE(gl) == EDGE_LIT)
//        return cl->watch_rank(FROM_GEN(gl));
//    lit l{cl->watcher(1)};
//    if(LTYPE(l) == EDGE_LIT)
//        return 0;
//    auto q{FROM_GEN(l)};
//    return constraint[q].l == caller.getBoundLiteral(q);
//}

template<typename T>
void ClauseBase<T>::unit_propagate(const lit gl) {

#ifdef DBG_WATCHERS
  verifyWatchers("before UP");
#endif

#ifdef DBG_TRACE
    if (DBG_TRACE & UNITPROPAGATION) {
        std::cout << "unit propagate false lit " << caller.prettyLiteral(gl) << "\n";
    }
#endif
    
    lit l{FROM_GEN(gl)};
    auto lt{LTYPE(gl)};
    
    
    
    //    if(lt == EDGE_LIT)
    //        for(auto p : binary[l])
    //        {
    //            assign(p,{this,l});
    //        }
    ////    else {
    ////
    ////    }
    
    //    auto wdix{(lt == EDGE_LIT ? NOT(l) : //NOT(constraint[c.l].l))};
    
    for (auto c{watch[lt][l].rbegin()}; c != watch[lt][l].rend(); ++c) {
        
        auto cl{*c};
        
#ifdef DBG_TRACE
        if (DBG_TRACE & UNITPROPAGATION) {
            std::cout << " watched by ";
            displayClause(std::cout, cl);
            std::cout << std::endl;
        }
#endif
        
        bool watch_rank{0}; // = cl->watch_rank(gl);
        index_t idx; //= cl->watcher_index[watch_rank];
        lit other; // = cl->watcher(1 - watch_rank);
        lit p;
        
        if(lt == BOUND_LIT) {
            BoundConstraint<T> lit_cons; //{constraint[other]};
            
            // the bound literal that triggered UP
            BoundConstraint<T> negated_cons{~(caller.getBound(l))};
            
            other = cl->watcher(1); // type is not bound lit, then the trigger corresponds to watch[0]
            
            if(LTYPE(other) == BOUND_LIT) {
                lit_cons = constraint[FROM_GEN(other)];
                if(lit_cons.l == negated_cons.l)
                {
                    // else check if it matches, and if so, swap
                    other = cl->watcher(0);
                    watch_rank = 1;
                } else {
                    // otherwise get the right constraint
                    lit_cons = constraint[cl->watcher(0)];
                }
            }
            
            // check that the literal in the clause is actually falsified
            if(not negated_cons.entails(lit_cons))
            {
#ifdef DBG_TRACE
                if (DBG_TRACE & UNITPROPAGATION) {
                    std::cout << " false trigger (" << lit_cons << " is not falsified)";
                }
#endif
                
                continue;
            }
        }  else {
            
            watch_rank = cl->watch_rank(gl);
            other = cl->watcher(1 - watch_rank);
            
        }
        
        idx = cl->watcher_index[watch_rank];
        
        if (satisfied(other)) {
            
#ifdef DBG_TRACE
            if (DBG_TRACE & UNITPROPAGATION) {
                std::cout << ": satisfied by " << prettyLiteral(other) << std::endl;
            }
#endif
            
            continue;
        }
        
        index_t i{idx};
        
        while (true) {
            
            if (++i == cl->size())
                i = 0;
            if (i == idx)
                break;
            
            // look for a replacement
            p = (*cl)[i];
            
#ifdef DBG_TRACE
            if (DBG_TRACE & UNITPROPAGATION) {
                std::cout << "  " << prettyLiteral(p);
            }
#endif
            
            if (p != other) {
                if (not falsified(p)) {
                    
#ifdef DBG_TRACE
                    if (DBG_TRACE & UNITPROPAGATION) {
                        std::cout << ": replace by " << prettyLiteral(p) << " (" << i << ") as "
                        << watch_rank << "-th watcher " ;
                    }
#endif
                    
                    set_watcher(watch_rank, i, cl);
                    
#ifdef DBG_TRACE
                    if (DBG_TRACE & UNITPROPAGATION) {
                        std::cout << " and rm from "  <<
                        prettyLiteral(gl) << "'s watches\n";
                    }
#endif
                    
                    // remove clause from l's watch list
                    swap(*c, *watch[lt][l].rbegin());
                    watch[lt][l].pop_back();
                    
                    break;
                }
            }
#ifdef DBG_TRACE
            else if (DBG_TRACE & UNITPROPAGATION) {
                std::cout << "*";
            }
#endif
        }
        
        if (i == idx) {
            
#ifdef DBG_TRACE
            if (DBG_TRACE & UNITPROPAGATION) {
                std::cout << ": new unit " << prettyLiteral(other) << std::endl;
            }
#endif
            
            // there was no replacement
            // unit_literals.push_back(other);
            // model[VAR(other)] = SIGN(other);
            //      assign(other, cl);
            //            assign(other,{cl,static_cast<hint>(cl->watch_index(NOT(watch_rank)))});
            
            
//            std::cout << "unit prop " << prettyLiteral(other) << " b/c cl[" << cl->id << "]: " ;
//            displayClause(std::cout, cl);
//            std::cout << "\n";

#ifdef DBG_WATCHERS
            verifyWatchers("at assign");
#endif
            assign(other,{this,cl->id});
            
            assert(cl == base[cl->id]);
            
            // for(var x{0}; x<num_variable(); ++x) {
            //   cout << x << ": " << model[x] << endl;
            // }
        }
#ifdef DBG_TRACE
        else if (DBG_TRACE & UNITPROPAGATION) {
            std::cout << " new watchers: ";
            displayClause(std::cout, cl);
            std::cout << std::endl;
        }
#endif
        
        // cout << "remaining: " << distance(c, watch[EDGE_LIT][l].rend()) << endl;
    }

#ifdef DBG_WATCHERS
    verifyWatchers("after UP");
#endif
}

template<typename T>
void ClauseBase<T>::set_watcher(const int r, const index_t i, Clause *cl) {
  //#ifdef DBG_WATCHERS
  //    verifyWatchers("before set watcher");
  //#endif

  cl->watcher_index[r] = i;

  //    std::cout << "watcher[" << r << "] = " << cl->watcher(r) << std::endl;
  // add it to p's watch list
  lit l{(*cl)[i]};

  if (LTYPE(l) == EDGE_LIT)
    watch[EDGE_LIT][FROM_GEN(l)].push_back(cl);
  else
    watch[BOUND_LIT][constraint[FROM_GEN(l)].l].push_back(cl);

  //#ifdef DBG_WATCHERS
  //    verifyWatchers("after set watcher");
  //#endif
}

//lit getVarInClause(const lit l) const {
//    if(LTYPE(l) == EDGE) {
//        return FROM_GEN(l);
//    }
//    return constraint[FROM_GEN(l)]
//}
//
//lit getVarInTrigger(const lit l) const;

template<typename T>
bool ClauseBase<T>::sameLit(const lit l_in_clause, const lit trigger_l) const {
    if(LTYPE(l_in_clause) != LTYPE(trigger_l))
        return false;
    if(LTYPE(l_in_clause) == EDGE_LIT)
        return l_in_clause == trigger_l;
    return constraint[FROM_GEN(l_in_clause)].l == caller.getBoundLiteral(FROM_GEN(trigger_l));
}


template<typename T>
lit ClauseBase<T>::getReasonLit(const lit l) const {
    if(LTYPE(l) == EDGE_LIT)
        return EDGE(NOT(FROM_GEN(l)));
    auto c{~(constraint[FROM_GEN(l)])};
    auto bl{caller.getImplicant(c)};
    return BOUND(bl);
}


template<typename T>
lit ClauseBase<T>::newNegLiteral(const lit l) {

#ifdef DBG_WATCHERS
  verifyWatchers("before new lit");
#endif

//      std::cout << l << std::endl;
  //    std::cout << "create a clause literal negation for " <<
  //    caller.prettyLiteral(l) << std::endl;
    
//    if(l == 15122) {
//        std::cout << *this << std::endl;
//    }

  if (LTYPE(l) == EDGE_LIT) {
    lit q{EDGE(NOT(FROM_GEN(l)))};
    //        std::cout << " => " << caller.prettyLiteral(q) << std::endl;
    return q;
  }

    
//    std::cout << "beg new lit " << (free_lit_indices.capacity() > 64 and free_lit_indices.has(64)) << std::endl;
    
  auto c{~(caller.getBound(FROM_GEN(l)))};

  //    std::cout << " => bound " << c << std::endl;

  lit li{NoLit};

  auto it{cons_list[c.l].begin()};
  for (; it != cons_list[c.l].end(); ++it) {
    if (constraint[*it].distance > c.distance)
      continue;
    if (constraint[*it].distance == c.distance)
        li = *it;
    break;
  }

  if (li != NoLit) {
      
//    std::cout << " already used (" << li << ")" << std::endl;

    ++cardinality[li];
      
//    std::cout << "reuse bound lit " << li << ": " << prettyLiteral(BOUND(li)) << " (" << cardinality[li] << ")" << std::endl;
  } else {
    if (not free_lit_indices.empty()) {

#ifdef DBG_WATCHERS
      verifyWatchers("here");
#endif
        
//        if(l == 15122) {
//            for (auto cl : watch[BOUND_LIT][84]) {
//                    displayClause(std::cout, cl);
//                std::cout << std::endl;
//            }
//        }

      //            std::cout << free_lit_indices << std::endl;

      li = free_lit_indices.back();
      //            free_lit_indices.remove_front(li);
      //            free_lit_indices.pop_back();
      free_lit_indices.remove_back(li);

      //            std::cout << free_lit_indices << std::endl;

      //            std::cout << constraint[li] << std::endl;
      //
      //            for(auto acl : watch[BOUND_LIT][constraint[li].l]) {
      //                std::cout << " watched by ";
      //                displayClause(std::cout, acl);
      //                std::cout << std::endl;
      //            }
      //
      //            assert(watch[BOUND_LIT][constraint[li].l].empty());
      //
      constraint[li] = c;
      cardinality[li] = 1;

//      std::cout << " recycle index (" << li << ") for " << c << std::endl;
        
        
//        if(l == 15122) {
//            for (auto cl : watch[BOUND_LIT][84]) {
//                    displayClause(std::cout, cl);
//                std::cout << std::endl;
//            }
//        }
        
//        std::cout << "new bound lit " << li << ": " << prettyLiteral(BOUND(li)) << " (recycled)" << std::endl;

#ifdef DBG_WATCHERS
      verifyWatchers("there");
#endif

    } else {
      li = static_cast<lit>(constraint.size());
      constraint.push_back(c);
      cardinality.push_back(1);
      free_lit_indices.reserve(constraint.size());
        
//        std::cout << "new bound lit " << li << ": " << prettyLiteral(BOUND(li)) << " (brand new)" << std::endl;

//      std::cout << " add new index (" << li << ")" << std::endl;
    }

    cons_list[c.l].insert(it, li);
  }

  //    std::cout << " return " << BOUND(li) << std::endl;

    
    
#ifdef DBG_WATCHERS
  verifyWatchers("after new lit");
#endif

  return BOUND(li);
}

template <typename T> void ClauseBase<T>::notifyRemoval(const lit l) {

  //    std::cout << l << " / " << cardinality.size() << std::endl;
  //
  //    std::cout << cardinality[l] << std::endl;

//if(l==64)
//
  if (--cardinality[l] == 0) {
      
//      std::cout << "remove bound lit " << l << ": " << prettyLiteral(BOUND(l)) << std::endl;
      
      if(free_lit_indices.has(l)) {
          std::cout << "ERROR\n";
          exit(1);
      }
    free_lit_indices.add(l);
      
      
      auto c{constraint[l]};
//      lit li{NoLit};

      auto it{cons_list[c.l].begin()};
      for (; it != cons_list[c.l].end(); ++it) {
        if (constraint[*it].distance > c.distance)
          continue;
        assert(l == *it);
          break;
      }
      
      cons_list[c.l].erase(it);

      
  }
//  else {
//      std::cout << "-- bound lit " << l << ": " << prettyLiteral(BOUND(l)) << " (" << cardinality[l] << ")" << std::endl;
//  }
}

template <typename T> void ClauseBase<T>::forget(Clause *cl) {

//      std::cout << "forget ";
//      displayClause(std::cout, cl);
//      std::cout << std::endl;
//    assert(base[cl->id] == cl);
    
//    std::cout << "beg forget " << (free_lit_indices.capacity() > 64 and free_lit_indices.has(64)) << std::endl;
//    if(base.size() > 43 and not base[43]->empty())
//    {
//        displayClause(std::cout, base[43]);
//        std::cout << std::endl;
//    }

  for (auto r{0}; r < 2; ++r) {
    auto w{cl->watcher(r)};
    auto wt{LTYPE(w)};
    auto wl{(wt == EDGE_LIT ? FROM_GEN(w) : constraint[FROM_GEN(w)].l)};

    //        std::cout << "remove from list of " << prettyLiteral(w) << " (" <<
    //        (wt == EDGE_LIT ? "edge" : "bound") << " / " << (wt == EDGE_LIT ?
    //        caller.prettyLiteral(EDGE(wl)) : prettyEventLit(wl)) << ")\n";

    for (auto c{watch[wt][wl].begin()}; c != watch[wt][wl].end(); ++c)
      if (*c == cl) {
        swap(*c, watch[wt][wl].back());
        watch[wt][wl].pop_back();
        break;
      }
  }

  for (auto l : *cl) {
    if (LTYPE(l) == BOUND_LIT) {
      notifyRemoval(FROM_GEN(l));
    }
  }
  free_cl_indices.add(cl->id);
    
  total_size -= cl->size();

  cl->clear();
    
//    if(base.size() > 43 and not base[43]->empty())
//    {
//        displayClause(std::cout, base[43]);
//        std::cout << std::endl;
//    }
//    std::cout << "end forget " << (free_lit_indices.capacity() > 64 and free_lit_indices.has(64)) << std::endl;
}

/// | |  1 2  6 5 3 4

template <typename T> void ClauseBase<T>::forget_worst() {
  forget(base[*(free_cl_indices.bbegin())]);

  //    std::cout << "forget ";
  //    displayClause(std::cout, cl);
  //    std::cout << std::endl;
  //
  //    for(auto l : *cl) {
  //        if(LTYPE(l) == BOUND_LIT) {
  //            notifyRemoval(FROM_GEN(l));
  //        }
  //    }
  //    free_cl_indices.push_back();
  //    total_size -= cl->size();
  //    cl->clear();
}

template <typename T> double ClauseBase<T>::loosenessOverActivity(const lit l) {
  if (LTYPE(l) == EDGE_LIT) {
    auto c{caller.getEdge(FROM_GEN(l))};
    return caller.looseness(c) / caller.getActivityMap()->get(c);
  } else {
    auto c{constraint[FROM_GEN(l)]};
    return caller.looseness(c) / caller.getActivityMap()->get(c);
  }
}

template <typename T> double ClauseBase<T>::inverseActivity(const lit l) {
  if (LTYPE(l) == EDGE_LIT) {
    auto c{caller.getEdge(FROM_GEN(l))};
    return 1.0 / caller.getActivityMap()->get(c);
  } else {
    auto c{constraint[FROM_GEN(l)]};
    return 1.0 / caller.getActivityMap()->get(c);
  }
}

template <typename T> double ClauseBase<T>::looseness(const lit l) {
  if (LTYPE(l) == EDGE_LIT) {
    auto c{caller.getEdge(FROM_GEN(l))};
    return caller.looseness(c);
  } else {
    auto c{constraint[FROM_GEN(l)]};
    return caller.looseness(c);
  }
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

  if (caller.getOptions().forget_strategy == Options::LiteralScore::Looseness or
      (caller.getActivityMap() == NULL and
       caller.getOptions().forget_strategy != Options::LiteralScore::Size)) {
    for (auto idx{free_cl_indices.bbegin()}; idx != free_cl_indices.bend();
         ++idx) {
      score[*idx] = 0;
      for (auto l : *base[*idx]) {
        score[*idx] += looseness(l);
      }
    }
  } else if (caller.getOptions().forget_strategy ==
             Options::LiteralScore::LoosenessOverActivity) {
    for (auto idx{free_cl_indices.bbegin()}; idx != free_cl_indices.bend();
         ++idx) {
      score[*idx] = 0;
      for (auto l : *base[*idx]) {
        score[*idx] += loosenessOverActivity(l);
      }
    }
  } else if (caller.getOptions().forget_strategy ==
             Options::LiteralScore::Activity) {
    for (auto idx{free_cl_indices.bbegin()}; idx != free_cl_indices.bend();
         ++idx) {
      score[*idx] = 0;
      for (auto l : *base[*idx]) {
        score[*idx] += inverseActivity(l);
      }
    }
  }


  if (caller.getOptions().forget_strategy == Options::LiteralScore::Size) {
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

  auto target_size = static_cast<size_t>(
      static_cast<double>(size()) * (1.0 - caller.getOptions().forgetfulness));

  while (size() > target_size) {
    forget_worst();
  }


    
    
#ifdef DBG_WATCHERS
  verifyWatchers("after forget");
#endif
}

template <typename T>
template <typename iter>
Clause *ClauseBase<T>::add(const iter first, const iter last,
                           const bool learnt) {
    
#ifdef DBG_WATCHERS
  verifyWatchers("before add");
#endif
  Clause *c{NULL};

  assert(first != last);

  if (first + 1 == last) {
    assign(*first, Constant::NoReason);
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
      c = new Clause(id);
      base.push_back(c);
      score.push_back(0);
      free_cl_indices.reserve(base.size());
      if (not learnt) {
        free_cl_indices.add(id);
        free_cl_indices.remove_front(id);
      }
    }

    for (auto l{first}; l != last; ++l) {
      c->push_back(*l);
    }

    set_watcher(0, 0, c);
    set_watcher(1, 1, c);

    lit l{(*c)[0]};
      
      total_size += c->size();

      
//      auto sbefore{caller.numBoundLiteral() + caller.numEdgeLiteral()};
//      std::cout << " literal ";
//      if(LTYPE(l) == EDGE_LIT) {
//                  std::cout << prettyLiteral(l) << std::endl;
//      } else {
//                  std::cout << prettyLiteral(l) << std::endl;
//      }
    assign(l, {this, c->id});
//      
//      auto safter{caller.numBoundLiteral() + caller.numEdgeLiteral()};
//      assert(sbefore < safter);
//      
//      std::cout << "\n\n\n" << safter - sbefore << "\n\n\n" << std::endl;
  }

#ifdef DBG_WATCHERS
  verifyWatchers("after add");
#endif

  return c;
}

template<typename T>
std::ostream& ClauseBase<T>::displayClause(std::ostream &os, const Clause* cl) const {
    //    os << "[" << cl->watch_index(0) << "|" << cl->watch_index(1) << "]";
    
    //    assert(not cl->empty());
    
//    if (cl->empty()) {
//        os << "(" << cl->id <<":)";
//    } else {
//        
//        os << "(";
//        os << cl->id <<":";
////        os << "[" << prettyLiteral((*cl)[0]) << "]";
////        if (0 == cl->watch_index(0) or 0 == cl->watch_index(1))
////            os << "*";
////        for (size_t i{1}; i < cl->size(); ++i) {
////            os << " or [" << prettyLiteral((*cl)[i]) << "]";
////            if (i == cl->watch_index(0) or i == cl->watch_index(1))
////                os << "*";
////        }
//        if(LTYPE((*cl)[0]) == BOUND_LIT)
//            os << FROM_GEN((*cl)[0]);
//        else
//            os << ".";
//        for (size_t i{1}; i < cl->size(); ++i) {
//            if(LTYPE((*cl)[i]) == BOUND_LIT)
//                os << " " << FROM_GEN((*cl)[i]);
//            else
//                os << " .";
//        }
//        os << ")";
//    }
    
    if (cl->empty()) {
        os << "()";
    } else {
        os << "(";
        os << "[" << prettyLiteral((*cl)[0]) << "]";
        if (0 == cl->watch_index(0) or 0 == cl->watch_index(1))
            os << "*";
        for (size_t i{1}; i < cl->size(); ++i) {
            os << " or [" << prettyLiteral((*cl)[i]) << "]";
            if (i == cl->watch_index(0) or i == cl->watch_index(1))
                os << "*";
        }
        os << ")";
    }
    return os;
}

template<typename T>
void ClauseBase<T>::xplain(const lit l, const hint h, std::vector<lit> &Cl) {

    Clause& reason(*(base[h]));
    
    if(l == NoLit) {
        for(auto p : reason)
            Cl.push_back(getReasonLit(p));
    } else {

        for(auto p : reason)
            if(not sameLit(p,l)) {
                Cl.push_back(getReasonLit(p));
            }
    }
}

template<typename T>
std::ostream &ClauseBase<T>::print_reason(std::ostream &os, const hint h) const {
    return displayClause(os, base[h]);
}

template<typename T>
int ClauseBase<T>::getType() const {
    return CLAUSEEXPL;
}

template <typename T>
std::ostream &ClauseBase<T>::display(
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

  for (size_t x{0}; x < caller.numVariable(); ++x) {
    if (not watch[EDGE_LIT][POS(x)].empty()) {
      os << caller.prettyLiteral(EDGE(POS(x))) << " is watched in";
      for (auto cl : watch[EDGE_LIT][POS(x)]) {
        //                os << *cl;
        os << " ";
        displayClause(os, cl);
      }
      os << std::endl;
    }
    if (not watch[EDGE_LIT][NEG(x)].empty()) {
      os << caller.prettyLiteral(EDGE(NEG(x))) << " is watched in";
      for (auto cl : watch[EDGE_LIT][NEG(x)]) {
        //                os << *cl;
        os << " ";
        displayClause(os, cl);
      }
      os << std::endl;
    }
  }
  for (size_t x{0}; x < caller.numEvent(); ++x) {
    if (not watch[BOUND_LIT][LOWERBOUND(x)].empty()) {
      os << prettyEventLit(LOWERBOUND(x)) << " is watched in";
      for (auto cl : watch[BOUND_LIT][LOWERBOUND(x)]) {
        //                os << *cl;
        os << " ";
        displayClause(os, cl);
      }
      os << std::endl;
    }
    if (not watch[BOUND_LIT][UPPERBOUND(x)].empty()) {
      os << prettyEventLit(UPPERBOUND(x)) << " is watched in";
      for (auto cl : watch[BOUND_LIT][UPPERBOUND(x)]) {
        //                os << *cl;
        os << " ";
        displayClause(os, cl);
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
    for(auto cl : base) {
        if(cl->id != i)
        {
            std::cout << msg << "indexing error" << std::endl;
            exit(1);
        }
        if(not free_cl_indices.has(cl->id)) {
            for(auto l : *cl) {
                if(LTYPE(l) == BOUND_LIT) {
                    if(free_lit_indices.has(FROM_GEN(l))) {
                        std::cout << msg << ": error on free index (" << FROM_GEN(l) << ") for " << prettyLiteral(l) << " in clause " << cl->id << std::endl;
                        displayClause(std::cout, cl);
                        displayClause(std::cout, base[cl->id]);
                        std::cout << std::endl;
                        exit(1);
                    }
                }
            }
        }
        ++i;
    }

  size_t num_watchers{0};
  for (size_t x{0}; x < caller.numVariable(); ++x) {
    num_watchers += watch[EDGE_LIT][POS(x)].size();
    if (not watch[EDGE_LIT][POS(x)].empty()) {
      for (auto cl : watch[EDGE_LIT][POS(x)]) {

        if (cl->size() < 2) {
          std::cout << msg << ": error empty clause watching "
                    << caller.prettyLiteral(EDGE(POS(x))) << std::endl;
          exit(1);
        }

        if (free_cl_indices.has(cl->id)) {
          displayClause(std::cout, cl);
          std::cout << "'s id is free " << std::endl;
          exit(1);
        }

        if (cl->watcher(0) != EDGE(POS(x)) and cl->watcher(1) != EDGE(POS(x))) {
          std::cout << msg << ": error on clause " << cl->id << " -- ";
          displayClause(std::cout, cl);
          std::cout << " watching " << caller.prettyLiteral(EDGE(POS(x)))
                    << std::endl;
          exit(1);
        }
      }
    }
    if (not watch[EDGE_LIT][NEG(x)].empty()) {
      num_watchers += watch[EDGE_LIT][NEG(x)].size();
      for (auto cl : watch[EDGE_LIT][NEG(x)]) {
        if (cl->size() < 2) {
          std::cout << msg << ": error empty clause watching "
                    << caller.prettyLiteral(EDGE(NEG(x))) << std::endl;
          exit(1);
        }

        if (free_cl_indices.has(cl->id)) {
          displayClause(std::cout, cl);
          std::cout << "'s id is free " << std::endl;
          exit(1);
        }

        if (cl->watcher(0) != EDGE(NEG(x)) and cl->watcher(1) != EDGE(NEG(x))) {
          std::cout << msg << ": error on clause " << cl->id << " -- ";
          displayClause(std::cout, cl);
          std::cout << " watching " << caller.prettyLiteral(EDGE(NEG(x)))
                    << std::endl;
          exit(1);
        }
      }
    }
  }
  for (size_t x{0}; x < caller.numEvent(); ++x) {
    num_watchers += watch[BOUND_LIT][LOWERBOUND(x)].size();
    if (not watch[BOUND_LIT][LOWERBOUND(x)].empty()) {
      for (auto cl : watch[BOUND_LIT][LOWERBOUND(x)]) {
        if (cl->size() < 2) {
          std::cout << msg << ": error empty clause watching "
                    << prettyEventLit(LOWERBOUND(x)) << std::endl;
          exit(1);
        }

        if (free_cl_indices.has(cl->id)) {
          displayClause(std::cout, cl);
          std::cout << "'s id is free " << std::endl;
          exit(1);
        }

        auto w0{(LTYPE(cl->watcher(0)) == BOUND_LIT
                     ? constraint[FROM_GEN(cl->watcher(0))].l
                     : NoLit)};
        auto w1{(LTYPE(cl->watcher(1)) == BOUND_LIT
                     ? constraint[FROM_GEN(cl->watcher(1))].l
                     : NoLit)};

        //                auto w0{constraint[FROM_GEN(cl->watcher(0))].l};
        //                auto w1{constraint[FROM_GEN(cl->watcher(1))].l};
        if (w0 != LOWERBOUND(x) and w1 != LOWERBOUND(x)) {
          std::cout << msg << ": error on clause " << cl->id << " -- ";
          displayClause(std::cout, cl);
          std::cout << " watching " << prettyEventLit(LOWERBOUND(x))
                    << std::endl;
          std::cout << (LTYPE(cl->watcher(0)) == BOUND_LIT
                            ? FROM_GEN(cl->watcher(0))
                            : NoLit)
                    << " / "
                    << (LTYPE(cl->watcher(1)) == BOUND_LIT
                            ? FROM_GEN(cl->watcher(1))
                            : NoLit)
                    << " / " << LOWERBOUND(x) << std::endl;
          exit(1);
        }
      }
    }
    if (not watch[BOUND_LIT][UPPERBOUND(x)].empty()) {
      num_watchers += watch[BOUND_LIT][UPPERBOUND(x)].size();
      for (auto cl : watch[BOUND_LIT][UPPERBOUND(x)]) {
        if (cl->size() < 2) {
          std::cout << msg << ": error empty clause watching "
                    << prettyEventLit(LOWERBOUND(x)) << std::endl;
          exit(1);
        }

        if (free_cl_indices.has(cl->id)) {
          displayClause(std::cout, cl);
          std::cout << "'s id is free " << std::endl;
          exit(1);
        }

        //                displayClause(std::cout, cl);
        //                std::cout << " watching " <<
        //                prettyEventLit(UPPERBOUND(x)) << std::endl;

        //                if(free_lit_inidices.has(FROM_GEN(cl->watcher(0)))) {
        //                    std::cout << prettyLiteral(cl->watcher(0)) << "'s
        //                    id is free " << std::endl;
        //                }

        auto w0{(LTYPE(cl->watcher(0)) == BOUND_LIT
                     ? constraint[FROM_GEN(cl->watcher(0))].l
                     : NoLit)};
        auto w1{(LTYPE(cl->watcher(1)) == BOUND_LIT
                     ? constraint[FROM_GEN(cl->watcher(1))].l
                     : NoLit)};

        if (w0 != UPPERBOUND(x) and w1 != UPPERBOUND(x)) {
          std::cout << msg << ": error on clause " << cl->id << " -- ";
          displayClause(std::cout, cl);
          std::cout << " watching " << prettyEventLit(UPPERBOUND(x))
                    << std::endl;
          std::cout << (LTYPE(cl->watcher(0)) == BOUND_LIT
                            ? FROM_GEN(cl->watcher(0))
                            : NoLit)
                    << " / "
                    << (LTYPE(cl->watcher(1)) == BOUND_LIT
                            ? FROM_GEN(cl->watcher(1))
                            : NoLit)
                    << " / " << UPPERBOUND(x) << std::endl;
          exit(1);
        }
      }
    }
  }

    if(num_watchers != 2*size()) {
        std::cout << msg << ": wrong number of watchers !\n";
        exit(1);
    }
  //    std::cout << msg << ": ok\n";
}
#endif

template<typename T>
std::ostream &operator<<(std::ostream &os, const tempo::ClauseBase<T> &x) {
  return x.display(os);
}


////// NEW CLAUSES
///

template <typename T>
NewClauseBase<T>::NewClauseBase(Solver<T> &c)
    : caller(c), handlerToken(caller.SearchRestarted.subscribe_handled(
                     [this]() { this->forget(); })) {}


template<typename T>
size_t NewClauseBase<T>::numLearnt() const {
    return free_cl_indices.backsize();
}

template <typename T> size_t NewClauseBase<T>::size() const {

  if ((base.size() - free_cl_indices.size()) != (free_cl_indices.backsize() + free_cl_indices.frontsize())) {
    std::cout << "what ?\n";
    exit(1);
  }

  return free_cl_indices.backsize() + free_cl_indices.frontsize();
}

template <typename T> size_t NewClauseBase<T>::volume() const {
  return total_size;
}

// template <typename T>
// void NewClauseBase<T>::resize(const size_t n, const size_t m) {
//   watch[BOOLEAN].resize(2 * n);
//   watch[NUMERIC].resize(2 * m);
// }

template <typename T> void NewClauseBase<T>::newBooleanVar(const var_t x) {
  watch[BOOLEAN].resize(static_cast<size_t>(2 * x + 2));
}

template <typename T> void NewClauseBase<T>::newNumericVar(const var_t x) {
  watch[NUMERIC].resize(static_cast<size_t>(2 * x + 2));
}

template <typename T>
bool NewClauseBase<T>::satisfied(const Literal<T> l) const {
  if (l.isNumeric())
    return caller.numeric.satisfied(l);
  else
    return caller.boolean.satisfied(l);
}

template <typename T>
bool NewClauseBase<T>::falsified(const Literal<T> l) const {
  if (l.isNumeric())
    return caller.numeric.falsified(l);
  else
    return caller.boolean.falsified(l);
}

// template <typename T>
// void NewClauseBase<T>::assign(const Literal<T> l, const NewExplanation<T> &e)
// {
//   caller.set(l, e);
// }

template <typename T>
void NewClauseBase<T>::unit_propagate(const Literal<T> l) {

#ifdef DBG_WATCHERS
  verifyWatchers("before UP");
#endif

#ifdef DBG_TRACE
  if (DBG_TRACE & UNITPROPAGATION) {
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
    if (DBG_TRACE & UNITPROPAGATION) {
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
              if (DBG_TRACE & UNITPROPAGATION) {
                  std::cout << " false trigger (" << c_lit << " is not falsified)\n";
              }
#endif
              
              continue;
          } else {
              
#ifdef DBG_TRACE
              if (DBG_TRACE & UNITPROPAGATION) {
                  std::cout << " true trigger (" << c_lit << " and " << l << " are contradictory)\n";
              }
#endif
              
          }
      }

    if (satisfied(other)) {

#ifdef DBG_TRACE
      if (DBG_TRACE & UNITPROPAGATION) {
        std::cout << ": satisfied by " << other << std::endl;
      }
#endif

      continue;
    }

    index_t i{idx};

#ifdef DBG_TRACE
      if (DBG_TRACE & UNITPROPAGATION) {
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
      if (DBG_TRACE & UNITPROPAGATION) {
        std::cout << "  " << p;
      }
#endif

      if (p != other) {
        if (not falsified(p)) {

#ifdef DBG_TRACE
          if (DBG_TRACE & UNITPROPAGATION) {
            std::cout << ": replace by " << p << " (" << i << ") as "
                      << watch_rank << "-th watcher ";
          }
#endif

          set_watcher(watch_rank, i, cl);

#ifdef DBG_TRACE
          if (DBG_TRACE & UNITPROPAGATION) {
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
      else if (DBG_TRACE & UNITPROPAGATION) {
        std::cout << "*";
      }
#endif
    }

    if (i == idx) {

#ifdef DBG_TRACE
      if (DBG_TRACE & UNITPROPAGATION) {
        std::cout << ": new unit " << other << std::endl;
      }
#endif

#ifdef DBG_WATCHERS
      verifyWatchers("at assign");
#endif

      caller.set(other, {this, cl->id});

      assert(cl == base[cl->id]);
    }
  }

#ifdef DBG_WATCHERS
  verifyWatchers("after UP");
#endif
}

template <typename T>
void NewClauseBase<T>::set_watcher(const int r, const index_t i,
                                   NewClause<T> *cl) {
    
    
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

template <typename T> void NewClauseBase<T>::forget(NewClause<T> *cl) {
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

template <typename T> void NewClauseBase<T>::forget_worst() {
    
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
double NewClauseBase<T>::loosenessOverActivity(const Literal<T> l) {
  return caller.looseness(l) / caller.getActivityMap()->get(l, caller);
}

template <typename T>
double NewClauseBase<T>::inverseActivity(const Literal<T> l) {
  return 1.0 / caller.getActivityMap()->get(l, caller);
}

template <typename T> double NewClauseBase<T>::looseness(const Literal<T> l) {
  return caller.looseness(l);
}

template <typename T> void NewClauseBase<T>::forgetAll() {
  while (free_cl_indices.backsize() > 0) {
    forget_worst();
  }
}

template <typename T> void NewClauseBase<T>::forget() {

  if (size() < 1000)
    return;

#ifdef DBG_WATCHERS
  verifyWatchers("before forget");
#endif

  if (caller.getOptions().forget_strategy == Options::LiteralScore::Looseness or
      ( // caller.getActivityMap() == NULL and
          caller.getOptions().forget_strategy != Options::LiteralScore::Size)) {
    for (auto idx{free_cl_indices.bbegin()}; idx != free_cl_indices.bend();
         ++idx) {
      score[*idx] = 0;
      for (auto l : *base[*idx]) {
        score[*idx] += looseness(l);
      }
    }
  } else if (caller.getOptions().forget_strategy ==
             Options::LiteralScore::LoosenessOverActivity) {
    for (auto idx{free_cl_indices.bbegin()}; idx != free_cl_indices.bend();
         ++idx) {
      score[*idx] = 0;
      for (auto l : *base[*idx]) {
        score[*idx] += loosenessOverActivity(l);
      }
    }
  } else if (caller.getOptions().forget_strategy ==
             Options::LiteralScore::Activity) {
    for (auto idx{free_cl_indices.bbegin()}; idx != free_cl_indices.bend();
         ++idx) {
      score[*idx] = 0;
      for (auto l : *base[*idx]) {
        score[*idx] += inverseActivity(l);
      }
    }
  }

  if (caller.getOptions().forget_strategy == Options::LiteralScore::Size) {
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

  auto target_size = static_cast<size_t>(
      static_cast<double>(numLearnt()) * (1.0 - caller.getOptions().forgetfulness));

  while (numLearnt() > target_size) {
    forget_worst();
  }

#ifdef DBG_WATCHERS
  verifyWatchers("after forget");
#endif
}

template <typename T>
template <typename iter>
NewClause<T> *NewClauseBase<T>::add(const iter first, const iter last,
                                    const bool learnt) {

#ifdef DBG_WATCHERS
  verifyWatchers("before add");
#endif

  //  std::cout << "add clause (" << (learnt ? "learnt" : "base" ) << ")";
  //  for (auto l{first}; l != last; ++l) {
  //    std::cout << " " << *l;
  //  }
  //  std::cout << std::endl;

  NewClause<T> *c{NULL};

  assert(first != last);

  if (first + 1 == last) {
    //    assign(*first, Constant::NewNoReason<T>);
    caller.set(*first, Constant::NewNoReason<T>);
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
      c = new NewClause<T>(id);
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
          caller.set(l, {this, c->id});
      }
  }

#ifdef DBG_WATCHERS
  verifyWatchers("after add");
#endif

  return c;
}

//template <typename T>
//std::ostream &NewClauseBase<T>::displayClause(std::ostream &os,
//                                              const NewClause<T> *cl) const {
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
void NewClauseBase<T>::xplain(const Literal<T> l, const hint h,
                              std::vector<Literal<T>> &Cl) {

  //    int i{0};
  //    for(auto cl : base) {
  //        assert(cl->id == i);
  //        ++i;
  //    }

  NewClause<T> &reason(*(base[h]));

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
      if (not l.sameVariable(p)) {
        Cl.push_back(~p);
      }
  }
}

template <typename T>
std::ostream &NewClauseBase<T>::print_reason(std::ostream &os,
                                             const hint h) const {
//  return displayClause(os, base[h]);
    os << *(base[h]);
    return os;
}

template <typename T> int NewClauseBase<T>::getType() const {
  return CLAUSEEXPL;
}

template <typename T>
std::ostream &NewClauseBase<T>::display(
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

  for (size_t x{0}; x < caller.boolean.size(); ++x) {
    auto l{caller.boolean.getLiteral(true, x)};
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
  for (size_t x{0}; x < caller.numeric.size(); ++x) {
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
void NewClauseBase<T>::verifyWatchers(const char *msg) const {
    
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
    for (size_t x{0}; x < caller.boolean.size(); ++x) {
        Literal<T> l{caller.boolean.getLiteral(true,x)};
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
    for (size_t x{0}; x < caller.numeric.size(); ++x) {
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
std::ostream &operator<<(std::ostream &os, const NewClauseBase<T> &x) {
  return x.display(os);
}

}

#endif // _TEMPO_CLAUSEBASE_HPP

