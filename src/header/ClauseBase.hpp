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
#include "heuristics/impl/ActivityMap.hpp"

//#define DBG_WATCHERS

#define NEW_WATCHERS

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
    ~ClauseBase();
    
    // notify that the solver has Boolean variable x
    void newBooleanVar();
    
    // notify that the solver has numeric variable x
    void newNumericVar();
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
    Clause<T> *operator[](const index_t i) const;
    // returns the last clause to have been added
    Clause<T> *back();
    
    // create and add a clause from a list of literals (and returns it) 'learnt'
    // flag distinguished constraints from cuts
    template <typename iter>
    Clause<T> *add(const iter first, const iter last, const bool learnt = false, const int glue_score=1);
    //    template <typename Iterable>
    Clause<T> *add(const std::vector<Literal<T>>& lits) {
        return add(lits.begin(), lits.end());
    }
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
    
    // return the list of clauses that watch literal l

    
#ifdef NEW_WATCHERS
    std::list<std::pair<T, index_t>>::iterator
    getNumWatchList(const Literal<T> l);
    index_t getWatchList(const Literal<T> l);
    index_t getOrCreateWatchList(const Literal<T> l);
#else
    std::vector<Clause<T> *> &getWatchList(const Literal<T> l);
    std::vector<Clause<T> *> &getOrCreateWatchList(const Literal<T> l);
#endif
    // set the 'r'-th watcher of clause 'cl' to be its 'i'-th literal
    void set_watcher(const int r, const index_t i, Clause<T> *cl);
    
    //  // set the 'r'-th watcher of clause 'cl' to be its 'i'-th literal (and
    //  order
    //  // the watch-list of this literal)
    //  void set_watcher_numeric(const int r, const index_t i, Clause<T> *cl);
    
    // works for both numeric and Boolean, does not order numeric watch-lists
    void unit_propagate(const Literal<T> l);
    
#ifdef NEW_WATCHERS
    // helper to avoid code duplication
    void up(const Literal<T> l, const index_t widx);//std::vector<Clause<T> *> &watches);
    
    // works only for Boolean literals
    void unit_propagate_boolean(const Literal<T> l);
    
    // works only for numeric literals
    void unit_propagate_numeric(const Literal<T> l);
#else
    // works only for Boolean literals
    void unit_propagate_boolean(const Literal<T> l);
    
    // works only for numeric literals, orders numeric watch-lists
    void unit_propagate_numeric(const Literal<T> l);
#endif
    
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
    std::ostream &displayClauses(std::ostream &os) const;
    std::ostream &displayWatchStruct(std::ostream &os) const;
    //@}
    
    /**
     * @name explanation
     */
    //@{
    void xplain(const Literal<T> l, const hint h, std::vector<Literal<T>> &Cl);
    std::ostream &print_reason(std::ostream &, const hint) const;
    //  int getType() const;
    //@}
    
    /**
     * @name clause forgetting
     */
    //@{
    // make a clause unforgettable
    void makeUnforgettable(Clause<T> *cl);
    // forgets a learnt clause
    void forget(Clause<T> *cl);
    // forgets the last learnt clause (worst if they have been sorted)
    void forget_worst();
    // forgets according to the forgetting policy in solver.getOptions()
    void forget();
    // forgets all learnt clauses
    void forgetAll();
    // remove all clauses
    void clear();
    // literal activity score
    double activity(const Literal<T> l);
    double learningRate(const Literal<T> l);
    // literal score based on its semantic
    double looseness(const Literal<T> l);
    // literal score based on its activity
    double inverseActivity(const Literal<T> l);
    double inverseLearningRate(const Literal<T> l);
    // literal score based on both its semantic and its activity
    double loosenessOverActivity(const Literal<T> l);
    double loosenessOverLearningRate(const Literal<T> l);
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
    std::vector<int> glue;
    // free indices in the vector base (so that we can remove clauses)
    // - contain the indices in 'base' that are available (no current clause has
    // this id)
    // - at the back of the sparse set are the indices used by learnt clauses
    // - at the front of the sparse set are the indices used by other clauses
    SparseSet<int> free_cl_indices;
    
   
    
#ifdef NEW_WATCHERS
    // a watch-list for every Boolean literal
    std::vector<std::vector<Clause<T> *>> watchers;
    SparseSet<> free_wl_indices;
    
    // for each signed numeric variable, a list of pairs <v,l> with v a numeric
    // bound and l a watch-list (ordered by increasing v)
    std::vector<std::list<std::pair<T, index_t>>> numeric_watch;
    std::vector<index_t> boolean_watch;
    
//    // the actual watch lists, not sorted
//    std::vector<std::vector<Clause<T> *>> unsorted_numeric_watch;

    
    //    //TODO: tmp debug thing
    //    std::vector<Clause<T> *>* watch_list;
#else
    // the watch lists for every literel (watch[BOOLEAN] and watch[NUMERIC])
    std::vector<std::vector<Clause<T> *>> watch[2];
#endif
    
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
                                                                   [this](const bool) { this->forget(); })) {
    
    Constraint<T>::priority = Priority::Low;
    
#ifdef NEW_WATCHERS
    watchers.resize(1);
    free_wl_indices.reserve(1);
#endif
    
}

template<class T>
ClauseBase<T>::~ClauseBase() {
    for (auto *c: base) {
        delete c;
    }
}

template <typename T>
void ClauseBase<T>::post(const int idx) {
    
    Constraint<T>::cons_id = idx;
    if (solver.getOptions().full_up) {
        //    for(var_t x{0}; x<solver.boolean.size(); ++x) {
        //        solver.wake_me_on(solver.boolean.getLiteral(true, x),
        //        Constraint<T>::cons_id);
        //        solver.wake_me_on(solver.boolean.getLiteral(false, x),
        //        Constraint<T>::cons_id);
        //    }
        for (var_t x{static_cast<var_t>(solver.getNumericScope(idx).size())}; x < solver.numeric.size(); ++x) {
            
//            std::cout << "wake me on " << x << std::endl;
            
            solver.wake_me_on(lb<T>(x), Constraint<T>::cons_id);
            solver.wake_me_on(ub<T>(x), Constraint<T>::cons_id);
        }
        triggered_bounds.reserve(2 * solver.numeric.size());
    }
//    
//    
//    
//    std::cout << " #vars=" << solver.getNumericScope(idx).size();
}

// propagate the constraint
template <typename T>
void ClauseBase<T>::propagate() {

    ++num_prop;
    
    for (auto b : triggered_bounds) {
      ++solver.num_unit_propagations;
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
    
//    if ((base.size() - free_cl_indices.size()) !=
//        (free_cl_indices.backsize() + free_cl_indices.frontsize())) {
//        std::cout << "what ?\n";
//        exit(1);
//    }
    
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

template <typename T> Clause<T> *ClauseBase<T>::operator[](const index_t i) const {
    if (i >= free_cl_indices.capacity() or free_cl_indices.has(i))
        return NULL;
    return base[i];
}

template <typename T> Clause<T> *ClauseBase<T>::back() {
    return base[*(free_cl_indices.bbegin())];
}

template <typename T> void ClauseBase<T>::newBooleanVar() {
#ifdef NEW_WATCHERS
//    auto p{static_cast<size_t>(Literal<T>::index(true, x))};
//    auto n{static_cast<size_t>(Literal<T>::index(false, x))};
//    auto m{std::max(p,n)};
//    boolean_watch.resize(m+1);
//    boolean_watch[n] = watchers.size();
//    boolean_watch[p] = watchers.size()+1;
//    watchers.resize(watchers.size()+2);
    boolean_watch.resize(2*solver.boolean.size());
    *boolean_watch.rbegin() = watchers.size()+1;
    *(boolean_watch.rbegin()+1) = watchers.size();
    watchers.resize(watchers.size()+2);
#else
    //  watch[BOOLEAN].resize(static_cast<size_t>(2 * x + 2));
//    watch[BOOLEAN].resize(static_cast<size_t>(Literal<T>::index(true, x) + 1));
    watch[BOOLEAN].resize(2*solver.boolean.size());
#endif
}

template <typename T> void ClauseBase<T>::newNumericVar() {
#ifdef NEW_WATCHERS
    //    std::vector<Clause<T>*> empty_watchlist;
    //    std::pair<T,std::vector<Clause<T>*>>
    //    empty_watchlists(Constant::Infinity<T>,empty_watchlist);
    
    
//    if(static_cast<size_t>(Literal<T>::index(true, x) + 1) != 2*solver.numeric.size()) {
//        std::cout << static_cast<size_t>(Literal<T>::index(true, x) + 1) << " != " << 2*solver.numeric.size() << std::endl;
//        exit(1);
//    }
    
    
    auto sz{numeric_watch.size()};
    numeric_watch.resize(2*solver.numeric.size());
    while (sz < numeric_watch.size()) {
        numeric_watch[sz].insert(numeric_watch[sz].end(),
                                 {Constant::Infinity<T>, 0});
        //         static_cast<index_t>(unsorted_numeric_watch.size())});
        //    unsorted_numeric_watch.resize(unsorted_numeric_watch.size() + 1);
        //    //            numeric_watch[sz].insert(numeric_watch[sz].end(),
        //    //            empty_watchlists);
        ++sz;
    }
#else
    //    watch[NUMERIC].resize(static_cast<size_t>(2 * x + 2));
    watch[NUMERIC].resize(2*solver.numeric.size());
#endif
    //    numeric_watch.resize(static_cast<size_t>(Literal<T>::index(x,true)+1));
    //    count.resize(static_cast<size_t>(2 * x + 2), 0);
    
    
//    if (solver.getOptions().full_up) {
//        for (var_t x{static_cast<var_t>(solver.getNumericScope(Constraint<T>::cons_id).size())}; x < solver.numeric.size(); ++x) {
//            
//            std::cout << "*wake me on " << x << std::endl;
//            
//            solver.wake_me_on(lb<T>(x), Constraint<T>::cons_id);
//            solver.wake_me_on(ub<T>(x), Constraint<T>::cons_id);
//        }
//        triggered_bounds.reserve(2 * solver.numeric.size());
//    }
    
    
    
//    std::cout << " #vars=" << solver.getNumericScope(idx).size();
    
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
//   if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
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

#ifdef NEW_WATCHERS

template <typename T>
void ClauseBase<T>::up(const Literal<T> l, const index_t widx){ //std::vector<Clause<T> *> &watches) {
    
#ifdef DBG_WATCHERS
    verifyWatchers("before up");
#endif

    //    if(free_cl_indices.capacity() > 3124 and (not
    //    free_cl_indices.has(3124))) {
    //        std::cout << base[3124]->size() << " / " <<
    //        solver.num_choicepoints << std::endl;
    //    }

    //  std::cout << "watches.size() = " << watches.size() << std::endl;
    
    //  for (auto c{watch[lt][l].rbegin()}; c != watch[lt][l].rend(); ++c) {
//    for (auto c{watchers[widx].rbegin()}; c != watchers[widx].rend(); ++c) {
    size_t k{watchers[widx].size()};
//    size_t i{n};
    while(k-->0) {
        
//        std::cout << k << "/" << n << std::endl;
        
//        auto cl{*c};
        auto cl{watchers[widx][k]};
        
#ifdef DBG_TRACE
        if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
            std::cout << " info_t = " << info_t(l) << " /watched by " << *cl
            << std::endl;
        }
#endif

        assert(cl->watch_index(1) < cl->size());
        //        if(cl->watch_index(1) >= cl->size()) {
        //            std::cout << "this bug!\n";
        //            exit(1);
        //        }
        //        if (cl->watch_index(1) >= cl->size()) {
        //          std::cout << "u*cl_" << cl->id << ": " << cl->watch_index(1)
        //                    << " >= " << cl->size()
        //                    << " #cp=" << solver.num_choicepoints << " widx="
        //                    << widx
        //                    << " k=" << k << std::endl;
        //          exit(1);
        //        }

        bool watch_rank{cl->watch_rank(l)};
        index_t idx{cl->watched_index[watch_rank]};
        Literal<T> other{cl->watched(1 - watch_rank)};
        //      Literal<T> c_lit{(*cl)[idx]};
        
        //      assert(c_lit.sameVariable(l));
        //      assert(c_lit.sign() != l.sign());
        
        if (satisfied(other)) {
            
#ifdef DBG_TRACE
            if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                std::cout << ": satisfied by " << other << std::endl;
            }
#endif
            
            continue;
        }
        
        index_t i{idx};
        
#ifdef DBG_TRACE
        if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
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
            if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                std::cout << "  " << p;
            }
#endif
            
            if (p != other) {
                if (not falsified(p)) {
                    
#ifdef DBG_TRACE
                    if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                        std::cout << ": replace by " << p << " (" << i << ") as "
                        << watch_rank << "-th watcher ";
                    }
#endif
                    
                    //          std::cout << "watches.size() = " << watches.size() <<
                    //          std::endl;
                    //
                    //            displayWatchStruct(std::cout);
                    
                    set_watcher(watch_rank, i, cl);
                    
#ifdef DBG_TRACE
                    if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                        std::cout << " and rm from " << l << "'s watches\n";
                    }
#endif
                    
                    //            displayWatchStruct(std::cout);
                    
                    // remove clause from l's watch list

                    assert(watchers[widx].size() >= 0);

                    //          std::cout << "watches.size() = " << watches.size() <<
                    //          std::endl; std::cout << "rm " << **c << " from \n"; for
                    //          (auto cla : watches) {
                    //            std::cout << *cla << std::endl;
                    //          }
                    //          std::cout << ".\n";
                    
//                    swap(*c, *watches.rbegin());
                    
//                    std::cout << "remove " << *cl << " from:\n";
//                    for(auto ocl : watchers[widx]) {
//                        std::cout << *ocl << std::endl;
//                    }
                    
//                    typename std::vector<Clause<T> *>::iterator kl{watchers[widx].begin()};
//                    kl += i;
                    
//                    std::cout << i << "/" << watchers[widx].size() << std::endl;
                    
                    assert((*(watchers[widx].begin()+k))->id == cl->id);
                    
                    
                    swap(*(watchers[widx].begin()+k), *watchers[widx].rbegin());
                    watchers[widx].pop_back();
                    
                    break;
                }
            }
#ifdef DBG_TRACE
            else if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                std::cout << "*";
            }
#endif
        }
        
        if (i == idx) {
            
#ifdef DBG_TRACE
            if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                std::cout << "-> new unit @lvl " << solver.level() << ": " << other
                << std::endl;
            }
#endif
            
#ifdef DBG_WATCHERS
            verifyWatchers("at assign");
#endif
            
#ifdef DBGP0
            if (solver.level() <= solver.init_level) {
                std::cout << "pruning @lvl " << solver.init_level << std::endl;
            }
#endif
            
            if (solver.level() <= solver.init_level) {
                makeUnforgettable(cl);
            }
            
            solver.set(other, {this, cl->id});
            
            assert(cl == base[cl->id]);
        }
    }
}

template <typename T>
void ClauseBase<T>::unit_propagate_boolean(const Literal<T> l) {
    
    // boolean_watch[l] contains all clauses watched by ~l.
    // unless l is a numeric literal x <= k, in which case watch[lt][l] contains
    // clauses watched by some literal -x <= v therefore the trigger is "real"
    // only if k+v < 0
    
//    std::cout << l << " (" << static_cast<int>(l) << ") / " << boolean_watch.size() << std::endl;
    
    up(l, boolean_watch[l]);
}

template <typename T>
void ClauseBase<T>::unit_propagate_numeric(const Literal<T> l) {
    
    //    std::cout << l << "'s watch-lists:";
    //    for(auto wl{getNumWatchList(l)}; wl!=numeric_watch[l].end(); ++wl) {
    //        std::cout << "<"<< wl->first << ":" << wl->second << ">";
    //    }
    //    std::cout << ".\n";
    
    auto wl{getNumWatchList(l)};
    while (wl->first < Constant::Infinity<T>) {
//        auto &watches{unsorted_numeric_watch[wl->second]};
        auto next{wl};
        ++next;
        
//        if (watches.empty()) {
//            std::cout << "\nunit-prop <" << wl->first << ":" << wl->second << "> "
//            << l << std::endl;
//            displayWatchStruct(std::cout);
//            
//            std::cout << "\nverify watchers\n";
//            
//            verifyWatchers("wtf");
//            
//            exit(1);
//        }
        
//        assert(not watches.empty());
        
        //    std::cout << "->\n";
        //    for (auto cla : watches) {
        //      std::cout << *cla << std::endl;
        //    }
        //    std::cout << ".\n";
        
        up(l, wl->second);
        if (watchers[wl->second].empty()) {
            
            //        std::cout << "erase <" << wl->first << ":" << unsorted_numeric<<
            //        ">" <<
            
            //        std::cout << "add " << wl->second << " to " << free_wl_indices
            //        << std::endl;
            
            free_wl_indices.add(wl->second);

            //            if(numeric_watch[l].empty()) {
            //                std::cout << "bug!\n";
            //                exit(1);
            //            }

            //            auto flag{numeric_watch[l].size() == 15};
            //            if(flag) {
            //                std::cout << solver.num_unit_propagations << ":
            //                erase <" << wl->first << ", " << wl->second << "
            //                (" << watchers[wl->second].size() << ")> from:\n";
            //                for(auto x : numeric_watch[l]) {
            //                    std::cout << "<" << x.first << ", " <<
            //                    x.second << " (" << watchers[x.second].size()
            //                    << ")>\n";
            //                }
            //            }
            ////            std::cout << numeric_watch[l].size() << std::endl;

            numeric_watch[l].erase(wl);

            //            std::cout << "->" << numeric_watch[l].size() <<
            //            std::endl;
        }
        
        wl = next;
    }
    
    //    for(auto wl{getNumWatchList(l)}; wl->first < Constant::Infinity<T>;
    //    ++wl) {
    //        auto& watches{unsorted_numeric_watch[wl->second]};
    //
    //        up(l, watches);
    //        if(watches.empty()) {
    //
    //            std::cout << "erase\n";
    //
    //            numeric_watch[l].erase(wl);
    //        }
    //
    //    }
}
#endif

template <typename T> void ClauseBase<T>::unit_propagate(const Literal<T> l) {
    
    ++num_up;
    
#ifdef DBG_WATCHERS
    verifyWatchers("before UP");
#endif
    
#ifdef DBG_TRACE
    if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
        std::cout << "\nunit propagate " << l << std::endl;
    }
#endif
    
#ifdef NEW_WATCHERS
    if (l.isNumeric()) {
        unit_propagate_numeric(l);
    } else {
        unit_propagate_boolean(l);
    }
#else
    
    auto lt{litType(l)};
    
    // watch[lt][l] contains all clauses watched by ~l.
    // unless l is a numeric literal x <= k, in which case watch[lt][l] contains
    // clauses watched by some literal -x <= v therefore the trigger is "real"
    // only if k+v < 0
    auto &watches{getWatchList(l)};
    
    //  for (auto c{watch[lt][l].rbegin()}; c != watch[lt][l].rend(); ++c) {
    for (auto c{watches.rbegin()}; c != watches.rend(); ++c) {
        
        auto cl{*c};
        
#ifdef DBG_TRACE
        if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
            std::cout << " watched by " << *cl << std::endl;
        }
#endif

        if (cl->watch_index(1) >= cl->size()) {
          std::cout << "*cl_" << cl->id << ": " << cl->watch_index(1)
                    << " >= " << cl->size() << std::endl;
          exit(1);
        }

        bool watch_rank{cl->watch_rank(l)};
        index_t idx{cl->watched_index[watch_rank]};
        Literal<T> other{cl->watched(1 - watch_rank)};
        Literal<T> c_lit{(*cl)[idx]};
        
        assert(c_lit.sameVariable(l));
        assert(c_lit.sign() != l.sign());
        
        if (lt == NUMERIC) {
            if (c_lit.value() + l.value() >= 0) {
                
#ifdef DBG_TRACE
                if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                    std::cout << " false trigger (" << c_lit << " is not falsified)\n";
                }
#endif
                
                ++num_miss;
                continue;
            } else {
                
#ifdef DBG_TRACE
                if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                    std::cout << " true trigger (" << c_lit << " and " << l
                    << " are contradictory)\n";
                }
#endif
            }
        }
        
        if (satisfied(other)) {
            
#ifdef DBG_TRACE
            if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                std::cout << ": satisfied by " << other << std::endl;
            }
#endif
            
            continue;
        }
        
        index_t i{idx};
        
#ifdef DBG_TRACE
        if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
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
            if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                std::cout << "  " << p;
            }
#endif
            
            if (p != other) {
                if (not falsified(p)) {
                    
#ifdef DBG_TRACE
                    if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                        std::cout << ": replace by " << p << " (" << i << ") as "
                        << watch_rank << "-th watcher ";
                    }
#endif
                    
                    set_watcher(watch_rank, i, cl);
                    
#ifdef DBG_TRACE
                    if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                        std::cout << " and rm from " << l << "'s watches\n";
                    }
#endif
                    
                    // remove clause from l's watch list
                    swap(*c, *watches.rbegin());
                    watches.pop_back();
                    
                    break;
                }
            }
#ifdef DBG_TRACE
            else if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                std::cout << "*";
            }
#endif
        }
        
        if (i == idx) {
            
#ifdef DBG_TRACE
            if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                std::cout << "-> new unit @lvl " << solver.level() << ": " << other
                << std::endl;
            }
#endif
            
#ifdef DBG_WATCHERS
            verifyWatchers("at assign");
#endif
            
#ifdef DBGP0
            if (solver.level() <= solver.init_level) {
                std::cout << "pruning @lvl " << solver.init_level << std::endl;
            }
#endif
            
            if (solver.level() <= solver.init_level) {
                makeUnforgettable(cl);
            }
            
            solver.set(other, {this, cl->id});
            
            assert(cl == base[cl->id]);
        }
    }
#endif
    
#ifdef DBG_WATCHERS
    verifyWatchers("after UP");
#endif
}

#ifndef NEW_WATCHERS
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
        
        // #ifdef DBG_TRACE
        //     if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
        //       for (index_t i{0}; i < watch[NUMERIC][l].size(); ++i) {
        //         std::cout << std::setw(3) << i << " "
        //                   << watch[NUMERIC][l][i]->watched(
        //                          watch[NUMERIC][l][i]->watch_rank(l))
        //                   << " " << *(watch[NUMERIC][l][i]) << std::endl;
        //       }
        //       std::cout << "stack:";
        //       for (auto i : search_stack) {
        //         std::cout << " " << i;
        //       }
        //       std::cout << std::endl;
        //     }
        // #endif
        
        auto k{search_stack.back()};
        
        if (static_cast<size_t>(k) >= watch[NUMERIC][l].size()) {
            search_stack.pop_back();
            continue;
        }
        
        auto cl{watch[NUMERIC][l][k]};

        if (cl->watch_index(1) >= cl->size()) {
          std::cout << solver.pretty(l) << " n*cl_" << cl->id << ": "
                    << cl->watch_index(1) << " >= " << cl->size() << " : "
                    << solver.num_choicepoints << std::endl;
          exit(1);
        }

        bool watch_rank{cl->watch_rank(l)};
        index_t idx{cl->watched_index[watch_rank]};
        Literal<T> other{cl->watched(1 - watch_rank)};
        Literal<T> c_lit{(*cl)[idx]};
        
#ifdef DBG_TRACE
        if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
            std::cout << "explore clause " << *cl << " (" << c_lit << ")\n";
        }
#endif
        
        if (c_lit.value() + l.value() >= 0) {
            
#ifdef DBG_TRACE
            if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
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
            if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
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
        if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
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
            if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                std::cout << "  " << p;
            }
#endif
            
            if (p != other) {
                if (not falsified(p)) {
                    
#ifdef DBG_TRACE
                    if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                        std::cout << ": replace by " << p << " (" << i << ") as "
                        << watch_rank << "-th watcher ";
                    }
#endif
                    
                    //            if(p.isNumeric())
                    //                set_watcher_numeric(watch_rank, i, cl);
                    //            else
                    set_watcher(watch_rank, i, cl);
                    
#ifdef DBG_TRACE
                    if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
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
            else if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                std::cout << "*";
            }
#endif
        }
        
        if (i == idx) {
            
#ifdef DBG_TRACE
            if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                std::cout << "-> new unit @lvl " << solver.level() << ": " << other
                << std::endl;
                //        std::cout << ": new unit " << other << std::endl;
            }
#endif
            
#ifdef DBG_WATCHERS
            verifyWatchers("at assign");
#endif
            
#ifdef DBGP0
            if (solver.level() <= solver.init_level) {
                std::cout << "pruning @lvl " << solver.init_level << std::endl;
            }
#endif
            if (solver.level() <= solver.init_level) {
                makeUnforgettable(cl);
            }
            
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
//    if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
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
//              if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
//                  std::cout << " false trigger (" << c_lit << " is not
//                  falsified)\n";
//              }
//#endif
//
//              continue;
//          } else {
//
//#ifdef DBG_TRACE
//              if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
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
//      if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
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
//      if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
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
//      if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
//        std::cout << "  " << p;
//      }
//#endif
//
//      if (p != other) {
//        if (not falsified(p)) {
//
//#ifdef DBG_TRACE
//          if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
//            std::cout << ": replace by " << p << " (" << i << ") as "
//                      << watch_rank << "-th watcher ";
//          }
//#endif
//
//          set_watcher(watch_rank, i, cl);
//
//#ifdef DBG_TRACE
//          if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
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
//      else if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
//        std::cout << "*";
//      }
//#endif
//    }
//
//    if (i == idx) {
//
//#ifdef DBG_TRACE
//      if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
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
        if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
            std::cout << " watched by " << *cl << std::endl;
        }
#endif

        if (cl->watch_index(1) >= cl->size()) {
          std::cout << "b*cl_" << cl->id << ": " << cl->watch_index(1)
                    << " >= " << cl->size() << std::endl;
          exit(1);
        }

        bool watch_rank{cl->watch_rank(l)};
        index_t idx{cl->watched_index[watch_rank]};
        Literal<T> other{cl->watched(1 - watch_rank)};
        //    Literal<T> c_lit{(*cl)[idx]};
        
        if (satisfied(other)) {
            
#ifdef DBG_TRACE
            if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                std::cout << ": satisfied by " << other << std::endl;
            }
#endif
            
            continue;
        }
        
        index_t i{idx};
        
#ifdef DBG_TRACE
        if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
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
            if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                std::cout << "  " << p;
            }
#endif
            
            if (p != other) {
                if (not falsified(p)) {
                    
#ifdef DBG_TRACE
                    if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                        std::cout << ": replace by " << p << " (" << i << ") as "
                        << watch_rank << "-th watcher ";
                    }
#endif
                    
                    set_watcher(watch_rank, i, cl);
                    
#ifdef DBG_TRACE
                    if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
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
            else if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                std::cout << "*";
            }
#endif
        }
        
        if (i == idx) {
            
#ifdef DBG_TRACE
            if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
                std::cout << "-> new unit @lvl " << solver.level() << ": " << other
                << std::endl;
                //        std::cout << ": new unit " << other << std::endl;
            }
#endif
            
#ifdef DBG_WATCHERS
            verifyWatchers("at assign");
#endif
            
#ifdef DBGP0
            if (solver.level() <= solver.init_level) {
                std::cout << "pruning @lvl " << solver.init_level << std::endl;
            }
#endif
            if (solver.level() <= solver.init_level) {
                makeUnforgettable(cl);
            }
            
            solver.set(other, {this, cl->id});
            
            assert(cl == base[cl->id]);
        }
    }
}
#endif

#ifdef NEW_WATCHERS
template <typename T>
std::list<std::pair<T, index_t>>::iterator
ClauseBase<T>::getNumWatchList(const Literal<T> l) {
    auto &wl(numeric_watch[l]);
    auto li{wl.begin()}; // list of watch-lists, ordered by
    
    //    std::cout << li->first << " / " << l << std::endl;
    while (
           /*li != watch_lists.end() and*/ li->first < l.value()) {
               // no need to check for the end because there is a sentinel with
               // <Constant::Infinity<T>,empty>
               ++li;
           }
    //    std::cout << "-> " << li->first << " / " << l << std::endl;
    
//    unsorted_numeric_watch.reserve(unsorted_numeric_watch.size() +
//                                   unsorted_numeric_watch[li->second].size());
    
    return li;
}
#endif



#ifdef NEW_WATCHERS
// return the watch-list of literal l, or the iterator after its expected
// position if l is numeric and has no watch-list yet
template <typename T>
index_t ClauseBase<T>::getWatchList(const Literal<T> l) {
    if (l.isNumeric()) {
        return getNumWatchList(l)->second;
    } else {
        return boolean_watch[l];
    }
}


// return the watch-list of literal l, or the iterator after its expected
// position if l is numeric and has no watch-list yet
template <typename T>
index_t
ClauseBase<T>::getOrCreateWatchList(const Literal<T> l) {
    if (l.isNumeric()) {
        assert(l.value() < Constant::Infinity<T>);
        auto wl = getNumWatchList(l);
        if (wl->first == l.value()) {
            return wl->second;
        } else {
            auto idx{static_cast<index_t>(watchers.size())};
            if (not free_wl_indices.empty()) {
                idx = free_wl_indices.front();
                free_wl_indices.remove_back(idx);
            } else {
                watchers.resize(watchers.size() + 1);
                free_wl_indices.reserve(watchers.size());
            }
            numeric_watch[l].insert(wl, {l.value(), idx});
            
            return idx;
        }
    } else {
        return boolean_watch[l];
    }
}
#else
// return the watch-list of literal l, or the iterator after its expected
// position if l is numeric and has no watch-list yet
template <typename T>
std::vector<Clause<T> *> &ClauseBase<T>::getWatchList(const Literal<T> l) {
    return watch[litType(l)][l];
}
#endif

template <typename T>
void ClauseBase<T>::set_watcher(const int r, const index_t i, Clause<T> *cl) {
    
    // #ifdef DBG_WATCHERS
    //             verifyWatchers("before set watcher");
    // #endif
    
    cl->watched_index[r] = i;
    Literal<T> l{~((*cl)[i])};
    
#ifdef NEW_WATCHERS
    auto &watches{watchers[getOrCreateWatchList(l)]};
#else
    auto &watches{getWatchList(l)};
#endif
    
    if (l.isNumeric()) {
        watches.push_back(cl);
        
        if (solver.getOptions().order_bound_watch) {
            heap::percolate_up(watches.begin(), watches.size() - 1,
                               [l](const Clause<T> *c1, const Clause<T> *c2) {
                return c1->watched(c1->watch_rank(l)).value() <
                c2->watched(c2->watch_rank(l)).value();
            });
        }
    } else {
        watches.push_back(cl);
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

template <typename T> void ClauseBase<T>::makeUnforgettable(Clause<T> *cl) {
    //  if (not free_cl_indices.isback(cl->id)) {
    //    std::cout << "clause " << *cl << " is not learnt\n";
    //
    //    if (free_cl_indices.isfront(cl->id)) {
    //      std::cout << "(marked)\n";
    //    } else {
    //      std::cout << "(unreferenced)\n";
    //    }
    //
    //    exit(1);
    //  }
    
    if (free_cl_indices.isback(cl->id)) {
        free_cl_indices.add(cl->id);
        free_cl_indices.remove_front(cl->id);
    }
}

template <typename T> void ClauseBase<T>::forget(Clause<T> *cl) {

  //    if(cl->id == 3124) {
  //        std::cout << "forget cl3124 @" << solver.num_choicepoints << "\n";
  //
  //        std::cout << " watching " << cl->watched(0) << " and " <<
  //        cl->watched(1) << std::endl;
  //
  //
  //    }

  for (auto r{0}; r < 2; ++r) {
    auto wl{~(cl->watched(r))};

#ifdef NEW_WATCHERS
        if(wl.isNumeric()) {
            auto watchlist{getNumWatchList(wl)};
            auto &watches{watchers[watchlist->second]};
            for (auto c{watches.begin()}; c != watches.end(); ++c)
                if (*c == cl) {
                    swap(*c, watches.back());
                    watches.pop_back();
                    break;
                }
            if(watches.empty()) {
                free_wl_indices.add(watchlist->second);
                numeric_watch[wl].erase(watchlist);
            }
        } else {
            auto &watches{watchers[getWatchList(wl)]};
            for (auto c{watches.begin()}; c != watches.end(); ++c)
                if (*c == cl) {
                    swap(*c, watches.back());
                    watches.pop_back();
                    break;
                }
        }
#else
        auto &watches{getWatchList(wl)};
        for (auto c{watches.begin()}; c != watches.end(); ++c)
            if (*c == cl) {
                swap(*c, watches.back());
                watches.pop_back();
                break;
            }
#endif
  }

    free_cl_indices.add(cl->id);
    
    total_size -= cl->size();
    
    cl->clear();

    //
    //    if(cl->id == 3124 and solver.num_choicepoints == 8576) {
    //        exit(1);
    //    }
}

/// | |  1 2  6 5 3 4

template <typename T> void ClauseBase<T>::forget_worst() {
    forget(base[*(free_cl_indices.bbegin())]);
}

template <typename T> double ClauseBase<T>::activity(const Literal<T> l) {
  return solver.getActivity(l);
}

template <typename T> double ClauseBase<T>::learningRate(const Literal<T> l) {
  return solver.getLearningRate(l);
}

template <typename T>
double ClauseBase<T>::loosenessOverActivity(const Literal<T> l) {
    return solver.looseness(l) * inverseActivity(l);
}

template <typename T>
double ClauseBase<T>::loosenessOverLearningRate(const Literal<T> l) {
  auto ilr{inverseLearningRate(l)};
  return (ilr == Constant::Infinity<double> ? ilr : solver.looseness(l) * ilr);
}

template <typename T>
double ClauseBase<T>::inverseActivity(const Literal<T> l) {
    return heuristics::impl::ActivityMap::baseIncrement / activity(l); // solver.getActivityMap()->get(l, solver);
}

template <typename T>
double ClauseBase<T>::inverseLearningRate(const Literal<T> l) {
  auto lr{learningRate(l)};
  return (lr != 0 ? 1 / lr : Constant::Infinity<double>);
}

template <typename T> double ClauseBase<T>::looseness(const Literal<T> l) {
    return solver.looseness(l);
}

template <typename T> void ClauseBase<T>::forgetAll() {
    
//    std::cout << "clauses: " << free_cl_indices << std::endl;
    
    while (free_cl_indices.backsize() > 0) {
        forget_worst();
//        std::cout << "clauses: " << free_cl_indices << std::endl;
    }
}

template <typename T> void ClauseBase<T>::clear() {
    while (free_cl_indices.backsize() > 0) {
        forget_worst();
    }
    while (free_cl_indices.frontsize() > 0) {
        forget(base[*(free_cl_indices.frbegin())]);
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

    if (solver.getOptions().forget_strategy ==
        Options::LiteralScore::GlueTimeActivity) {
      for (auto idx{free_cl_indices.bbegin()}; idx != free_cl_indices.bend();
           ++idx) {
        score[*idx] = 0;
        for (auto l : *base[*idx]) {
          score[*idx] += inverseActivity(l);
        }
        score[*idx] *= glue[*idx];
      }
    } else if (solver.getOptions().forget_strategy ==
               Options::LiteralScore::GlueTimeLearningRate) {
      for (auto idx{free_cl_indices.bbegin()}; idx != free_cl_indices.bend();
           ++idx) {
        score[*idx] = 0;
        for (auto l : *base[*idx]) {
          score[*idx] += inverseLearningRate(l);
        }
        score[*idx] *= glue[*idx];
      }
    } else if (solver.getOptions().forget_strategy ==
               Options::LiteralScore::Looseness) {
      for (auto idx{free_cl_indices.bbegin()}; idx != free_cl_indices.bend();
           ++idx) {
        score[*idx] = 0;

        //        std::cout << "compute (1) score for " << *base[*idx] <<
        //        std::endl;

        for (auto l : *base[*idx]) {
          score[*idx] += looseness(l);
        }
      }
    } else if (solver.getOptions().forget_strategy ==
               Options::LiteralScore::LoosenessOverActivity) {

      //  std::cout << 3124 with l/a\n";

      for (auto idx{free_cl_indices.bbegin()}; idx != free_cl_indices.bend();
           ++idx) {
        score[*idx] = 0;

        //        std::cout << "compute (2) score for " << *base[*idx] <<
        //        std::endl;

        for (auto l : *base[*idx]) {
          score[*idx] += loosenessOverActivity(l);
        }
      }
    } else if (solver.getOptions().forget_strategy ==
               Options::LiteralScore::LoosenessOverLearningRate) {

      //  std::cout << 3124 with l/a\n";

      for (auto idx{free_cl_indices.bbegin()}; idx != free_cl_indices.bend();
           ++idx) {
        score[*idx] = 0;

        //        std::cout << "compute (2) score for " << *base[*idx] <<
        //        std::endl;

        for (auto l : *base[*idx]) {
          score[*idx] += loosenessOverLearningRate(l);
        }
      }
    } else if (solver.getOptions().forget_strategy ==
               Options::LiteralScore::Activity) {
      for (auto idx{free_cl_indices.bbegin()}; idx != free_cl_indices.bend();
           ++idx) {
        score[*idx] = 0;

        //        std::cout << "compute (2) score for " << *base[*idx] <<
        //        std::endl;

          double bool_score{0};
          double num_score{0};
          int n_bool{0};
          int n_num{0};
        for (auto l : *base[*idx]) {
          score[*idx] += inverseActivity(l);
            if(l.isNumeric()) {
                ++n_num;
                num_score += inverseActivity(l);
            } else {
                ++n_bool;
                bool_score += inverseActivity(l);
            }
        }
          
//          if(n_bool > 0 and n_num > 0) {
//              std::cout << "score factor: " << ((bool_score * static_cast<double>(n_num)) / (static_cast<double>(n_bool) * num_score)) ;
//              
//              std::cout << " count ratio: " << static_cast<double>(n_bool) / static_cast<double>(n_num) << std::endl;
//          }
          
//          std::cout << "b: " ;
//          if(n_bool > 0)
//              std::cout << bool_score / static_cast<double>(n_bool) << " (" << n_bool << ") || n: " ;
//          else
//              std::cout << "0 || n: ";
//          if(n_num > 0)
//              std::cout << num_score / static_cast<double>(n_num) << " (" << n_num << ")\n";
//          else
//              std::cout << "0\n";
      }
    } else if (solver.getOptions().forget_strategy ==
               Options::LiteralScore::LearningRate) {
      for (auto idx{free_cl_indices.bbegin()}; idx != free_cl_indices.bend();
           ++idx) {
        score[*idx] = 0;

        //        std::cout << "compute (2) score for " << *base[*idx] <<
        //        std::endl;

        for (auto l : *base[*idx]) {
          score[*idx] += inverseLearningRate(l);
        }
      }
    }

    if (solver.getOptions().forget_strategy == Options::LiteralScore::Size) {
        std::sort(free_cl_indices.bbegin(), free_cl_indices.bend(),
                  [&](const int i, const int j) {
            return (base[i]->size() > base[j]->size());
        });
    } else if (solver.getOptions().forget_strategy == Options::LiteralScore::Glue) {
        std::sort(free_cl_indices.bbegin(), free_cl_indices.bend(),
                  [&](const int i, const int j) {
            return (glue[base[i]->id] > glue[base[j]->id] || (glue[base[i]->id] == glue[base[j]->id] and base[i]->size() > base[j]->size()));
        });
    } else {
        std::sort(free_cl_indices.bbegin(), free_cl_indices.bend(),
                  [&](const int i, const int j) {
            return (score[base[i]->id] > score[base[j]->id]);
        });
    }
    
    free_cl_indices.re_index(free_cl_indices.bbegin(), free_cl_indices.bend());
    
    
//        std::cout << "ordering, from worst to best:\n";
//        for(auto i{free_cl_indices.bbegin()}; i!=free_cl_indices.bend(); ++i) {
//            std::cout << std::setw(4) << glue[*i] << " " << std::setw(8) << std::setprecision(3) << score[*i] << " " << *base[*i] << std::endl;
//        }
    
    auto target_size =
    static_cast<size_t>(static_cast<double>(numLearnt()) *
                        (1.0 - solver.getOptions().forgetfulness));
    
    while (numLearnt() > target_size) {
        forget_worst();
    }
    
//        std::cout << "\n===>\n";
//        for(auto i{free_cl_indices.bbegin()}; i!=free_cl_indices.bend(); ++i) {
//            std::cout << std::setw(4) << glue[*i] << " " << std::setw(8) << std::setprecision(3) << score[*i] << " " << *base[*i] << std::endl;
//        }
    
#ifdef DBG_WATCHERS
    verifyWatchers("after forget");
#endif
}

template <typename T>
template <typename iter>
Clause<T> *ClauseBase<T>::add(const iter first, const iter last,
                              const bool learnt, const int glue_score) {
    
#ifdef DBG_WATCHERS
    verifyWatchers("before add");
#endif
    
//    if(learnt == false) {
//        std::cout << "ADD TRUE CLAUSE!\n";
//    }
    
    //  std::cout << "add clause (" << (learnt ? "learnt" : "base" ) << ")";
    //  for (auto l{first}; l != last; ++l) {
    //    std::cout << " " << *l;
    //  }
    //  std::cout << std::endl;
    
    Clause<T> *c = nullptr;
    
    if (first == last) {
        throw Failure<T>({this, Constant::NoHint});
    }
    
    if (first + 1 == last) {
        //    assign(*first, Constant::NewNoReason<T>);
        solver.set(*first); //, Constant::GroundFact<T>);
    } else {
        if (not free_cl_indices.empty()) {
            int id{free_cl_indices.back()};
            if (learnt)
                free_cl_indices.remove_back(id);
            else
                free_cl_indices.remove_front(id);
            c = base[id];
            glue[id] = glue_score;
            assert(c->empty());
        } else {
            int id{static_cast<int>(base.size())};
            c = new Clause<T>(id);
            base.push_back(c);
            score.push_back(0);
            glue.push_back(glue_score);
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

        ////        std::cout << c->id << std::endl;
        //        if(c->id == 3124) {
        //            std::cout << "(re)use id 3124 for a clause of size " <<
        //            c->size() << " @" << solver.num_choicepoints << std::endl;
        //
        //            for(auto q : *c) {
        //                std::cout << q << std::endl;
        //            }
        //            std::cout << std::endl;
        //
        //            exit(1);
        //
        //        }

        total_size += c->size();
        
        if (learnt) {
            Literal<T> l{(*c)[0]};
            //          assign(l, {this, c->id});
            
#ifdef DBGP0
            if (solver.level() <= solver.init_level) {
                std::cout << "pruning @lvl " << solver.init_level << std::endl;
            }
#endif
            if (solver.level() <= solver.init_level) {
                makeUnforgettable(c);
            }
            
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
    
    if (l == Contradiction<T>) {
        
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

//template <typename T> int ClauseBase<T>::getType() const {
//  return CLAUSEEXPL;
//}

template <typename T>
std::ostream &ClauseBase<T>::display(std::ostream &os) const {
    os << "clause base";
    return os;
}

template <typename T>
std::ostream &ClauseBase<T>::displayClauses(std::ostream &os) const {
    for (auto i{free_cl_indices.fbegin()}; i != free_cl_indices.fend(); ++i) {
        std::cout << *base[*i] << std::endl;
    }
    std::cout << "++\n";
    for (auto i{free_cl_indices.bbegin()}; i != free_cl_indices.bend(); ++i) {
        std::cout << *base[*i] << std::endl;
    }
    std::cout << "--\n";
    
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
#ifdef NEW_WATCHERS
    os << "bool watchers:\n";
    for (size_t x{0}; x < solver.boolean.size(); ++x) {
        auto l{solver.boolean.getLiteral(true, x)};
        if (not boolean_watch[l].empty()) {
            os << l << " is watched in";
            for (auto cl : boolean_watch[l]) {
                os << " " << *cl;
            }
            os << std::endl;
        }
        if (not boolean_watch[~l].empty()) {
            os << ~l << " is watched in";
            for (auto cl : boolean_watch[~l]) {
                os << " " << *cl;
            }
            os << std::endl;
        }
    }
    os << "numeric watch-lists (not sorted:" << watchers.capacity()
    << "):\n";
    auto i{0};
    for (auto &wl : watchers) {
        os << i++ << ": ";
        for (auto cl : wl) {
            os << " " << *cl; // cl->id();
        }
        os << std::endl;
    }
    os << "numeric watch-lists (sorted):\n";
    for (size_t x{0}; x < solver.numeric.size(); ++x) {
        auto l = lb<T>(x);
        if (numeric_watch[l].size() > 1) {
            os << l << " is watched in";
            for (auto wl : numeric_watch[l]) {
                os << " <" << wl.first << ":" << wl.second << ">";
                //        os << "[" << wl.first << "]:";
                //        for (auto cl : unsorted_numeric_watch[wl.second]) {
                //          os << " " << *cl;
                //        }
            }
            os << std::endl;
        }
        auto u = ub<T>(x);
        if (numeric_watch[u].size() > 1) {
            os << u << " is watched in";
            for (auto wl : numeric_watch[u]) {
                os << " <" << wl.first << ":" << wl.second << ">";
                //        os << "[" << wl.first << "]:";
                //        for (auto cl : unsorted_numeric_watch[wl.second]) {
                //          os << " " << *cl;
                //        }
            }
            os << std::endl;
        }
    }
#else
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
        auto u = ub<T>(x);
        if (not watch[NUMERIC][u].empty()) {
            os << u << " is watched in";
            for (auto cl : watch[NUMERIC][u]) {
                os << " " << *cl;
            }
            os << std::endl;
        }
    }
#endif
    
    return os;
}

#ifdef DBG_WATCHERS
template <typename T>
void ClauseBase<T>::verifyWatchers(const char *msg) const {
    
    //    std::vector<>
    
    //    std::cout << std::endl << "watches:\n";
    //    displayWatchStruct(std::cout);
    //    std::cout << std::endl;
    
    int i{0};
    for (auto cl : base) {
        if (cl->id != i) {
          std::cout << msg << "indexing error @" << solver.num_choicepoints
                    << std::endl;
          exit(1);
        }
        ++i;
    }
    
    size_t num_watchers{0};
    for (size_t x{1}; x < solver.boolean.size(); ++x) {
        Literal<T> l{solver.boolean.getLiteral(true, x)};
#ifdef NEW_WATCHERS
        num_watchers += watchers[boolean_watch[l]].size();
        if (not watchers[boolean_watch[l]].empty()) {
            for (auto cl : watchers[boolean_watch[l]]) {
                
                if (cl->size() < 2) {
                  std::cout << msg << "(*): error empty clause watching " << l
                            << " @" << solver.num_choicepoints << std::endl;
                  exit(1);
                }
                
                if (free_cl_indices.has(cl->id)) {
                  std::cout << msg << "(*): " << *cl << "'s id is free @"
                            << solver.num_choicepoints << std::endl;
                  exit(1);
                }
                
                if (cl->watched(0) != ~l and cl->watched(1) != ~l) {
                  std::cout << msg << "(*): error on clause " << cl->id
                            << " -- " << *cl << " on " << ~l
                            << "'s watch-list but marked as watching "
                            << cl->watched(0) << " and " << cl->watched(1)
                            << " @" << solver.num_choicepoints << std::endl;
                  exit(1);
                }
            }
        }
        if (not watchers[boolean_watch[~l]].empty()) {
            num_watchers += watchers[boolean_watch[~l]].size();
            for (auto cl : watchers[boolean_watch[~l]]) {
                if (cl->size() < 2) {
                  std::cout << msg << "(*): error empty clause watching " << l
                            << " @" << solver.num_choicepoints << std::endl;
                  exit(1);
                }
                
                if (free_cl_indices.has(cl->id)) {
                  std::cout << msg << "(*): " << *cl << "'s id is free @"
                            << solver.num_choicepoints << std::endl;
                  exit(1);
                }
                
                if (cl->watched(0) != l and cl->watched(1) != l) {
                  std::cout << msg << "(*): error on clause " << cl->id
                            << " -- " << *cl << " on " << l
                            << "'s watch-list but marked as watching "
                            << cl->watched(0) << " and " << cl->watched(1)
                            << " @" << solver.num_choicepoints << std::endl;
                  exit(1);
                }
            }
        }
#else
        num_watchers += watch[BOOLEAN][l].size();
        if (not watch[BOOLEAN][l].empty()) {
            for (auto cl : watch[BOOLEAN][l]) {
                
                if (cl->size() < 2) {
                  std::cout << msg << ": error empty clause watching " << l
                            << " @" << solver.num_choicepoints << std::endl;
                  exit(1);
                }
                
                if (free_cl_indices.has(cl->id)) {
                  std::cout << msg << ": " << *cl << "'s id is free @"
                            << solver.num_choicepoints << std::endl;
                  exit(1);
                }
                
                if (cl->watched(0) != ~l and cl->watched(1) != ~l) {
                  std::cout << msg << ": error on clause " << cl->id << " -- "
                            << *cl << " watching " << ~l << " @"
                            << solver.num_choicepoints << std::endl;
                  exit(1);
                }
            }
        }
        if (not watch[BOOLEAN][~l].empty()) {
            num_watchers += watch[BOOLEAN][~l].size();
            for (auto cl : watch[BOOLEAN][~l]) {
                if (cl->size() < 2) {
                  std::cout << msg << ": error empty clause watching " << l
                            << " @" << solver.num_choicepoints << std::endl;
                  exit(1);
                }
                
                if (free_cl_indices.has(cl->id)) {
                  std::cout << msg << ": " << *cl << "'s id is free @"
                            << solver.num_choicepoints << std::endl;
                  exit(1);
                }
                
                if (cl->watched(0) != l and cl->watched(1) != l) {
                  std::cout << msg << ": error on clause " << cl->id << " -- "
                            << *cl << " watching " << l << " @"
                            << solver.num_choicepoints << std::endl;
                  exit(1);
                }
            }
        }
#endif
    }
    for (size_t x{0}; x < solver.numeric.size(); ++x) {
        auto l{lb<T>(x)};
        auto p{ub<T>(x)};
        
#ifdef NEW_WATCHERS
        T prev{-Constant::Infinity<T>};
        for (auto &wl : numeric_watch[l]) {
            if (wl.first <= prev) {
              std::cout << msg << "(*): error numeric watches not sorted @"
                        << solver.num_choicepoints << std::endl;
              exit(1);
            }
            prev = wl.first;
            
            auto &watches{watchers[wl.second]};
            num_watchers += watches.size();
            if (not watches.empty()) {
                for (auto cl : watches) {
                    if (cl->size() < 2) {
                      std::cout << msg << "(*): error empty clause watching "
                                << l << " @" << solver.num_choicepoints
                                << std::endl;
                      exit(1);
                    }
                    
                    if (free_cl_indices.has(cl->id)) {
                      std::cout << msg << "(*): " << *cl << "'s id is free @"
                                << solver.num_choicepoints << std::endl;
                      exit(1);
                    }
                    
                    if (not((cl->watched(0).sameVariable(l) and
                             cl->watched(0).sign() != l.sign()) or
                            (cl->watched(1).sameVariable(l) and
                             cl->watched(1).sign() != l.sign()))) {
                      std::cout << msg << "(*): error on clause " << cl->id
                                << " -- " << *cl << " on " << l
                                << "'s watch-list but marked as watching "
                                << cl->watched(0) << " and " << cl->watched(1)
                                << " @" << solver.num_choicepoints << std::endl;
                      exit(1);
                    }
                }
            } else if (wl.first != Constant::Infinity<T>) {
              std::cout << msg
                        << "(*): empty watch list guarded by a finite value "
                        << l << " @" << solver.num_choicepoints << std::endl;
              exit(1);
            }
        }
        prev = -Constant::Infinity<T>;
        for (auto &wl : numeric_watch[p]) {
            if (wl.first <= prev) {
              std::cout << msg << "(*): error numeric watches not sorted "
                        << " @" << solver.num_choicepoints << std::endl;
              exit(1);
            }
            prev = wl.first;
            
            auto &watches{watchers[wl.second]};
            num_watchers += watches.size();
            if (not watches.empty()) {
                for (auto cl : watches) {
                    if (cl->size() < 2) {
                      std::cout << msg << "(*): error empty clause watching "
                                << l << " @" << solver.num_choicepoints
                                << std::endl;
                      exit(1);
                    }
                    
                    if (free_cl_indices.has(cl->id)) {
                      std::cout << msg << "(*): " << *cl << "'s id is free @"
                                << solver.num_choicepoints << std::endl;
                      exit(1);
                    }
                    
                    if (not((cl->watched(0).sameVariable(p) and
                             cl->watched(0).sign() != p.sign()) or
                            (cl->watched(1).sameVariable(p) and
                             cl->watched(1).sign() != p.sign()))) {
                      std::cout << msg << "(*): error on clause " << cl->id
                                << " -- " << *cl << " on " << p
                                << "'s watch-list but marked as watching "
                                << cl->watched(0) << " and " << cl->watched(1)
                                << " @" << solver.num_choicepoints << std::endl;
                      exit(1);
                    }
                }
            } else if (wl.first != Constant::Infinity<T>) {
              std::cout << msg
                        << "(*): empty watch list guarded by a finite value "
                        << l << " @" << solver.num_choicepoints << std::endl;
              exit(1);
            }
        }
#else
        num_watchers += watch[NUMERIC][l].size();
        if (not watch[NUMERIC][l].empty()) {
            for (auto cl : watch[NUMERIC][l]) {
                if (cl->size() < 2) {
                  std::cout << msg << ": error empty clause watching " << l
                            << " @" << solver.num_choicepoints << std::endl;
                  exit(1);
                }
                
                if (free_cl_indices.has(cl->id)) {
                  std::cout << *cl << "'s id is free @"
                            << solver.num_choicepoints << std::endl;
                  exit(1);
                }
                
                if (not((cl->watched(0).sameVariable(l) and
                         cl->watched(0).sign() != l.sign()) or
                        (cl->watched(1).sameVariable(l) and
                         cl->watched(1).sign() != l.sign()))) {
                  std::cout << msg << ": error on clause " << cl->id << " -- "
                            << *cl << " watching " << l << " (" << info_t(l)
                            << ")"
                            << " @" << solver.num_choicepoints << std::endl;
                  exit(1);
                }
            }
        }
        
        if (not watch[NUMERIC][p].empty()) {
            num_watchers += watch[NUMERIC][p].size();
            for (auto cl : watch[NUMERIC][p]) {
                if (cl->size() < 2) {
                  std::cout << msg << ": error empty clause watching " << p
                            << " @" << solver.num_choicepoints << std::endl;
                  exit(1);
                }
                
                if (free_cl_indices.has(cl->id)) {
                  std::cout << *cl << "'s id is free  @"
                            << solver.num_choicepoints << std::endl;
                  exit(1);
                }
                
                if (not((cl->watched(0).sameVariable(p) and
                         cl->watched(0).sign() != p.sign()) or
                        (cl->watched(1).sameVariable(p) and
                         cl->watched(1).sign() != p.sign()))) {
                  std::cout << msg << ": error on clause " << cl->id << " -- "
                            << *cl << " watching " << p << " (" << info_t(p)
                            << ")"
                            << " @" << solver.num_choicepoints << std::endl;
                  exit(1);
                }
            }
        }
#endif
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

