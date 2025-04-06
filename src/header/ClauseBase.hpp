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
//#define DBGLRBF

//#define NEW_WATCHERS
#define CORRECT_ORDER

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

    // #ifdef NEW_WATCHERS
    std::list<std::pair<T, index_t>>::iterator
    getNumWatchList(const Literal<T> l);
//    index_t getWatchList(const Literal<T> l);
    index_t getOrCreateWatchList(const Literal<T> l);
    // #else
    //     std::vector<Clause<T> *> &getWatchList(const Literal<T> l);
    //     std::vector<Clause<T> *> &getOrCreateWatchList(const Literal<T> l);
    // #endif
    //  set the 'r'-th watcher of clause 'cl' to be its 'i'-th literal
    void set_watcher(const int r, const index_t i, Clause<T> *cl);
    
    //  // set the 'r'-th watcher of clause 'cl' to be its 'i'-th literal (and
    //  order
    //  // the watch-list of this literal)
    //  void set_watcher_numeric(const int r, const index_t i, Clause<T> *cl);
    
    // works for both numeric and Boolean, does not order numeric watch-lists
//    void unit_propagate(const Literal<T> l);

    // #ifdef NEW_WATCHERS
    //  helper to avoid code duplication
    void up(const Literal<T> l, const index_t widx);//std::vector<Clause<T> *> &watches);
    
    // works only for Boolean literals
    void unit_propagate_boolean(const Literal<T> l);
    
    // works only for numeric literals
    void unit_propagate_numeric(const Literal<T> l);
    // #else
    //     // works only for Boolean literals
    //     void unit_propagate_boolean(const Literal<T> l);
    //
    //     // works only for numeric literals, orders numeric watch-lists
    //     void unit_propagate_numeric(const Literal<T> l);
    // #endif

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

    // #ifdef NEW_WATCHERS
    //  a watch-list for every Boolean literal
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
    // #else
    //     // the watch lists for every literel (watch[BOOLEAN] and
    //     watch[NUMERIC]) std::vector<std::vector<Clause<T> *>> watch[2];
    // #endif

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
    size_t num_boolean{0};
    size_t num_numeric{0};
    //@}
    
    // indexing helper for the watch lists
    static const bool BOOLEAN{false};
    static const bool NUMERIC{true};
    
    // type helper for the literals -> return l.isNumeric() == NUMERIC is l is
    // numeric
    static bool litType(Literal<T> l) { return l.isNumeric(); }
    
#ifdef DBG_WATCHERS
    void verifyWatchers(const char *msg) const;
    void checkNumWatchers(const char *msg) const;
#endif
};


////// NEW CLAUSES
///

template <typename T>
ClauseBase<T>::ClauseBase(Solver<T> &c)
: solver(c), handlerToken(solver.SearchRestarted.subscribe_handled(
                                                                   [this](const bool) { this->forget(); })) {
    
    Constraint<T>::priority = Priority::Low;

    // #ifdef NEW_WATCHERS
    watchers.resize(1);
    free_wl_indices.reserve(1);
    // #endif
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
        unit_propagate_numeric(p);
    }
    clearTriggers();
}

template <typename T> void ClauseBase<T>::clearTriggers() {
    triggered_bounds.clear();
}

template <typename T>
bool ClauseBase<T>::notify(const Literal<T> l, const int) {
    
    assert(l.isNumeric());
   
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
    return free_cl_indices.backsize() + free_cl_indices.frontsize();
}

template <typename T> size_t ClauseBase<T>::volume() const {
    return total_size;
}

template <typename T> Clause<T> *ClauseBase<T>::operator[](const index_t i) const {
    if (i >= free_cl_indices.capacity() or free_cl_indices.has(i))
        return NULL;
    return base[i];
}

template <typename T> Clause<T> *ClauseBase<T>::back() {
    return base[*(free_cl_indices.bbegin())];
}

template <typename T> void ClauseBase<T>::newBooleanVar() {
  boolean_watch.resize(2 * solver.boolean.size());
  *boolean_watch.rbegin() = watchers.size() + 1;
  *(boolean_watch.rbegin() + 1) = watchers.size();
  watchers.resize(watchers.size() + 2);
}

template <typename T> void ClauseBase<T>::newNumericVar() {

  auto sz{numeric_watch.size()};
  numeric_watch.resize(2 * solver.numeric.size());
  while (sz < numeric_watch.size()) {
      
#ifdef CORRECT_ORDER
//      numeric_watch[sz].insert(numeric_watch[sz].end(),
//                               {-Constant::Infinity<T>, 0});
#else
    numeric_watch[sz].insert(numeric_watch[sz].end(),
                             {Constant::Infinity<T>, 0});
#endif
    ++sz;
  }
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


template <typename T> Clause<T> *ClauseBase<T>::consistent() {
    
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
void ClauseBase<T>::up(const Literal<T> l, const index_t widx){ //std::vector<Clause<T> *> &watches) {
    
#ifdef DBG_WATCHERS
    verifyWatchers("before up");
#endif

    size_t k{watchers[widx].size()};
    while(k-->0) {

        auto cl{watchers[widx][k]};
        
#ifdef DBG_TRACE
        if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
            std::cout << " info_t = " << info_t(l) << " /watched by " << *cl
            << std::endl;
        }
#endif

        assert(cl->watch_index(1) < cl->size());
 
        bool watch_rank{cl->watch_rank(l)};
        index_t idx{cl->watched_index[watch_rank]};
        Literal<T> other{cl->watched(1 - watch_rank)};
         
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
 
                    assert(watchers[widx].size() >= 0);
   
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
    up(l, boolean_watch[l]);
}

template <typename T>
void ClauseBase<T>::unit_propagate_numeric(const Literal<T> l) {
    
    
//    std::cout << "UP " << l << std::endl;
    
    
#ifdef CORRECT_ORDER
    
    auto wl{numeric_watch[l].begin()};
    auto the_end{numeric_watch[l].end()};
    
//    std::cout << " (" << wl->first << ")\n";
    
    while(wl != the_end and l.value() <= wl->first) {
        
//        std::cout << " *" << wl->first << std::endl;
        
        auto next{wl};
        ++next;
        up(l, wl->second);
        if (watchers[wl->second].empty()) {
          free_wl_indices.add(wl->second);
          numeric_watch[l].erase(wl);
        }
        wl = next;
    }
    
//    std::cout << "ignore";
//    while(wl != numeric_watch[l].end()) {
//        std::cout << " " << wl->first ;
//        ++wl;
//    }
//    std::cout << std::endl;
    
#else
    
    auto wl{getNumWatchList(l)};
    while (wl->first < Constant::Infinity<T>) {
        
//        std::cout << " *" << wl->first << std::endl;
        
        auto next{wl};
        ++next;
        up(l, wl->second);
        if (watchers[wl->second].empty()) {
          free_wl_indices.add(wl->second);
          numeric_watch[l].erase(wl);
        }
        wl = next;
    }
    
#endif
    
}
// #endif

//template <typename T> void ClauseBase<T>::unit_propagate(const Literal<T> l) {
//    
//    ++num_up;
//    
//#ifdef DBG_WATCHERS
//    verifyWatchers("before UP");
//#endif
//    
//#ifdef DBG_TRACE
//    if (DBG_CLBOUND and (DBG_TRACE & UNITPROPAGATION)) {
//        std::cout << "\nunit propagate " << l << std::endl;
//    }
//#endif
//
//    // #ifdef NEW_WATCHERS
//    if (l.isNumeric()) {
//        unit_propagate_numeric(l);
//    } else {
//        unit_propagate_boolean(l);
//    }
//
//#ifdef DBG_WATCHERS
//    verifyWatchers("after UP");
//#endif
//}


// #ifdef NEW_WATCHERS
template <typename T>
std::list<std::pair<T, index_t>>::iterator
ClauseBase<T>::getNumWatchList(const Literal<T> l) {
    
    auto &wl(numeric_watch[l]);
    auto li{wl.begin()}; // list of watch-lists, ordered by
          // no need to check for the end because there is a sentinel with
    
#ifdef CORRECT_ORDER
    
    while (li->first > l.value()) { ++li; if(li == wl.end()) { std::cout << "wtf?\n"; exit(1); } }
    
#else
    
    while (li->first < l.value()) { ++li; }
   
#endif
    
    return li;
}
// #endif

//// #ifdef NEW_WATCHERS
////  return the watch-list of literal l, or the iterator after its expected
////  position if l is numeric and has no watch-list yet
//template <typename T>
//index_t ClauseBase<T>::getWatchList(const Literal<T> l) {
//    if (l.isNumeric()) {
//        return getNumWatchList(l)->second;
//    } else {
//        return boolean_watch[l];
//    }
//}


// return the watch-list of literal l, or the iterator after its expected
// position if l is numeric and has no watch-list yet
template <typename T>
index_t
ClauseBase<T>::getOrCreateWatchList(const Literal<T> l) {
    if (l.isNumeric()) {
        assert(l.value() < Constant::Infinity<T>);
        
#ifdef CORRECT_ORDER
        
//        std::cout << " add " << l << " into";
//        
//        for(auto wl : numeric_watch[l]) {
//            std::cout << " (" << wl.first << "|" << wl.second << ")";
//        }
//        std::cout << std::endl;
        
        auto wl{numeric_watch[l].begin()};
        auto the_end{numeric_watch[l].end()};
        
//        std::cout << " 1=> [" << std::distance(numeric_watch[l].begin(), wl) << "]\n";
        
        while(wl != the_end and wl->first > l.value()) {
            ++wl;
        }
        
//        std::cout << " 2=> [" << std::distance(numeric_watch[l].begin(), wl) << "]\n";
        
        if (wl != the_end and wl->first == l.value()) {
            
//            std::cout << " --> (" << wl->first << "|" << wl->second << ")\n";
            
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
             
//            std::cout << " ==> [" << std::distance(numeric_watch[l].begin(), wl) << "] (" << wl->first << "|" << wl->second << ")\n";
            
            numeric_watch[l].insert(wl, {l.value(), idx});
            
//            std::cout << " :: ";
//            for(auto wli : numeric_watch[l]) {
//                std::cout << " (" << wli.first << "|" << wli.second << ")";
//            }
//            std::cout << std::endl;
            
            return idx;
        }
   
#else
        
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
#endif
        
    } else {
        return boolean_watch[l];
    }
}


template <typename T>
void ClauseBase<T>::set_watcher(const int r, const index_t i, Clause<T> *cl) {
    
    cl->watched_index[r] = i;
    Literal<T> l{~((*cl)[i])};

    
//    if(solver.num_cons_propagations >= 4748414) {
//        if(l.isNumeric()) {
//            std::cout << "set watch " << r << " at idx " << i << ": " << l << "\ncurrent WL";
//            for(auto& wl : numeric_watch[l]) {
//                std::cout << " (" << wl.first << "|" << wl.second << "|" << watchers[wl.second].size() << ")";
//            }
//            std::cout << std::endl;
//        }
//    }
    
    auto widx{getOrCreateWatchList(l)};
    
//    if(solver.num_cons_propagations >= 4748414) {
//        if(l.isNumeric()) {
//            std::cout << "after goc WL";
//            for(auto& wl : numeric_watch[l]) {
//                std::cout << " (" << wl.first << "|" << wl.second << "|" << watchers[wl.second].size() << ")";
//            }
//            std::cout << std::endl << " ==> " << widx << std::endl;
//        }
//    }
    
    auto &watches{watchers[widx]};

    if (l.isNumeric()) {
        watches.push_back(cl);
    } else {
        watches.push_back(cl);
    }
}


template <typename T> void ClauseBase<T>::makeUnforgettable(Clause<T> *cl) {
    if (free_cl_indices.isback(cl->id)) {
        free_cl_indices.add(cl->id);
        free_cl_indices.remove_front(cl->id);
    }
}

template <typename T> void ClauseBase<T>::forget(Clause<T> *cl) {
  for (auto r{0}; r < 2; ++r) {
    auto wl{~(cl->watched(r))};

    if (wl.isNumeric()) {
      auto watchlist{getNumWatchList(wl)};
      auto &watches{watchers[watchlist->second]};
      for (auto c{watches.begin()}; c != watches.end(); ++c)
        if (*c == cl) {
          swap(*c, watches.back());
          watches.pop_back();
          break;
        }
      if (watches.empty()) {
        free_wl_indices.add(watchlist->second);
        numeric_watch[wl].erase(watchlist);
      }
    } else {
//      auto &watches{watchers[getWatchList(wl)]}; boolean_watch[l]
        auto &watches{watchers[boolean_watch[wl]]};
      for (auto c{watches.begin()}; c != watches.end(); ++c)
        if (*c == cl) {
          swap(*c, watches.back());
          watches.pop_back();
          break;
        }
    }
  }

    free_cl_indices.add(cl->id);
    
    total_size -= cl->size();
    for(auto l : *cl) {
        if(l.isNumeric())
            --num_numeric;
        else
            --num_boolean;
    }
    
    cl->clear();
}

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
    while (free_cl_indices.backsize() > 0) {
        forget_worst();
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
    
    auto bias{1.0 / solver.getOptions().literal_bias};
    
#ifdef DBG_WATCHERS
    verifyWatchers("before forget");
#endif

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
        for (auto l : *base[*idx]) {
          score[*idx] += looseness(l);
        }
      }
    } else if (solver.getOptions().forget_strategy ==
               Options::LiteralScore::LoosenessOverActivity) {
      for (auto idx{free_cl_indices.bbegin()}; idx != free_cl_indices.bend();
           ++idx) {
        score[*idx] = 0;
        for (auto l : *base[*idx]) {
          score[*idx] += loosenessOverActivity(l);
        }
      }
    } else if (solver.getOptions().forget_strategy ==
               Options::LiteralScore::LoosenessOverLearningRate) {
      for (auto idx{free_cl_indices.bbegin()}; idx != free_cl_indices.bend();
           ++idx) {
        score[*idx] = 0;
        for (auto l : *base[*idx]) {
          score[*idx] += loosenessOverLearningRate(l);
        }
      }
    } else if (solver.getOptions().forget_strategy ==
               Options::LiteralScore::Activity) {
      for (auto idx{free_cl_indices.bbegin()}; idx != free_cl_indices.bend();
           ++idx) {
        score[*idx] = 0;
        for (auto l : *base[*idx]) {
            auto s{inverseActivity(l)};
            
            if(not l.isNumeric())
                s *= bias;
            
            score[*idx] += s;
        }
      }
    } else if (solver.getOptions().forget_strategy ==
               Options::LiteralScore::LearningRate) {
      for (auto idx{free_cl_indices.bbegin()}; idx != free_cl_indices.bend();
           ++idx) {
        score[*idx] = 0;

#ifdef DBGLRBF
                    double bool_score{0};
                    double num_score{0};
                    int n_bool{0};
                    int n_num{0};
#endif
          
        for (auto l : *base[*idx]) {
            
            auto s{inverseLearningRate(l)};
            
#ifdef DBGLRBF
            if(l.isNumeric()) {
                ++n_num;
                num_score += s;
            } else {
                ++n_bool;
                bool_score += s;
            }
#endif
            
            if(not l.isNumeric())
                s *= bias;
            
          score[*idx] += s;
            
            
        }
          
#ifdef DBGLRBF
                    std::cout << "b: " ;
                    if(n_bool > 0)
                        std::cout << bool_score / static_cast<double>(n_bool) << " (" << n_bool << ") || n: " ;
                    else
                        std::cout << "0 || n: ";
                    if(n_num > 0)
                        std::cout << num_score / static_cast<double>(n_num) << " (" << n_num << ")\n";
                    else
                        std::cout << "0\n";
#endif
      }
    } else if (solver.getOptions().forget_strategy == Options::LiteralScore::Size) {
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

//    if(solver.num_cons_propagations >= 4748414) {
//        std::cout << "add clause (" << (learnt ? "learnt" : "base" ) << ")";
//        for (auto l{first}; l != last; ++l) {
//            std::cout << " " << *l;
//        }
//        std::cout << std::endl;
//    }
    
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
            if(l->isNumeric()) {
                ++num_numeric;
            } else {
                ++num_boolean;
            }
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
  // #ifdef NEW_WATCHERS
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
    // #else
    //     for (size_t x{0}; x < solver.boolean.size(); ++x) {
    //         auto l{solver.boolean.getLiteral(true, x)};
    //         if (not watch[BOOLEAN][l].empty()) {
    //             os << l << " is watched in";
    //             for (auto cl : watch[BOOLEAN][l]) {
    //                 os << " " << *cl;
    //             }
    //             os << std::endl;
    //         }
    //         if (not watch[BOOLEAN][~l].empty()) {
    //             os << ~l << " is watched in";
    //             for (auto cl : watch[BOOLEAN][~l]) {
    //                 os << " " << *cl;
    //             }
    //             os << std::endl;
    //         }
    //     }
    //     for (size_t x{0}; x < solver.numeric.size(); ++x) {
    //         auto l = lb<T>(x);
    //         if (not watch[NUMERIC][l].empty()) {
    //             os << l << " is watched in";
    //             for (auto cl : watch[NUMERIC][l]) {
    //                 os << " " << *cl;
    //             }
    //             os << std::endl;
    //         }
    //         auto u = ub<T>(x);
    //         if (not watch[NUMERIC][u].empty()) {
    //             os << u << " is watched in";
    //             for (auto cl : watch[NUMERIC][u]) {
    //                 os << " " << *cl;
    //             }
    //             os << std::endl;
    //         }
    //     }
    // #endif

    return os;
}

#ifdef DBG_WATCHERS
template <typename T>
void ClauseBase<T>::checkNumWatchers(const char *msg) const {
    size_t num_watchers{0};
    for (size_t x{1}; x < solver.boolean.size(); ++x) {
        Literal<T> l{solver.boolean.getLiteral(true, x)};
        num_watchers += watchers[boolean_watch[l]].size();
        num_watchers += watchers[boolean_watch[~l]].size();
    }
    for (size_t x{0}; x < solver.numeric.size(); ++x) {
        auto l{lb<T>(x)};
        auto p{ub<T>(x)};
        for (auto &wl : numeric_watch[l]) {
            auto &watches{watchers[wl.second]};
            num_watchers += watches.size();
        }
        for (auto &wl : numeric_watch[p]) {
            auto &watches{watchers[wl.second]};
            num_watchers += watches.size();
        }
    }
    if (num_watchers != 2 * size()) {
        std::cout << msg << ": wrong number of watchers (" << solver.num_cons_propagations << ") !\n";
        exit(1);
    }
}

template <typename T>
void ClauseBase<T>::verifyWatchers(const char *msg) const {
    checkNumWatchers(msg);
 
//    
//
//    int i{0};
//    for (auto cl : base) {
//        if (cl->id != i) {
//          std::cout << msg << "indexing error @" << solver.num_choicepoints
//                    << std::endl;
//          exit(1);
//        }
//        ++i;
//    }
//    
//    size_t num_watchers{0};
//    for (size_t x{1}; x < solver.boolean.size(); ++x) {
//        Literal<T> l{solver.boolean.getLiteral(true, x)};
//        // #ifdef NEW_WATCHERS
//        num_watchers += watchers[boolean_watch[l]].size();
//        if (not watchers[boolean_watch[l]].empty()) {
//            for (auto cl : watchers[boolean_watch[l]]) {
//                
//                if (cl->size() < 2) {
//                  std::cout << msg << "(*): error empty clause watching " << l
//                            << " @" << solver.num_choicepoints << std::endl;
//                  exit(1);
//                }
//                
//                if (free_cl_indices.has(cl->id)) {
//                  std::cout << msg << "(*): " << *cl << "'s id is free @"
//                            << solver.num_choicepoints << std::endl;
//                  exit(1);
//                }
//                
//                if (cl->watched(0) != ~l and cl->watched(1) != ~l) {
//                  std::cout << msg << "(*): error on clause " << cl->id
//                            << " -- " << *cl << " on " << ~l
//                            << "'s watch-list but marked as watching "
//                            << cl->watched(0) << " and " << cl->watched(1)
//                            << " @" << solver.num_choicepoints << std::endl;
//                  exit(1);
//                }
//            }
//        }
//        if (not watchers[boolean_watch[~l]].empty()) {
//            num_watchers += watchers[boolean_watch[~l]].size();
//            for (auto cl : watchers[boolean_watch[~l]]) {
//                if (cl->size() < 2) {
//                  std::cout << msg << "(*): error empty clause watching " << l
//                            << " @" << solver.num_choicepoints << std::endl;
//                  exit(1);
//                }
//                
//                if (free_cl_indices.has(cl->id)) {
//                  std::cout << msg << "(*): " << *cl << "'s id is free @"
//                            << solver.num_choicepoints << std::endl;
//                  exit(1);
//                }
//                
//                if (cl->watched(0) != l and cl->watched(1) != l) {
//                  std::cout << msg << "(*): error on clause " << cl->id
//                            << " -- " << *cl << " on " << l
//                            << "'s watch-list but marked as watching "
//                            << cl->watched(0) << " and " << cl->watched(1)
//                            << " @" << solver.num_choicepoints << std::endl;
//                  exit(1);
//                }
//            }
//        }
//    }
//    for (size_t x{0}; x < solver.numeric.size(); ++x) {
//        auto l{lb<T>(x)};
//        auto p{ub<T>(x)};
//
//#ifdef CORRECT_ORDER
//        T prev{Constant::Infinity<T>};
//#else
//        T prev{-Constant::Infinity<T>};
//#endif
//        for (auto &wl : numeric_watch[l]) {
//#ifdef CORRECT_ORDER
//            if (wl.first >= prev) {
//#else
//            if (wl.first <= prev) {
//#endif
//              std::cout << msg << "(*): error numeric watches not sorted @"
//                        << solver.num_choicepoints << std::endl;
//              exit(1);
//            }
//            prev = wl.first;
//            
//            auto &watches{watchers[wl.second]};
//            num_watchers += watches.size();
//            if (not watches.empty()) {
//                for (auto cl : watches) {
//                    if (cl->size() < 2) {
//                      std::cout << msg << "(*): error empty clause watching "
//                                << l << " @" << solver.num_choicepoints
//                                << std::endl;
//                      exit(1);
//                    }
//                    
//                    if (free_cl_indices.has(cl->id)) {
//                      std::cout << msg << "(*): " << *cl << "'s id is free @"
//                                << solver.num_choicepoints << std::endl;
//                      exit(1);
//                    }
//                    
//                    if (not((cl->watched(0).sameVariable(l) and
//                             cl->watched(0).sign() != l.sign()) or
//                            (cl->watched(1).sameVariable(l) and
//                             cl->watched(1).sign() != l.sign()))) {
//                      std::cout << msg << "(*): error on clause " << cl->id
//                                << " -- " << *cl << " on " << l
//                                << "'s watch-list but marked as watching "
//                                << cl->watched(0) << " and " << cl->watched(1)
//                                << " @" << solver.num_choicepoints << std::endl;
//                      exit(1);
//                    }
//                }
//            } else if (wl.first != Constant::Infinity<T>) {
//              std::cout << msg
//                        << "(*): empty watch list guarded by a finite value "
//                        << l << " @" << solver.num_choicepoints << std::endl;
//              exit(1);
//            }
//        }
//#ifdef CORRECT_ORDER
//            prev = Constant::Infinity<T>;
//#else
//            prev = -Constant::Infinity<T>;
//#endif
//
//        for (auto &wl : numeric_watch[p]) {
//#ifdef CORRECT_ORDER
//            if (wl.first >= prev) {
//#else
//            if (wl.first <= prev) {
//#endif
//              std::cout << msg << "(*): error numeric watches not sorted "
//                        << " @" << solver.num_choicepoints << std::endl;
//              exit(1);
//            }
//            prev = wl.first;
//            
//            auto &watches{watchers[wl.second]};
//            num_watchers += watches.size();
//            if (not watches.empty()) {
//                for (auto cl : watches) {
//                    if (cl->size() < 2) {
//                      std::cout << msg << "(*): error empty clause watching "
//                                << l << " @" << solver.num_choicepoints
//                                << std::endl;
//                      exit(1);
//                    }
//                    
//                    if (free_cl_indices.has(cl->id)) {
//                      std::cout << msg << "(*): " << *cl << "'s id is free @"
//                                << solver.num_choicepoints << std::endl;
//                      exit(1);
//                    }
//                    
//                    if (not((cl->watched(0).sameVariable(p) and
//                             cl->watched(0).sign() != p.sign()) or
//                            (cl->watched(1).sameVariable(p) and
//                             cl->watched(1).sign() != p.sign()))) {
//                      std::cout << msg << "(*): error on clause " << cl->id
//                                << " -- " << *cl << " on " << p
//                                << "'s watch-list but marked as watching "
//                                << cl->watched(0) << " and " << cl->watched(1)
//                                << " @" << solver.num_choicepoints << std::endl;
//                      exit(1);
//                    }
//                }
//            } else if (wl.first != Constant::Infinity<T>) {
//              std::cout << msg
//                        << "(*): empty watch list guarded by a finite value "
//                        << l << " @" << solver.num_choicepoints << std::endl;
//              exit(1);
//            }
//        }
//    }
//    
//    if (num_watchers != 2 * size()) {
//        std::cout << msg << ": wrong number of watchers !\n";
//        exit(1);
//    }
//    
//    //    std::cout << msg << ": ok\n";
}
#endif

template <typename T>
std::ostream &operator<<(std::ostream &os, const ClauseBase<T> &x) {
    return x.display(os);
}

}

#endif // _TEMPO_CLAUSEBASE_HPP

