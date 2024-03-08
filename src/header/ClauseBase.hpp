

#ifndef _TEMPO_CLAUSEBASE_HPP
#define _TEMPO_CLAUSEBASE_HPP

#include <vector>
#include <assert.h>

#include "Global.hpp"
#include "Clause.hpp"
#include "Failure.hpp"
#include "SparseSet.hpp"


namespace tempo {

template<typename T> class Scheduler;
template<typename T> class BoundConstraint;

template<class T>
class ClauseBase : Explainer {
    
public:
    
    std::vector<Clause *> base;
    //    std::vector<Clause *> learnt;
    
    ClauseBase(Scheduler<T>&);
    ~ClauseBase() = default;
    
    void resize(const size_t n, const size_t m);
    size_t size() const;
    
    //    void unit_propagate(const var l);
    void set_watcher(const int r, const index_t i, Clause *cl);
    
    

    template <typename iter>
    void add(const iter first, const iter last);
    
    template <typename iter>
    void learn(const iter first, const iter last);
    
    void unit_propagate(const lit l);
    
    void assign(const lit l, Explanation e);
    
    
    std::ostream& display(std::ostream &os) const;
    
    std::ostream& displayClause(std::ostream &os, const Clause* c) const;
    
    
    bool satisfied(const genlit) const;
    bool falsified(const genlit) const;
    
    
    lit newNegLiteral(const lit l);
    lit getReasonLit(const lit l) const;
    
    
    bool sameLit(const lit, const lit) const;
//    lit getVarInClause(const lit l) const;
//    lit getVarInTrigger(const lit l) const;
    
    // explanation stuff
    void xplain(const lit l, const hint h, std::vector<lit> &Cl) ;
    std::ostream &print_reason(std::ostream &, const hint) const;
    int getType() const;
    
    
    std::string prettyLiteral(const genlit l) const;
    
private:
    
    Scheduler<T>& caller;
    
    // clauses, watch struct watch[BOUND_LIT] for bound literals, watch[EDGE_LIT] for edges
    std::vector<std::vector<Clause*>> watch[2];
    
    std::vector<std::vector<lit>> binary;
    
    // store the set of bound constraints appearing in a literal
    std::vector<BoundConstraint<T>> constraint;
    // the number of clauses in which this bound-constraint appears
    std::vector<int> cardinality;
    
    // free indices in the vector constraints (so that we can remove constraints when they are not used)
    SparseSet<int> free_indices;
    
    // for every event-lit, the set of (pointers to) constraints ordered by bound
    std::vector<std::vector<int>> cons_list;
    
    // statistics
    size_t total_size{0};
    int num_units{0};
    //    int num_restarts{0};
    
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

template<typename T>
ClauseBase<T>::ClauseBase(Scheduler<T>& c) : caller(c) {}

template<typename T>
size_t ClauseBase<T>::size() const {
    return base.size(); // + learnt.size();
    //    return var_level.size();
}

template<typename T>
void ClauseBase<T>::resize(const size_t n, const size_t m) {
    //    var_level.resize(n,-1);
    //    activity.resize(n, 0);
    
    binary.resize(2*m);
    watch[EDGE_LIT].resize(2*m);
    
    //    std::cout << 2*n << std::endl;
    
    watch[BOUND_LIT].resize(2*n);
    
    cons_list.resize(2*n);
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
}

template<typename T>
void ClauseBase<T>::set_watcher(const int r, const index_t i, Clause *cl) {
    cl->watcher_index[r] = i;
    
    //    std::cout << "watcher[" << r << "] = " << cl->watcher(r) << std::endl;
    // add it to p's watch list
    lit l{(*cl)[i]};
    
    if(LTYPE(l) == EDGE_LIT)
        watch[EDGE_LIT][FROM_GEN(l)].push_back(cl);
    else
        watch[BOUND_LIT][constraint[FROM_GEN(l)].l].push_back(cl);
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
    
    //    std::cout << l << std::endl;
    //    std::cout << "create a clause literal negation for " << caller.prettyLiteral(l) << std::endl;
    
    if(LTYPE(l) == EDGE_LIT) {
        lit q{EDGE(NOT(FROM_GEN(l)))};
        //        std::cout << " => " << caller.prettyLiteral(q) << std::endl;
        return q;
    }
    
    
    auto c{~(caller.getBound(FROM_GEN(l)))};
    
    
    //    std::cout << " => bound " << c << std::endl;
    
    lit li{NoLit};
    
    auto it{cons_list[c.l].begin()};
    for(; it!=cons_list[c.l].end(); ++it) {
        if(constraint[*it].distance > c.distance)
            continue;
        if(constraint[*it].distance < c.distance)
            break;
        li = *it;
    }
    
    if(li != NoLit) {
        
        //        std::cout << " already used (" << li << ")" << std::endl;
        
        ++cardinality[li];
    } else {
        if(not free_indices.empty()) {
            li = free_indices.front();
            free_indices.remove_front(li);
            constraint[li] = c;
            cardinality[li] = 1;
            
            //            std::cout << " recycle index (" << li << ")" << std::endl;
        } else {
            li = static_cast<lit>(constraint.size());
            constraint.push_back(c);
            cardinality.push_back(1);
            free_indices.reserve(constraint.size());
            
            //            std::cout << " add new index (" << li << ")" << std::endl;
        }
        
        cons_list[c.l].insert(it, li);
    }
    
    //    std::cout << " return " << BOUND(li) << std::endl;
    
    return BOUND(li);
}

//template<typename T>
//template <typename iter>
//void ClauseBase<T>::add_conflict_to(const iter first, const iter last,
//                          std::vector<Clause *> &clbase) {//}, const int id) {
//    int id{static_cast<int>(clbase.size())};
//  if (first == last) {
//
//#ifdef DBG_TRACE
//     if (DBG_TRACE & UNITPROPAGATION) {
//      std::cout << "FAIL on adding the empty clause!\n";
//     }
//#endif
//
//      throw Failure();
//  } else if (first + 1 == last) {
//    assign(NOT(*first), Constant::NoReason);
//  }
//  else if (first + 2 == last and LTYPE(*first) == EDGE_LIT and LTYPE(*(first+1)) == EDGE_LIT) {
//      lit l1{FROM_GEN(*first)};
//      lit l2{FROM_GEN(*(first+1))};
//      binary[l1].push_back(l2);
//      binary[l2].push_back(l1);
//  }
//  else {
//
//    Clause *c = new Clause(id);
//
////      assert(caller.numVariable() > static_cast<size_t>(VAR(*std::max_element(first, last))));
//
//      auto n{static_cast<size_t>(last-first)};
//    size_t not_falsified[2] = {n, n};
//    int k{0};
//    size_t max_level_lit_index{n};
//    int max_level{std::numeric_limits<int>::min()};
//    for (auto l{first}; l != last; ++l) {
//        lit z{newNegLiteral(*l)};
//
//        std::cout << " push " << z << std::endl;
//
//      c->push_back(z);
////        caller.incrementActivity(VAR(*l));
//      if (k < 2 and not_falsified[k] == n) {
//        if (not falsified(c->back())) {
//          not_falsified[k++] = c->size()-1;
////        } else if (caller.level(VAR(FROM_GEN(*l))) > max_level) {
////          max_level = caller.level(VAR(FROM_GEN(*l)));
////          max_level_lit = l;
////        }
//        } else if (caller.decisionLevel(c->back()) > max_level) {
//          max_level = caller.decisionLevel(c->back());
//          max_level_lit_index = c->size()-1;
//        }
//      }
//    }
//
//      std::cout << k << std::endl;
//
////      exit(1);
//
//    if (k == 0) {
//
//#ifdef DBG_TRACE
//     if (DBG_TRACE & UNITPROPAGATION) {
//      std::cout << "FAIL on adding a falsified clause!\n";
//     }
//#endif
//
//        throw Failure(); //Constant::NoReason);
//    }
//
//    for (auto i{0}; i < k; ++i) {
//
//        std::cout << "literal @" << not_falsified[i] << " = " << prettyLiteral((*c)[not_falsified[i]]) << " is not falsified (watcher)\n";
//
//        set_watcher(i, not_falsified[i], c);
//    }
//
//    if (k < 2) {
//
//        std::cout << "literal @" << max_level_lit_index << " = " << prettyLiteral((*c)[max_level_lit_index]) << " is the max-level falsified lit (watcher)\n";
//
//
//      set_watcher(1, max_level_lit_index, c);
//
//        lit l{(*c)[not_falsified[0]]};
//
//        assign(l, {c,l});
//
//    }
//
//
//
//
//    clbase.push_back(c);
//
//    total_size += c->size();
//
//      std::cout << "added ";
//      displayClause(std::cout, c);
//      std::cout << std::endl;
//
//  }
//
////    return max_level;
//}


template<typename T>
template <typename iter>
void ClauseBase<T>::add(const iter first, const iter last) {//}, const int id) {
    int id{static_cast<int>(size())};
    if (first == last) {
        
#ifdef DBG_TRACE
        if (DBG_TRACE & UNITPROPAGATION) {
            std::cout << "FAIL on adding the empty clause!\n";
        }
#endif
        
        throw Failure();
    } else if (first + 1 == last) {
        
//        std::cout << "unit clause " << prettyLiteral(*first) << " b/c cl[" << cl->id << "]: " ;
//        displayClause(std::cout, cl);
//        std::cout << "\n";
        
        assign(*first, Constant::NoReason);
    }
    //  else if (first + 2 == last and LTYPE(*first) == EDGE_LIT and LTYPE(*(first+1)) == EDGE_LIT) {
    //      lit l1{FROM_GEN(*first)};
    //      lit l2{FROM_GEN(*(first+1))};
    //      binary[l1].push_back(l2);
    //      binary[l2].push_back(l1);
    //  }
    else {
        
        Clause *c = new Clause(id);
        
        //      assert(caller.numVariable() > static_cast<size_t>(VAR(*std::max_element(first, last))));
        
        auto n{static_cast<size_t>(last-first)};
        size_t not_falsified[2] = {n, n};
        int k{0};
        size_t max_level_lit_index{n};
        int max_level{std::numeric_limits<int>::min()};
        for (auto l{first}; l != last; ++l) {
            //        lit z{newNegLiteral(*l)};
            lit z{*l};
            
            //        std::cout << " push " << z << std::endl;
            
            c->push_back(z);
            //        caller.incrementActivity(VAR(*l));
            if (k < 2 and not_falsified[k] == n) {
                if (not falsified(c->back())) {
                    not_falsified[k++] = c->size()-1;
                    //        } else if (caller.level(VAR(FROM_GEN(*l))) > max_level) {
                    //          max_level = caller.level(VAR(FROM_GEN(*l)));
                    //          max_level_lit = l;
                    //        }
                } else if (caller.decisionLevel(c->back()) > max_level) {
                    max_level = caller.decisionLevel(c->back());
                    max_level_lit_index = c->size()-1;
                }
            }
        }
        
        //      std::cout << k << std::endl;
        
        //      exit(1);
        
        if (k == 0) {
            
#ifdef DBG_TRACE
            if (DBG_TRACE & UNITPROPAGATION) {
                std::cout << "FAIL on adding a falsified clause!\n";
            }
#endif
            
            throw Failure(); //Constant::NoReason);
        }
        
        for (auto i{0}; i < k; ++i) {
            
            //#ifdef DBG_TRACE
            //     if (DBG_TRACE & UNITPROPAGATION) {
            //        std::cout << "literal @" << not_falsified[i] << " = " << prettyLiteral((*c)[not_falsified[i]]) << " is not falsified (watcher)\n";
            //     }
            //#endif
            
            set_watcher(i, not_falsified[i], c);
        }
        
        if (k < 2) {
            
            //#ifdef DBG_TRACE
            //     if (DBG_TRACE & UNITPROPAGATION) {
            //        std::cout << "literal @" << max_level_lit_index << " = " << prettyLiteral((*c)[max_level_lit_index]) << " is the max-level falsified lit (watcher)\n";
            //     }
            //#endif
            
            set_watcher(1, max_level_lit_index, c);
            
            lit l{(*c)[not_falsified[0]]};
            
            
            std::cout << "branch right " << prettyLiteral(l) << " b/c cl[" << c->id << "]: " ;
            displayClause(std::cout, c);
            std::cout << "\n";
            
            assign(l, {this,id});
            
        }
        
        
        
        
        base.push_back(c);
        
        total_size += c->size();
        
        //      std::cout << "added ";
        //      displayClause(std::cout, c);
        //      std::cout << std::endl;
        
    }
    
    //    return max_level;
}



template<typename T>
template <typename iter>
void ClauseBase<T>::learn(const iter first, const iter last) {
    
    int id{static_cast<int>(size())};
    assert(first != last);
    
    if (first + 1 == last) {
        assign(*first, Constant::NoReason);
    } else {

        Clause *c = new Clause(id);
        
        for (auto l{first}; l != last; ++l) {
            c->push_back(*l);
        }
        
        set_watcher(0, 0, c);
        set_watcher(1, 1, c);
        
        lit l{(*c)[0]};
        
        base.push_back(c);
        
        assign(l, {this,id});
        
        total_size += c->size();
    }
}

template<typename T>
std::ostream& ClauseBase<T>::displayClause(std::ostream &os, const Clause* cl) const {
//    os << "[" << cl->watch_index(0) << "|" << cl->watch_index(1) << "]";
    os << "(" << "[" << prettyLiteral((*cl)[0]) << "]";
    if(0 == cl->watch_index(0) or 0 == cl->watch_index(1))
        os << "*";
    for(size_t i{1}; i<cl->size(); ++i) {
        os << " or [" << prettyLiteral((*cl)[i]) << "]";
        if(i == cl->watch_index(0) or i == cl->watch_index(1))
            os << "*";
    }
    os << ")";
    return os;
}

template<typename T>
void ClauseBase<T>::xplain(const lit l, const hint h, std::vector<lit> &Cl) {
//    assert((*this)[h] == l);
//    for(auto p : *this)
//        if(p != l)
//            Cl.push_back(p);
    
//    std::cout << "explain " << caller.prettyLiteral(l) << " by cl[" << h << "]: " ;
//    displayClause(std::cout, base[h]);
//    std::cout << "\n";
    
    Clause& reason(*(base[h]));
    
    if(l == NoLit) {
        for(auto p : reason)
            Cl.push_back(getReasonLit(p));
    } else {
        
        bool in_clause{false};
        for(auto p : reason)
            //            if(p != l)
            if(sameLit(p,l)) {
                in_clause = true;
            } else {
                Cl.push_back(getReasonLit(p));
            }
        
        assert(in_clause);
    }
//    auto  n{size()};
//    for(size_t i=1; i<n; ++i) {
//        Cl.push_back(this->operator[]((h+i)%n));
//    }
}

template<typename T>
std::ostream &ClauseBase<T>::print_reason(std::ostream &os, const hint h) const {
    return displayClause(os, base[h]);
}

template<typename T>
int ClauseBase<T>::getType() const {
    return CLAUSEEXPL;
}

template<typename T>
std::ostream &ClauseBase<T>::display(std::ostream &os
//                                            , const std::function<int>& f=TODIMACS
                                            ) const {
    os << "base:";
    for(auto cl : base)
        os << " " << *cl;
                                                os << std::endl;
//    os << "learnt:";
//    for(auto cl : learnt)
//        os << " " << *cl;
//    os << "\nbinary:";
//    for(size_t x{0}; x<caller.numVariable(); ++x) {
//        if(not binary[NEG(x)].empty()) {
//            os << " " << POS(x) << " ->";
//            for(auto p : binary[NEG(x)]) {
//                os << " " << p;
//            }
//            //            os << std::endl;
//        }
//        if(not binary[POS(x)].empty()) {
//            os << " " << NEG(x) << " ->";
//            for(auto p : binary[POS(x)]) {
//                os << " " << p;
//            }
//            //            os << std::endl;
//        }
//    }
    for(size_t x{0}; x<caller.numVariable(); ++x) {
        if(not watch[EDGE_LIT][POS(x)].empty()) {
            os << POS(x) << " is watched in ";
            for(auto cl : watch[EDGE_LIT][POS(x)]) {
                os << *cl;
            }
            os << std::endl;
        }
        if(not watch[EDGE_LIT][NEG(x)].empty()) {
            os << NEG(x) << " is watched in ";
            for(auto cl : watch[EDGE_LIT][NEG(x)]) {
                os << *cl;
            }
            os << std::endl;
        }
    }
  return os;
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const tempo::ClauseBase<T> &x) {
  return x.display(os);
}


}

#endif // _TEMPO_CLAUSEBASE_HPP

