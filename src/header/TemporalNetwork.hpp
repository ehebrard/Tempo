
#ifndef _TEMPO_TEMPORALNETWORK_HPP
#define _TEMPO_TEMPORALNETWORK_HPP

#include <iomanip>
#include <iostream>
#include <vector>

#include "Global.hpp"
#include "BoundSystem.hpp"
#include "Explanation.hpp"
#include "DirectedGraph.hpp"

namespace tempo {


template <typename T> std::string pretty(T d) {
  if (not finite(d))
    return std::string("oo");
  return std::to_string(d);
}





template<typename T> class Scheduler;

template <typename T>
class TemporalNetwork : public Explainer {
    
public:
    
    TemporalNetwork(Scheduler<T>& s); //, BacktrackEnvironment *e=ReversibleObject::env);
    
    ~TemporalNetwork() = default;
    
    void resize(const size_t n);
    
    void newEdge(const event x, const event y, const T d
#ifdef RECOVER
                 , Explanation e=Constant::NoReason
#endif
    );

    
    void newBound(const BoundConstraint<T>& c, Explanation e=Constant::NoReason);
//    task newTask(const T min_dur, const T max_dur);
    
//    T lowerBound(const event x) const;
//    T upperBound(const event x) const;
    
    size_t size() const;
    
    std::ostream &display(std::ostream& os, const bool print_graph=false) const;
    
    BoundSystem<T> bounds;
    
//    lit getLiteral(const size_t idx) const;
    
    
    
    void xplain(const lit, const hint, std::vector<lit> &) override;// {}
    std::ostream &print_reason(std::ostream &os, const hint) const override;// { return os; }
    int getType() const override;// { return CYCLEEXPL; }
    void findExplanationPath(const event x, const event y, std::vector<lit>& path_cl,  const lit time_stamp) ;

    size_t arcCount() const { return core.arcCount(); }

  private:
    Scheduler<T>& sched;
    
    DirectedGraph<StampedLabeledEdge<T,lit>> core;
   
    SparseSet<> changed;

    ////public:
    std::vector<event> path;
    std::vector<lit> elit;
    std::vector<size_t> length;
    
//    private:
    std::vector<T> distance;
   
#ifdef RECOVER
    void boundClosure(const event x, const event y, const T d, Explanation e);
#else
    void boundClosure(const event x, const event y, const T d, const lit e_id);
#endif
    
    template <typename G>
    void update(const bool bt, const int s,
                                    const G &neighbors,
                                     const std::vector<T>& shortest_path
#ifdef RECOVER
                , Explanation e
#endif
                );
    
};


template <typename T>
TemporalNetwork<T>::TemporalNetwork(Scheduler<T>& s)
//, BacktrackEnvironment *e)
: bounds(s,0), sched(s), core(&(s.getEnv())) {
//    resize(2); // origin and end events are always created
}

template <typename T>
size_t TemporalNetwork<T>::size() const {
    return bounds.size();
}

template <typename T>
void TemporalNetwork<T>::resize(const size_t n) {
    if(n != size()) {
        bounds.resize(n);
        core.resize(n);
        changed.reserve(n);
        path.resize(n, NOEVENT);
        elit.resize(n, NoLit);
        length.resize(n, 0);
        distance.resize(n, INFTY);
    }
}


template<typename T>
void TemporalNetwork<T>::xplain(const lit l, const hint h, std::vector<lit> &Cl) {
    if(l == NoLit) {
//                std::cout << "negative cycle: " << prettyEvent(h) ;
        event x{EVENT(h)};
        auto s{SIGN(h)};
        
//        do {
//            
//            if(elit[x] >= 0) {
//                Cl.push_back(EDGE(sched.getEdgeLiteral(elit[x])));
//            }
//        
//            x = path[x];
//        
//        } while (EVENT(h) != x);
        
        do {
            
            auto p{bounds.reason[bounds.getIndex(LIT(x,s))].the_hint};
            
            if(LTYPE(p) == EDGE_LIT) {
                auto el{sched.getEdgeLiteral(FROM_GEN(p))};
                Cl.push_back(EDGE(el));
                x = (s==LOWER ? x = sched.getEdge(el).to : x = sched.getEdge(el).from);
            } else {
                x = EVENT(bounds.getConstraint(FROM_GEN(p)).l);
            }
            
        } while (EVENT(h) != x);
            
//            if(elit[x] >= 0) {
//                
//                assert(LTYPE(bounds.reason[bounds.getIndex(LIT(x,s))].the_hint) == EDGE_LIT);
//                
////                Cl.push_back(EDGE(sched.getEdgeLiteral(elit[x])));
////                                std::cout << " <= " ;
////                std::cout << elit[x] << " / " << FROM_GEN(bounds.reason[bounds.getIndex(LIT(x,s))].the_hint) << std::endl;
//                assert(elit[x] == FROM_GEN(bounds.reason[bounds.getIndex(LIT(x,s))].the_hint));
//                
//                auto el{sched.getEdgeLiteral(FROM_GEN(bounds.reason[bounds.getIndex(LIT(x,s))].the_hint))};
//                
//                Cl.push_back(EDGE(el));
//                
//                if(s==LOWER) {
//                    assert(sched.getEdge(el).to == path[x]);
//                    x = sched.getEdge(el).to;
//                } else {
//                    assert(sched.getEdge(el).from == path[x]);
//                    x = sched.getEdge(el).from;
//                }
//                
//                
//            } else {
//                
//                assert(LTYPE(bounds.reason[bounds.getIndex(LIT(x,s))].the_hint) == BOUND_LIT);
//                
////                bounds.getIndex(LIT(u,bt))
//                assert(path[x] == EVENT(bounds.getConstraint(FROM_GEN(bounds.reason[bounds.getIndex(LIT(x,s))].the_hint)).l));
//                
//                x = EVENT(bounds.getConstraint(FROM_GEN(bounds.reason[bounds.getIndex(LIT(x,s))].the_hint)).l);
//            }
////                        else {
////                            std::cout << " <- " ;
////                        }
////            x = path[x];
////            x = EVENT(bounds.getConstraint(FROM_GEN(bounds.reason[bounds.getIndex(LIT(x,s))].the_hint)).l);
////            std::cout << prettyEvent(x) ;
//            
//        } while (EVENT(h) != x);
////                std::cout << std::endl;
    } else {
        
        
#ifdef RECOVER
        
        assert(LTYPE(l) == BOUND_LIT);
        
        Cl.push_back(BOUND(h));
        
        auto cx{bounds.getConstraint(FROM_GEN(l))};
        
        auto cy{bounds.getConstraint(h)};
        
        if(SIGN(FROM_GEN(l)) == LOWER) {
            
            std::cout << "find path from " << prettyEvent(EVENT(cx.l)) << " to " << prettyEvent(EVENT(cy.l)) << std::endl;
            
            findExplanationPath(EVENT(cx.l), EVENT(cy.l), Cl, static_cast<lit>(bounds.getStamp(FROM_GEN(l))));
        } else {
            
            std::cout << "find path from " << prettyEvent(EVENT(cy.l)) << " to " << prettyEvent(EVENT(cx.l)) << std::endl;
            
            findExplanationPath(EVENT(cy.l), EVENT(cx.l), Cl, static_cast<lit>(bounds.getStamp(FROM_GEN(l))));
        }
#else
        if(LTYPE(h) == EDGE_LIT) {
            
            assert(FROM_GEN(h) <= bounds.getStamp(FROM_GEN(l)));
            
            // h is the responsible edge
            auto el{sched.getEdgeLiteral(FROM_GEN(h))};
            auto ec{sched.getEdge(el)};
            auto cx{bounds.getConstraint(FROM_GEN(l))};
            
            //            auto bl{FROM_GEN(l)};
            if(SIGN(cx.l) == LOWER) {
                // evt(cx.l) = ec.from -> ec.to
                assert(EVENT(cx.l) == ec.from);
                
                auto r{BOUND(bounds.getImplicant({LIT(ec.to, LOWER), cx.distance+ec.distance}))};
                
                assert(r < l);
                
                Cl.push_back(r);
            } else {
                // evt(cx.l) = ec.to <- ec.from
                assert(EVENT(cx.l) == ec.to);
                
                auto r{BOUND(bounds.getImplicant({LIT(ec.from, UPPER), cx.distance+ec.distance}))};
                
                assert(r < l);
                
                Cl.push_back(r);
            }
            Cl.push_back(EDGE(el));
        } else {
            Cl.push_back(h);
        }
#endif
        
    }
}

template<typename T>
std::ostream &TemporalNetwork<T>::print_reason(std::ostream &os, const hint h) const {
    if(h<0)
        os << "negative cycle via " << prettyEvent(h) ;
    else {
        if(LTYPE(h) == EDGE_LIT) {
            // h is the responsible edge
            auto el{sched.getEdgeLiteral(FROM_GEN(h))};
            auto ec{sched.getEdge(el)};
            os << "transitivity via " << ec ;
        } else {
            os << "transitivity via " << sched.prettyLiteral(h) << " (and ground edge)" ;
        }
    }
//    os << "precedences";
    return os;
}

template<typename T>
int TemporalNetwork<T>::getType() const { return CYCLEEXPL; }

template <typename T>
void TemporalNetwork<T>::newEdge(const event x, const event y, const T d
#ifdef RECOVER
                                 , Explanation e
#endif
                                 )
{
    
//    assert(e.the_hint == (static_cast<int>(sched.numEdgeLiteral())-1));
    auto e_id{static_cast<lit>(sched.numEdgeLiteral())-1};
    
    core.emplace_edge(x,y,d,e_id); //sched.numEdgeLiteral());
    
#ifdef RECOVER
    boundClosure(x,y,d,e);
#else
    boundClosure(x,y,d,EDGE(e_id));
#endif
}

#ifdef RECOVER
template <typename T>
void TemporalNetwork<T>::boundClosure(const event x, const event y, const T d, Explanation e) {
    // closure w.r.t. 0 (0 -> x -(d)-> y -> 0)
    
    
    // reduce the lower bound of x and precessor
    if(bounds.set(LOWER, x, bounds.lower()[y] + d, e)) {
        update(LOWER, x, core.backward(), bounds.lower(), e);
    }

    std::cout << 66 << std::endl;
    // increase the upper bound of y and successors
    if(bounds.set(UPPER, y, bounds.upper()[x] + d, e)) {
        update(UPPER, y, core, bounds.upper(), e);
    }
    
}
#else
template <typename T>
void TemporalNetwork<T>::boundClosure(const event x, const event y, const T d, const lit e_id) {
    // closure w.r.t. 0 (0 -> x -(d)-> y -> 0)

    Explanation e{this, e_id};
    if(e_id < 0)
        e = Constant::NoReason;

    // reduce the lower bound of x and precessor
    if(bounds.set(LOWER, x, bounds.lower()[y] + d, e)) {
        update(LOWER, x, core.backward(), bounds.lower());
    }

    // increase the upper bound of y and successors
    if(bounds.set(UPPER, y, bounds.upper()[x] + d, e)) {
        update(UPPER, y, core, bounds.upper());
    }
    
}
#endif

template <typename T>
void TemporalNetwork<T>::newBound(const BoundConstraint<T>& c, Explanation e)
{
    event x{EVENT(c.l)};
    if(SIGN(c.l) == LOWER) {
        if(bounds.set(LOWER, x, c.distance, e)) {
#ifdef RECOVER
            update(LOWER, x, core.backward(), bounds.lower()), {this,static_cast<hint>(bounds.numLiteral())-1});
#else
            update(LOWER, x, core.backward(), bounds.lower());
#endif
        }
    } else {
        if(bounds.set(UPPER, x, c.distance, e)) {
#ifdef RECOVER
            update(UPPER, x, core, bounds.upper(), {this,static_cast<hint>(bounds.numLiteral())-1});
#else
            update(UPPER, x, core, bounds.upper());
#endif
        }
    }
}


// find a an explanation path from x to y (before lit 'time_stamp') and push it onto Cl
template <typename T>
void TemporalNetwork<T>::findExplanationPath(const event x, const event y, std::vector<lit>& path_cl,  const lit time_stamp) {
    
    changed.clear();
    changed.add(x);
    distance[x] = 0;
    
#ifdef DBG_BELLMAN_EXPL
    SparseSet<event> visited(size());
    int iter_max{10000};
    if(DBG_BELLMAN_EXPL) {
        core.display(std::cout, [](const event e) {return prettyEvent(e);}, [](const StampedLabeledEdge<T,lit>& e) {return prettyEvent(static_cast<int>(e))+"|"+std::to_string(e.label());});
        std::cout << "\nstart explore from " << prettyEvent(x) << std::endl ;
    }
#endif
    
    auto critical_path_root{y};
    
//
    while (not changed.empty()) {
    
        
        auto u{changed.front()};
        changed.pop_front();
        
#ifdef DBG_BELLMAN_EXPL
        if(iter_max-- == 0)
                    exit(1);
        
        if(DBG_BELLMAN_EXPL) {
            std::cout << "pop " << prettyEvent(u) << " q=(";
            for(auto evt : changed)
                std::cout << " " << prettyEvent(evt);
            std::cout << " )" << std::endl ;
        }
#endif
        
        assert(u != ORIGIN);
 
        for (auto edge : core[u]) {
            
            
            int v{edge};
            auto w{edge.label()};
            auto s{edge.stamp()};
            if(s > time_stamp) {
                continue;
            }
            
            if (distance[u] + w < distance[v]) {
                
#ifdef DBG_BELLMAN_EXPL
                if(DBG_BELLMAN_EXPL) {
                    std::cout << " * shorter path to " << prettyEvent(v) << " (" << length[u]+1 << "/" << size() << ") "<< std::endl ;
                }
#endif

                path[v] = u;
                elit[v] = s;
                length[v] = length[u]+1;
                
                if(length[v] > size()) {
//                    std::cout << "negative cycle detected!\n";
//                    exit(1);
                    critical_path_root = v;
                    changed.setStart(changed.end_idx());
                    break;
                }
                
                assert (v != x);
                
                distance[v] = distance[u] + w;
                
                if (v != ORIGIN and not changed.has(v))
                    changed.add(v);
#ifdef DBG_BELLMAN_EXPL
            else if(DBG_BELLMAN_EXPL) {
                std::cout << " does not push " << prettyEvent(v) << std::endl ;
            }
#endif
                

                
//#ifdef DBG_BELLMAN_EXPL
//                if(DBG_BELLMAN_EXPL) {
//                    auto z{v};
//                    visited.clear();
//                    visited.add(z);
////                    std::cout << "d[" << prettyEvent(v) << "]=" << distance[v] << " / d[" << prettyEvent(181) << "]=" << distance[181] << std::endl;
//                    std::cout << prettyEvent(z) ;
//                    while(z != x) {
//                        z = path[z];
//                        if(visited.has(z)) {
//                            //                        break;
//                            std::cout << std::endl;
//                            exit(1);
//                        }
//                        std::cout << " <-|- " << prettyEvent(z);
//                        visited.add(z);
//                    };
//                    std::cout << std::endl;
//                }
//#endif

                
            }
#ifdef DBG_BELLMAN_EXPL
            else if(DBG_BELLMAN_EXPL) {
                std::cout << " ignore " << prettyEvent(v) << std::endl ;
            }
#endif
        }
        
//        changed.pop_front();
    }
    
    for(auto i{changed.fbegin()}; i!=changed.fend(); ++i) {
        distance[*i] = INFTY;
        length[*i] = 0;
        path[*i] = NOEVENT;
        elit[*i] = NoLit;
    }
    
    auto z{critical_path_root};
    
    changed.add(z);
    
#ifdef DBG_BELLMAN_EXPL
    if(DBG_BELLMAN_EXPL) {
        std::cout << prettyEvent(z) ;
    }
                        visited.clear();
                        visited.add(z);
#endif
    
    while(z != x) {
        auto p{elit[z]};
//        std::cout << " [" << p << "/" << prettyEvent(path[z]) << "] ";
        
        z = path[z];
        if(p >= 0) {
            auto l{sched.getEdgeLiteral(p)};
            
#ifdef DBG_BELLMAN_EXPL
            if(DBG_BELLMAN_EXPL) {
                
                std::cout << " <= " << prettyEvent(z) ;
            }
#endif
            
            path_cl.push_back(EDGE(l));
            assert(sched.getEdge(l).from == z);
        } 
#ifdef DBG_BELLMAN_EXPL
        else if(DBG_BELLMAN_EXPL) {
            std::cout << " <- " << prettyEvent(z);
            if(visited.has(z)) {
                                std::cout << std::endl;
                                exit(1);
                            }
            visited.add(z);
        }
#endif
        
        if(changed.has(z))
            break;
        changed.add(z);
        
//        if(z == critical_path_root)
//            break;
    }
    
#ifdef DBG_BELLMAN_EXPL
    if(DBG_BELLMAN_EXPL) {
        std::cout << std::endl;
    }
#endif
    //                        exit(1);
    
    

  
}

template <typename T>
template <typename G>
void TemporalNetwork<T>::update(const bool bt,
                                const event s,
                                const G &neighbors,
                                const std::vector<T>& shortest_path
#ifdef RECOVER
                                , Explanation e
#endif
                                ) {
    
//                                    int max_iter{1000};
                                    
    changed.clear();
    changed.add(s);
    
#ifdef DBG_BELLMAN
    core.display(std::cout, [](const event e) {return prettyEvent(e);}, [](const StampedLabeledEdge<T,lit>& e) {return prettyEvent(static_cast<int>(e));});
    std::cout << "\nstart explore from " << prettyEvent(s) << (bt==LOWER ? " (backward)" : " (forward)") << std::endl ;
#endif
    
    while (not changed.empty()) {
        
//        if(--max_iter < 0)
//            exit(1);
        
        
        auto u{changed.front()};
        changed.pop_front();
        
#ifdef DBG_BELLMAN
        std::cout << "pop " << prettyEvent(u) << " q=(";
        for(auto evt : changed)
            std::cout << " " << prettyEvent(evt);
        std::cout << " )" << std::endl ;
#endif
        
        for (auto edge : neighbors[u]) {
            int v{edge};
            auto w{edge.label()};
            
            if (shortest_path[u] + w < shortest_path[v]) {
                    
#ifdef DBG_BELLMAN
                    std::cout << " * shorter path " << (bt==LOWER ? "from " : "to ") << prettyEvent(v) << std::endl ;
#endif
                    
//                    path[v] = u;
//                elit[v] = edge.stamp();
                    
                    if (v == s) {
                        
                        
                        
#ifdef DBG_TRACE
                        if (DBG_TRACE & PROPAGATION) {
                            std::cout << "FAIL on negative cycle!";
                            event evt{v};
                            std::cout << " " << prettyEvent(evt) ;
                            //                            do {
                            //                                if(bt==LOWER)
                            //                                    std::cout << "
                            //                                    -> ";
                            //                                else
                            //                                    std::cout << "
                            //                                    <- ";
                            //                                evt = path[evt];
                            //                                std::cout <<
                            //                                prettyEvent(evt) ;
                            //                            } while (evt != v);
                            std::cout << std::endl;
                        }
#endif
                        
                        
                        throw Failure({this,LIT(s,bt)});//Constant::NoReason);
                    }
#ifdef RECOVER
                    bounds.set(bt, v, shortest_path[u] + w, e);
#else
                    bounds.set(bt, v, shortest_path[u] + w, {this, (edge.stamp() >= 0 ? EDGE(edge.stamp()) : BOUND(bounds.getIndex(LIT(u,bt))))});
#endif
                    
                    if (not changed.has(v))
                        changed.add(v);
                }
#ifdef DBG_BELLMAN
            else
                std::cout << " ignore " << prettyEvent(v) << std::endl ;
#endif
        }
    }
}



template <typename T>
std::ostream &TemporalNetwork<T>::display(std::ostream& os, const bool print_graph) const {
    os << bounds << std::endl;
    if(print_graph) core.display(os, [](const event e) {return prettyEvent(e);}, [](const StampedLabeledEdge<T,lit>& e) {return prettyEvent(static_cast<int>(e));});
    return os;
}


template <typename T>
std::ostream &operator<<(std::ostream &os, const TemporalNetwork<T> &x) {
  return x.display(os);
}

} // namespace tempo

#endif
