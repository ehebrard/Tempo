
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

    void newEdge(const event x, const event y, const T d);

    void newBound(const BoundConstraint<T>& c, Explanation e=Constant::NoReason);
//    task newTask(const T min_dur, const T max_dur);
    
//    T lowerBound(const event x) const;
//    T upperBound(const event x) const;
    
    size_t size() const;
    
    std::ostream &display(std::ostream& os, const bool print_graph=false) const;
    std::ostream &displayConstraint(std::ostream &os, const DistanceConstraint<T>&) const;
    
    BoundSystem<T> bounds;
    
//    lit getLiteral(const size_t idx) const;
    
    
    
    void xplain(const lit, const hint, std::vector<lit> &) override;// {}
    std::ostream &print_reason(std::ostream &os, const hint) const override;// { return os; }
    int getType() const override;// { return CYCLEEXPL; }
    void findExplanationPath(const event x, const event y, std::vector<lit>& path_cl,  const lit time_stamp) ;

    size_t arcCount() const { return core.arcCount(); }

    const auto &getForwardGraph() const;
    const auto &getBackwardGraph() const;

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

    void boundClosure(const event x, const event y, const T d, const lit e_id);

    template <typename G>
    void update(const bool bt, const int s, const G &neighbors,
                const std::vector<T> &shortest_path);
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

template <typename T> const auto &TemporalNetwork<T>::getForwardGraph() const {
  return core.forward();
}

template <typename T> const auto &TemporalNetwork<T>::getBackwardGraph() const {
  return core.backward();
}

template<typename T>
void TemporalNetwork<T>::xplain(const lit l, const hint h, std::vector<lit> &Cl) {
    if(l == NoLit) {
        
//        std::cout << "xplain failure\n";

        event x{EVENT(h)};
        auto s{SIGN(h)};

        do {
            
            auto p{bounds.reason[bounds.getIndex(LIT(x,s))].the_hint};
            
            if(LTYPE(p) == EDGE_LIT) {
                auto el{sched.getEdgeLiteral(FROM_GEN(p))};
                Cl.push_back(EDGE(el));
                x = (s == LOWER ? sched.getEdge(el).to
                                : sched.getEdge(el).from);
            } else {
                x = EVENT(bounds.getConstraint(FROM_GEN(p)).l);
            }
            
        } while (EVENT(h) != x);

    } else {

        if(LTYPE(h) == EDGE_LIT) {
            
            
            
            assert(FROM_GEN(h) <= bounds.getStamp(FROM_GEN(l)));
            
            // h is the responsible edge
            auto el{sched.getEdgeLiteral(FROM_GEN(h))};
            auto ec{sched.getEdge(el)};
            auto cx{bounds.getConstraint(FROM_GEN(l))};
            
//            std::cout << "\nxplain " << sched.prettyLiteral(l) << " b/c " ;
//            displayConstraint(std::cout, ec);
//            std::cout << std::endl;
            
            //            auto bl{FROM_GEN(l)};
            if(SIGN(cx.l) == LOWER) {
                // evt(cx.l) = ec.from -> ec.to
                assert(EVENT(cx.l) == ec.from);
                
                
                BoundConstraint<T> eb{LIT(ec.to, LOWER), cx.distance-ec.distance};
                
//                std::cout << " AND ";
//                bounds.displayBound(std::cout,eb);
//                std::cout << std::endl;
                
                auto r{BOUND(bounds.getImplicant(eb))};
                
//                std::cout << " ==> " << sched.prettyLiteral(r) << std::endl;
                
                assert(r < l);
                
                Cl.push_back(r);
            } else {
                // evt(cx.l) = ec.to <- ec.from

                assert(EVENT(cx.l) == ec.to);
                
                BoundConstraint<T> eb{LIT(ec.from, UPPER), cx.distance-ec.distance};
                
//                std::cout << " AND ";
//                bounds.displayBound(std::cout,eb);
//                std::cout << std::endl;
                
                auto r{BOUND(bounds.getImplicant(eb))};
                
//                std::cout << " ==> " << sched.prettyLiteral(r) << std::endl;
                
                assert(r < l);
                
                Cl.push_back(r);
            }
            Cl.push_back(EDGE(el));
        } else {
            
//            std::cout << "xplain bound\n";
            
            Cl.push_back(h);
        }
    }
}

template<typename T>
std::ostream &TemporalNetwork<T>::print_reason(std::ostream &os, const hint h) const {
    if(h<0)
        os << "negative cycle via " << bounds.getLabel(h) ;
    else {
        if(LTYPE(h) == EDGE_LIT) {
            // h is the responsible edge
            auto el{sched.getEdgeLiteral(FROM_GEN(h))};
            auto ec{sched.getEdge(el)};
            os << "shortest path via " ;
            displayConstraint(os, ec);
        } else {
          os << "shortest path via " << sched.prettyLiteral(h)
             << " (and ground edge)";
        }
    }
//    os << "precedences";
    return os;
}

template<typename T>
int TemporalNetwork<T>::getType() const { return CYCLEEXPL; }

template <typename T>
void TemporalNetwork<T>::newEdge(const event x, const event y, const T d) {

  //    assert(e.the_hint == (static_cast<int>(sched.numEdgeLiteral())-1));
  auto e_id{static_cast<lit>(sched.numEdgeLiteral()) - 1};

  core.emplace_edge(x, y, d, e_id); // sched.numEdgeLiteral());
  boundClosure(x, y, d, EDGE(e_id));
}

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

template <typename T>
void TemporalNetwork<T>::newBound(const BoundConstraint<T>& c, Explanation e)
{
    event x{EVENT(c.l)};
    if(SIGN(c.l) == LOWER) {
      if (bounds.set(LOWER, x, c.distance, e)) {
        update(LOWER, x, core.backward(), bounds.lower());
      }
    } else {
      if (bounds.set(UPPER, x, c.distance, e)) {
        update(UPPER, x, core, bounds.upper());
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
                    std::cout << " * shorter path to " << bounds.getLabel(v) << " (" << length[u]+1 << "/" << size() << ") "<< std::endl ;
                }
#endif

                path[v] = u;
                elit[v] = s;
                length[v] = length[u]+1;

                if (length[v] > size()) {
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
                std::cout << " does not push " << bounds.getLabel(v) << std::endl ;
            }
#endif

            }
#ifdef DBG_BELLMAN_EXPL
            else if(DBG_BELLMAN_EXPL) {
                std::cout << " ignore " << bounds.getLabel(v) << std::endl ;
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
        std::cout << bounds.getLabel(z) ;
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
                
                std::cout << " <= " << bounds.getLabel(z) ;
            }
#endif
            
            path_cl.push_back(EDGE(l));
            assert(sched.getEdge(l).from == z);
        } 
#ifdef DBG_BELLMAN_EXPL
        else if(DBG_BELLMAN_EXPL) {
            std::cout << " <- " << bounds.getLabel(z);
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
}

template <typename T>
template <typename G>
void TemporalNetwork<T>::update(const bool bt, const event s,
                                const G &neighbors,
                                const std::vector<T> &shortest_path) {

  //                                    int max_iter{1000};

  changed.clear();
  changed.add(s);

#ifdef DBG_BELLMAN
    if(DBG_BELLMAN) {
        core.display(std::cout, [this](const event e) {return bounds.getLabel(e);}, [this](const StampedLabeledEdge<T,lit>& e) {return bounds.getLabel(static_cast<int>(e)) + " (" + std::to_string(e.label()) + ")";});
        std::cout << "\nstart explore from " << bounds.getLabel(s) << (bt==LOWER ? " (backward)" : " (forward)") << std::endl ;
    }
#endif
    
    while (not changed.empty()) {
    
        auto u{changed.front()};
        changed.pop_front();
        
#ifdef DBG_BELLMAN
        if(DBG_BELLMAN) {
            std::cout << "pop " << bounds.getLabel(u) << " q=(";
            for(auto evt : changed)
                std::cout << " " << bounds.getLabel(evt);
            std::cout << " )" << std::endl ;
        }
#endif
        
        for (auto edge : neighbors[u]) {
            int v{edge};
            auto w{edge.label()};
            
            if (shortest_path[u] + w < shortest_path[v]) {
                    
#ifdef DBG_BELLMAN
                if(DBG_BELLMAN) {
                    std::cout << " * shorter path " << (bt==LOWER ? "from " : "to ") << bounds.getLabel(v) << std::endl ;
                }
#endif
                    
                    if (v == s) {
                        
                        
                        
#ifdef DBG_TRACE
//                        if(DBG_BELLMAN) {
                            if (DBG_TRACE & PROPAGATION) {
                                std::cout << "FAIL on negative cycle!";
                                event evt{v};
                                std::cout << " " << bounds.getLabel(evt) ;
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
//                        }
#endif
                        
                        
                        throw Failure({this,LIT(s,bt)});//Constant::NoReason);
                    }
                    bounds.set(
                        bt, v, shortest_path[u] + w,
                        {this, (edge.stamp() >= 0
                                    ? EDGE(edge.stamp())
                                    : BOUND(bounds.getIndex(LIT(u, bt))))});

                    if (not changed.has(v))
                        changed.add(v);
                }
#ifdef DBG_BELLMAN
            else if(DBG_BELLMAN) {
                std::cout << " ignore " << bounds.getLabel(v) << std::endl ;
            }
#endif
        }
    }
}

template <typename T>
std::ostream &TemporalNetwork<T>::display(std::ostream& os, const bool print_graph) const {
    os << bounds << std::endl;
//    if(print_graph) core.display(os, [](const event e) {return prettyEvent(e);}, [](const StampedLabeledEdge<T,lit>& e) {return prettyEvent(static_cast<int>(e));});
    
    if(print_graph) core.display(os, [this](const event e) {return bounds.getLabel(e);}, [this](const StampedLabeledEdge<T,lit>& e) {return bounds.getLabel(static_cast<int>(e));});
    return os;
}

template <typename T>
std::ostream &TemporalNetwork<T>::displayConstraint(std::ostream &os, const DistanceConstraint<T>& c) const {
    os << bounds.getLabel(c.to) << " - " << bounds.getLabel(c.from)
       << " <= " << c.distance;
    return os;
}


template <typename T>
std::ostream &operator<<(std::ostream &os, const TemporalNetwork<T> &x) {
  return x.display(os);
}

} // namespace tempo

#endif
