/************************************************
 * Tempo FullTransitivity.hpp
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

#ifndef TEMPO_FULLTRANSITIVITY_HPP
#define TEMPO_FULLTRANSITIVITY_HPP

#include <cassert>
#include <map>
#include <vector>

#include "Explanation.hpp"
#include "Global.hpp"
#include "ReversibleObject.hpp"
#include "constraints/Constraint.hpp"
#include "util/SparseSet.hpp"

namespace tempo {

template<typename T>
class Solver;


template<typename T>
struct PathExplanation {
    Literal<T> literal{Solver<T>::Contradiction};
    index_t prefix{Constant::NoIndex};
    index_t suffix{Constant::NoIndex};
};


template <typename T> class FullTransitivity : public Constraint<T>
//, public ReversibleObject<T> 
{
private:
  Solver<T> &m_solver;

    // all the edges, in chronological order
    std::vector<DistanceConstraint<T>> edges;
    // the reason for each edge[: a literal, and the leading and following sub-paths]
    std::vector<PathExplanation<T>> reason;
    // for each pair x,y, the index of the last literal
    std::vector<std::vector<index_t>> _index_;
    // for each edge, the index of the previous edge
    std::vector<index_t> previous;
    
    
    // the distance from
    std::vector<std::vector<T>> distance_from;
    // the distance to (transpose of 'distance_from')
    std::vector<std::vector<T>> distance_to;
    
    // for each pair x,y, the corresponding disjunct literal, or Solver::Contradiction<T> if there is none
    std::vector<std::vector<Literal<T>>> literal;
    
    // util for Bellman-Ford
    SparseSet<> changed;

public:
  FullTransitivity(Solver<T> &solver);
  virtual ~FullTransitivity();

    // add a set of disjuncts
    template <typename Iter>
    void addResource(const Iter beg_disjunct, const Iter end_disjunct);

  bool notify(const Literal<T>, const int rank) override;
  void post(const int idx) override;
  void propagate() override;
    
    T distance(const int x, const int y) const;
    bool addEdge(const int x, const int y, const T d, const PathExplanation<T> r=NoReason);
    void undo();
    
    void propagateForward(const DistanceConstraint<T> e, const Literal<T> l=Solver<T>::Contradiction);
    void propagateBackward(const DistanceConstraint<T> e, const Literal<T> l=Solver<T>::Contradiction);
    
    template <typename G>
    void update(const int x, const int y, const G &neighbors, std::vector<T>& distance);

  void xplain(const Literal<T> l, const hint h, std::vector<Literal<T>> &Cl) override;
//  int getType() const override;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;
    
    void printMatrix() const;
    
    
    static PathExplanation<T> NoReason;

};

template <typename T>
PathExplanation<T> FullTransitivity<T>::NoReason = PathExplanation<T>();


template <typename T>
T FullTransitivity<T>::distance(const int x, const int y) const {
    return edges[_index_[x][y]].distance;
}

template <typename T>
bool FullTransitivity<T>::addEdge(const int x, const int y, const T d, const PathExplanation<T> r) {
    
    if(distance(x,y) <= d)
        return false;
    
#ifdef DBG_FTRANS
    if (DBG_FTRANS) {
        std::cout << "  - x" << x << " -> x" << y << " (" << distance(x,y) << "/" << d <<  ")\n";
    }
#endif
    
    auto i{static_cast<index_t>(edges.size())};
    
    distance_from[x][y] = distance_to[y][x] = d;
    
    previous.push_back(_index_[x][y]);
    
    _index_[x][y] = i;
    edges.emplace_back(x,y,d);
    reason.push_back(r);
    

    if(literal[y][x] != Solver<T>::Contradiction and m_solver.boolean.getEdge(literal[y][x]).distance + distance_from[x][y] < 0) {
        
#ifdef DBG_FTRANS
    if (DBG_FTRANS) {
        std::cout << "  --> infer literal " << m_solver.pretty(literal[x][y]) << "\n";
    }
#endif
        
        m_solver.set(literal[x][y], {this, static_cast<hint>(i)});
    }
    
    return true;
}


template <typename T>
void FullTransitivity<T>::undo() {
    auto r{reason.back()};
    auto e{edges.back()};
    auto i{previous.back()};

    assert(i == 0 or edges[i].from = e.from and edges[i].to = e.to);
    
    distance_from[e.from][e.to] = distance_to[e.to][e.from] = edges[i].distance;

    reason.pop_back();
    edges.pop_back();
    previous.pop_back();
}



template <typename T>
FullTransitivity<T>::FullTransitivity(Solver<T> &solver)
    : m_solver(solver) {

  Constraint<T>::priority = Priority::Low;
        
        auto n{m_solver.numeric.size()};
        
        _index_.resize(n);
        distance_from.resize(n);
        distance_to.resize(n);
        literal.resize(n);
        
        changed.reserve(n);
        
        edges.push_back({Constant::NoVar, Constant::NoVar, Constant::Infinity<T>});
        reason.push_back(NoReason);
        for(auto &row : _index_)
            while(row.size() < n)
                row.emplace_back(0);
        
        for(auto &row : distance_from)
            while(row.size() < n)
                row.emplace_back(Constant::Infinity<T>);
        for(auto &row : distance_to)
            while(row.size() < n)
                row.emplace_back(Constant::Infinity<T>);
        
        
        for(auto &row : literal)
            row.resize(n, Solver<T>::Contradiction);
        
        for(size_t x{0}; x<_index_.size(); ++x) {
            _index_[x][x] = Constant::NoIndex;
            distance_from[x][x] = 0;
            distance_to[x][x] = 0;
        }
        
//        std::cout << "\nforward:\n";
//        for(auto x : m_solver.core) {
//        std::cout << x;
//            for(auto e : m_solver.core[x]) {
//                int y{e};
//                std::cout << " " << y;
//            }
//            std::cout << std::endl;
//        }
//        
//        std::cout << "\nbackward:\n";
//        for(auto x : m_solver.core) {
//        std::cout << x;
//            for(auto e : m_solver.core.backward()[x]) {
//                int y{e};
//                std::cout << " " << y;
//            }
//            std::cout << std::endl;
//        }
        printMatrix();
        
        for(auto x : m_solver.core) {
            
            
            
            if(x != 0) {
                if(m_solver.numeric.upper(x) != Constant::Infinity<T>) {
                    
                    std::cout << "\n**upper bound of " << x << std::endl;
                    
                    if(addEdge(0,x,m_solver.numeric.upper(x))) {
                        propagateForward(edges.back());
                        propagateBackward(edges.back());
                    }
//                    propagateForward({0,x,m_solver.numeric.upper(x)});
//                    propagateBackward({0,x,m_solver.numeric.upper(x)});
                    
//                    update(0,x,m_solver.core,distance_from[0]);
//                    update(x,0,m_solver.core.backward(),distance_to[x]);
                    
                    printMatrix();
                }
                
                if(m_solver.numeric.lower(x) != -Constant::Infinity<T>) {
                    
                    std::cout << "\n**lower bound of " << x << std::endl;
                    
                    if(addEdge(x,0,-m_solver.numeric.lower(x))) {
                        propagateForward(edges.back());
                        propagateBackward(edges.back());
                    }
//                    propagateForward({x,0,-m_solver.numeric.lower(x)});
//                    propagateBackward({x,0,-m_solver.numeric.lower(x)});
                    
//                    update(x,0,m_solver.core,distance_from[x]);
//                    update(0,x,m_solver.core.backward(),distance_to[0]);
                    
                    printMatrix();
                }
            }
            
            for(auto e : m_solver.core[x]) {
                
                int y{e};
                
                std::cout << "\n**edge " << x << " -> " << y << std::endl;
                
                if(addEdge(x,y,e.label())) {
                    propagateForward(edges.back());
                    propagateBackward(edges.back());
                }
                
                
//                
//                update(x,y,m_solver.core,distance_from[x]);
//                
//                
//                for(auto b : changed) {
//                    addEdge(x,b,distance_from[x][b]);
//                    for(int a{0}; a<n; ++a) {
//                        if(distance_from[a][x] + distance_from[y][b] < distance_from[a][b]) {
//                            addEdge(a,b,distance_from[a][y] + distance_from[y][b]);
//                        }
//                    }
//                }
//                
//                update(y,x,m_solver.core.backward(),distance_to[y]);
//                
//                
//                for(auto a : changed) {
//                    addEdge(a,y,distance_to[y][a]);
//                    for(int b{0}; b<n; ++b) {
//                        if(distance_from[a][x] + distance_from[y][b] < distance_from[a][b]) {
//                            addEdge(a,b,distance_from[a][y] + distance_from[y][b]);
//                        }
//                    }
//                }
                
                

                printMatrix();
            }
        }

//        printMatrix();
        exit(1);
}

template <typename T> FullTransitivity<T>::~FullTransitivity() {}


//template <typename T>
//void FullTransitivity<T>::addResource(const NoOverlapExpression<T>& res) {
//    for(auto disjunct{res.begDisjunct()}; disjunct != res.endDisjunct(); ++disjunct) {
//        auto l{m_solver.boolean.getLiteral(true, disjunct->id())};
//        auto prec_true{m_solver.boolean.getEdge(l)};
//        auto prec_false{m_solver.boolean.getEdge(~l)};
//        resource_graph[prec_true.from].add(prec_true.to);
//        resource_graph[prec_true.to].add(prec_true.from);
//        literal[prec_true.from][prec_true.to] = l;
//        literal[prec_false.from][prec_false.to] = ~l;
//        if(prec_true.from != prec_false.to or prec_true.to != prec_false.from) {
//            resource_graph[prec_false.from].add(prec_false.to);
//            resource_graph[prec_false.to].add(prec_false.from);
//        }
//    }
//}


template <typename T>
template <typename Iter>
void FullTransitivity<T>::addResource(const Iter beg_disjunct, const Iter end_disjunct) {
    for(auto disjunct{beg_disjunct}; disjunct != end_disjunct; ++disjunct) {
        auto l{m_solver.boolean.getLiteral(true, disjunct->id())};
        auto prec_true{m_solver.boolean.getEdge(l)};
        auto prec_false{m_solver.boolean.getEdge(~l)};
//        resource_graph[prec_true.from].add(prec_true.to);
//        resource_graph[prec_true.to].add(prec_true.from);
        literal[prec_true.from][prec_true.to] = l;
        literal[prec_false.from][prec_false.to] = ~l;
//        if(prec_true.from != prec_false.to or prec_true.to != prec_false.from) {
//            resource_graph[prec_false.from].add(prec_false.to);
//            resource_graph[prec_false.to].add(prec_false.from);
//        }
    }
}



template <typename T> void FullTransitivity<T>::post(const int idx) {

    Constraint<T>::cons_id = idx;
    Constraint<T>::idempotent = false;

#ifdef DEBUG_CONSTRAINT
  if (debug_flag > 0) {
    std::cout << "post " << *this << std::endl;
  }
#endif

//    for(auto x : resource_graph) {
//        for(auto y : resource_graph[x]) {
//            m_solver.wake_me_on(literal[x][y]);
//        }
//    }
    
    int n{static_cast<int>(_index_.size())};
    for(int x{0}; x<n; ++x) {
//        m_solver.wake_me_on(lb<T>(x), this->id());
//        m_solver.wake_me_on(ub<T>(x), this->id());
        for(int y{0}; y<n; ++y) {
            if(literal[x][y] != Solver<T>::Contradiction)
                m_solver.wake_me_on(literal[x][y], this->id());
        }
    }
}



template <typename T>
bool FullTransitivity<T>::notify(const Literal<T> l, const int) {

#ifdef DBG_FTRANS
  if (DBG_FTRANS) {
    std::cout << "\nnotify " << m_solver.pretty(l) << std::endl;
  }
#endif
    
//    if(l.isNumeric()) {
//        if(l.sign() == bound::lower) {
//            if(not changed_lb.has(l.variable()))
//                changed_lb.add(l.variable());
//        } else {
//            if(not changed_ub.has(l.variable()))
//                changed_ub.add(l.variable());
//        }
//    }
    
    
//    int n{static_cast<int>(_index_.size())};
    auto e{m_solver.boolean.getEdge(l)};
    
    if(addEdge(e.from, e.to, e.distance)) {
        propagateForward(edges.back());
        propagateBackward(edges.back());
    }
    
//    
//    update(e.from, e.to, m_solver.core, distance_from[e.from]);
//    
//    
//#ifdef DBG_FTRANS
//  if (DBG_FTRANS) {
//      std::cout << "new distances from x" << e.from << ":\n";
//  }
//#endif
//    
//    
//    for(auto y : changed) {
//        addEdge(e.from,y,distance_from[e.from][y],{l,Constant::NoIndex,_index_[e.to][y]});
//        for(int x{0}; x<n; ++x) {
//            if(distance_from[x][e.from] + distance_from[e.from][y] < distance_from[x][y]) {
//                addEdge(x,y,distance_from[x][e.from] + distance_from[e.from][y],{l,_index_[x][e.from],_index_[e.to][y]});
//            }
//        }
//    }
//    
//    
//    update(e.to, e.from, m_solver.core.backward(), distance_to[e.to]);
//    
//
//#ifdef DBG_FTRANS
//  if (DBG_FTRANS) {
//      std::cout << "new distances to x" << e.from << ":\n";
//  }
//#endif
//    
//    for(auto x : changed) {
//        addEdge(x,e.to,distance_to[e.to][x],{l,_index_[x][e.from],Constant::NoIndex});
//        for(int y{0}; y<n; ++y) {
//            if(distance_from[x][e.from] + distance_from[e.from][y] < distance_from[x][y]) {
//                addEdge(x,y,distance_from[x][e.from] + distance_from[e.from][y],{l,_index_[x][e.from],_index_[e.to][y]});
//            }
//        }
//    }
    
  return false;
}



template <typename T>
void FullTransitivity<T>::propagateForward(const DistanceConstraint<T> e, const Literal<T> l) {
    update(e.from, e.to, m_solver.core, distance_from[e.from]);
    
    
#ifdef DBG_FTRANS
  if (DBG_FTRANS) {
      std::cout << "new distances from x" << e.from << ":\n";
  }
#endif
    
    int n{static_cast<int>(_index_.size())};
    for(auto y : changed) {
        if(l == Solver<T>::Contradiction)
            addEdge(e.from,y,distance_from[e.from][y]);
        else
            addEdge(e.from,y,distance_from[e.from][y],{l,Constant::NoIndex,_index_[e.to][y]});
        for(int x{0}; x<n; ++x) {
            if(distance_from[x][e.from] + distance_from[e.from][y] < distance_from[x][y]) {
                if(l == Solver<T>::Contradiction)
                    addEdge(x,y,distance_from[x][e.from] + distance_from[e.from][y]);
                else
                    addEdge(x,y,distance_from[x][e.from] + distance_from[e.from][y],{l,_index_[x][e.from],_index_[e.to][y]});
            }
        }
    }
}


template <typename T>
void FullTransitivity<T>::propagateBackward(const DistanceConstraint<T> e, const Literal<T> l) {
    
    update(e.to, e.from, m_solver.core.backward(), distance_to[e.to]);
    

#ifdef DBG_FTRANS
  if (DBG_FTRANS) {
      std::cout << "new distances to x" << e.from << ":\n";
  }
#endif
    
    int n{static_cast<int>(_index_.size())};
    for(auto x : changed) {
        if(l == Solver<T>::Contradiction)
            addEdge(x,e.to,distance_to[e.to][x]);
        else
            addEdge(x,e.to,distance_to[e.to][x],{l,_index_[x][e.from],Constant::NoIndex});
        for(int y{0}; y<n; ++y) {
            if(distance_from[x][e.from] + distance_from[e.from][y] < distance_from[x][y]) {
                if(l == Solver<T>::Contradiction)
                    addEdge(x,y,distance_from[x][e.from] + distance_from[e.from][y]);
                else
                    addEdge(x,y,distance_from[x][e.from] + distance_from[e.from][y],{l,_index_[x][e.from],_index_[e.to][y]});
            }
        }
    }
}


// update the distance
template <typename T>
template <typename G>
void FullTransitivity<T>::update(const int x, const int y, const G &neighbors, std::vector<T>& shortest_path) {

  changed.clear();
  changed.add(y);

#ifdef DBG_BELLMAN_FT
  int max_iter = 1000;
  if (DBG_BELLMAN_FT) {
    std::cout << m_solver.core << "\nstart explore from " << y << std::endl;
  }
#endif

  while (not changed.empty()) {

    auto u{changed.front()};
    changed.pop_front();

#ifdef DBG_BELLMAN_FT
    if (max_iter-- < 0)
      exit(1);
    if (DBG_BELLMAN_FT) {
      std::cout << "pop " << u << " d=" << shortest_path[u] << ", q=(";
      for (auto evt : changed)
        std::cout << " " << evt;
      std::cout << " )" << std::endl;
    }
#endif

    for (auto edge : neighbors[u]) {
      int v{edge};
      auto w{edge.label()};

      if (shortest_path[u] + w < shortest_path[v]) {

#ifdef DBG_BELLMAN_FT
        if (DBG_BELLMAN_FT) {
          std::cout << " * shorter path -> "
                    << v << "("
                    << (shortest_path[u] + w) << "/" << shortest_path[v] << ")"
                    << std::endl;
        }
#endif

        if (v == x) {

#ifdef DBG_FAIL
          if (DBG_FAIL) {
            std::cout << " negative cyle\n";
          }
#endif

          throw Failure<T>(
              {this, static_cast<hint>(_index_[y][x])});
        }
        shortest_path[v] = shortest_path[u] + w;

        if (not changed.has(v))
          changed.add(v);
      }
#ifdef DBG_BELLMAN_FT
      else if (DBG_BELLMAN_FT) {
        std::cout << " ignore " << v << std::endl;
      }
#endif
    }
  }
}



template <typename T> void FullTransitivity<T>::propagate() {

    
#ifdef DBG_FTRANS
    if (DBG_FTRANS) {
        std::cout << "\npropagate " << *this << std::endl;
    }
#endif

}

//template <typename T> int FullTransitivity<T>::getType() const {
//  return EXPL;
//}

template <typename T>
void FullTransitivity<T>::xplain(const Literal<T>, const hint,
                             std::vector<Literal<T>> &) {

#ifdef DBG_EXPL_FTRANS
    std::cout << "explain " << l << " with FullTransitivity constraint (hint=" << h
              << ")\n";
#endif

}

template <typename T>
std::ostream &FullTransitivity<T>::display(std::ostream &os) const {
  os << "FullTransitivity";

#ifdef DEBUG_CONSTRAINT
  os << "[" << cons_id << "]";
#endif


  return os;
}

template <typename T>
void FullTransitivity<T>::printMatrix() const {
    var_t endv{static_cast<var_t>(_index_.size())};
    std::cout << "     ";
    for(var_t x{0}; x<endv; ++x) {
        std::cout << "x" << std::setw(4) << std::left << x;
    }
    std::cout << std::endl;
    for(var_t x{0}; x<endv; ++x) {
        std::cout << "x" << std::setw(4) << std::left << x;
        for(var_t y{0}; y<endv; ++y) {
            assert(distance_from[x][y] == distance_to[y][x]);
            assert(distance_from[x][y] == distance(x,y));
            if(distance_from[x][y] == Constant::Infinity<T>)
                std::cout << "    .";
            else if(distance_from[x][y] > Constant::Infinity<T>/2)
                std::cout << "    *";
                //std::cout << "i-" << std::setw(2) << (Constant::Infinity<T> - distance_from[x][y]);
            else
                std::cout << std::right << std::setw(5) << distance_from[x][y];
        }
        std::cout << std::endl;
    }
}

template <typename T>
std::ostream &FullTransitivity<T>::print_reason(std::ostream &os,
                                            const hint) const {
  os << "FullTransitivity";
  return os;
}


} // namespace tempo

#endif
