/************************************************
 * Tempo Solver.hpp
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

#ifndef _TEMPO_DIRECTEDGRAPH_HPP
#define _TEMPO_DIRECTEDGRAPH_HPP

#include "Global.hpp"
#include "Constant.hpp"
#include "ReversibleObject.hpp"
#include "util/SparseSet.hpp"
#include <boost/dynamic_bitset.hpp>
#include <vector>

namespace tempo {

//! Arc with a label info (of type T)
template<typename T>
class LabeledEdge {

public:
  LabeledEdge(const int v, const T l) : endpoint_(v), label_(l) {}
  virtual ~LabeledEdge() = default;
  LabeledEdge(const LabeledEdge<T> &e) = default;
  LabeledEdge(LabeledEdge<T> &&e) noexcept = default;

  LabeledEdge<T> &operator=(const LabeledEdge<T> &e) {
    endpoint_ = static_cast<int>(e);
    label_ = e.label();
    return *this;
  }

  inline operator int() const { return endpoint_; }
  inline T label() const { return label_; }
  inline void operator=(const int x) { endpoint_ = x; }

  template <class printer>
  std::ostream &display(std::ostream &os, printer p) const {
    os << p(*this);
    return os;
  }

protected:
  int endpoint_{-1};
  T label_{Constant::Infinity<T>};
};

//! Arc with a label info (of type T), and another field (of type S)
template<typename T, typename S>
class StampedLabeledEdge {

public:
  StampedLabeledEdge(const int v, const T l, const S s = 0)
      : endpoint_(v), label_(l), stamp_(s) {}
  virtual ~StampedLabeledEdge() = default;
  StampedLabeledEdge(const StampedLabeledEdge<T, S> &e) = default;
  StampedLabeledEdge(StampedLabeledEdge<T, S> &&e) noexcept = default;

  StampedLabeledEdge<T, S> &operator=(const StampedLabeledEdge<T, S> &e) {
    endpoint_ = static_cast<int>(e);
    label_ = e.label();
    return *this;
  }

  inline operator int() const { return endpoint_; }
  inline T label() const { return label_; }
  inline S stamp() const { return stamp_; }
  inline void operator=(const int x) { endpoint_ = x; }

  template <class printer>
  std::ostream &display(std::ostream &os, printer p) const {
    os << p(*this);
    return os;
  }

protected:
  int endpoint_{-1};
  T label_{Constant::Infinity<T>};
  S stamp_{0};
};

//! Directed Graph
/*
Implementation of a DirectedGraph
- It supports
* [adding new Arcs/not really];
* reducing Arcs' weights;
* and merging vertices
- It is Reversible
- Arc are parameter
An arc is in a neighborhood list of a vertex x. Arc (x,y) MUST implement
'operator int()' which should return 'y'. Therefore, in a DirectedGraph<int> arc
are simply encoded as the end-point
*/
template<class Arc>
class DirectedGraph : public ReversibleObject {

private:
  // currently active vertex
  SparseSet<int> vertices;

  // the adjency matrices (one for OUTs, and one for INs)
  std::vector<std::vector<Arc>> neighbor[2];
  // for each vertex and for IN/OUT : is it currently active
  std::vector<bool> active[2];

public:
  const std::vector<Arc> &operator[](const int u) const {
    return neighbor[OUT][u];
  }
  const std::vector<unsigned> &rank(const int u) const {
    return neighbor_rank[OUT][u];
  }
  std::vector<int>::const_iterator begin() const { return vertices.begin(); }
  std::vector<int>::const_iterator end() const { return vertices.end(); }

  const std::vector<std::vector<Arc>> &forward() const { return neighbor[OUT]; }
  const std::vector<std::vector<Arc>> &backward() const { return neighbor[IN]; }

  //    const std::vector<std::vector<Arc>> &backward() const { return
  //    neighbor[IN]; }
  bool has(const int x) const { return vertices.has(x); }
  bool is_active(const int x, const bool l) const { return active[l][x]; }

  size_t arcCount() const;
  size_t vertexCount() const;
  size_t size() const;
  size_t outdegree(const int u) const;
  size_t indegree(const int u) const;

  DirectedGraph(BacktrackEnvironment *e = ReversibleObject::env);
  DirectedGraph(const int n_vertices,
                BacktrackEnvironment *e = ReversibleObject::env);
  void resize(const int n);
  void newVertex(const int x);

  // add an Arc from origin (reversible)
  template <typename... T> void emplace_edge(const int, const T... args);
  void add(const int, const Arc &);

  // remove a vertex (reversible)
  void remove(const int);

  // remove only from IN/OUT
  void remove(const int, const bool);

  // merge the second vertex to the first (reversible)
  void merge(const int, const int);

  // restore to the last saved state
  void undo() override;

  Arc reverse(const Arc &, const int);

  template <class vprinter, class eprinter>
  std::ostream &display(std::ostream &os, vprinter vp, eprinter ep) const;

protected:
  // the rank of vertex v in its neighbor's adjency list
  std::vector<std::vector<unsigned>> neighbor_rank[2];

  // number of arcs
  size_t numArc{0};

  // The removed vertex is the first one out of vertices. The first int in
  // trail
  // gives the number of new forward Arcs for its representant x, the
  // second
  // the number of backward Arcs. They must be last in x's neighbor list
  std::vector<int> trail;

  // add and
  void add_(const int, const Arc &);

  // remove the arc between two vertices
  void undo_(const int, const int);

  // remove i-th neighbor (OUT or IN depending on l) from u's
  void removeIthNeighbor(const int, const int, const int);

  // re-add i-th neighbor (OUT or IN depending on l) from u's
  void addIthNeighbor(const int, const int, const int);

  // re-add vertex u previously removed
  void recall(const int);

  void recall(const int, const bool);

  void removeVertex(const int, const bool);

private:
  boost::dynamic_bitset<> buffer;

#ifdef DEBUG_SG
  void verify(const char *msg);
#endif
};


// declare a new time point in the DirectedGraph
template<class Arc>
DirectedGraph<Arc>::DirectedGraph(BacktrackEnvironment *e) : ReversibleObject(e) {}

template<class Arc>
DirectedGraph<Arc>::DirectedGraph(const int n_vertices, BacktrackEnvironment *e) : ReversibleObject(e) {
  resize(n_vertices);
}

template<class Arc>
Arc DirectedGraph<Arc>::reverse(const Arc& a, const int origin) {
  Arc r{a};
  r = origin;
  return r;
}

template<class Arc>
template <typename... T>
void DirectedGraph<Arc>::emplace_edge(const int origin, const T ...args) {
  Arc a(args...);
  add(origin, a);
}

template<class Arc>
void DirectedGraph<Arc>::add(const int origin, const Arc& a) {

  ReversibleObject::save();
  trail.push_back(origin);
  trail.push_back(a);

  add_(origin, a);

#ifdef DEBUG_SG
  verify("add arc");
#endif

}

template <class Arc>
void DirectedGraph<Arc>::add_(const int origin, const Arc &a) {

  int destination = a;

  neighbor_rank[OUT][origin].push_back(
      static_cast<int>(neighbor[IN][destination].size()));

  neighbor_rank[IN][a].push_back(
      static_cast<int>(neighbor[OUT][origin].size()));

  neighbor[OUT][origin].emplace_back(a);
  neighbor[IN][destination].emplace_back(reverse(a, origin));

  ++numArc;
}

template<class Arc>
void DirectedGraph<Arc>::remove(const int u) {

ReversibleObject::save();
trail.push_back(u);

for (auto l{OUT}; l <= IN; ++l) {
  removeVertex(u, l);
}

vertices.remove_front(u);

#ifdef DEBUG_SG
verify("remove vertex");
#endif
}

template<class Arc>
void DirectedGraph<Arc>::remove(const int u, const bool l) {

ReversibleObject::save();
trail.push_back(u);

removeVertex(u, l);

active[l][u] = false;
//  vertices.remove_front(u);

#ifdef DEBUG_SG
verify("remove vertex");
#endif
}

template<class Arc>
void DirectedGraph<Arc>::removeVertex(const int u, const bool l) {

//    std::cout << "remove " << u << " (" << l << ")\n";
//    for(auto v : neighbor[l][u]) {
//        std::cout << v << ":";
//        for(auto w : neighbor[1-l][v])
//            std::cout << " " << w;
//        std::cout << std::endl;
//    }

for (auto i{neighbor[l][u].size()}; i-- > 0;) {
  removeIthNeighbor(u, l, i);
  --numArc;
}

//    std::cout << "\n==>\n";
//    for(auto v : neighbor[l][u]) {
//        std::cout << v << ":";
//        for(auto w : neighbor[1-l][v])
//            std::cout << " " << w;
//        std::cout << std::endl;
//    }
}

template<class Arc>
void DirectedGraph<Arc>::removeIthNeighbor(const int u, const int l, const int ith) {
  // let l=OUT, we want to remove the IN arc corresponding to the ith OUT arc
  // from u

  auto a{neighbor[l][u][ith]};
  int v{a};
  auto kth{neighbor_rank[l][u][ith]};
  // the arc is (u -> v) and is kth in v's IN neighbors

  assert(static_cast<unsigned>(ith) == neighbor_rank[1 - l][v][kth]);
  assert(u == static_cast<int>(neighbor[1 - l][v][kth]));

  if (kth < (neighbor[1 - l][v].size() - 1)) {
    // if not already at the back, we need to move the back arc into (u -> v)'s
    // place

    auto back_arc{neighbor[1 - l][v].back()};
    int w{back_arc};
    auto jth{neighbor_rank[1 - l][v].back()};
    // back arc is (w -> v) and is jth in w's neighbors

    assert(neighbor_rank[l][w][jth] == neighbor[1 - l][v].size() - 1);

    // copy the back arc at rank k
    neighbor[1 - l][v][kth] = back_arc;
    neighbor_rank[1 - l][v][kth] = jth;
    neighbor_rank[l][w][jth] = kth;
  }

  // pop the back arc
  neighbor[1 - l][v].pop_back();
  neighbor_rank[1 - l][v].pop_back();
}

template<class Arc>
void DirectedGraph<Arc>::addIthNeighbor(const int u, const int l, const int ith) {
  // let l=OUT, we want to add the IN arc corresponding to the ith OUT arc from
  // u, at the rank it was before
  auto a{neighbor[l][u][ith]};
  int v{a};
  auto kth{neighbor_rank[l][u][ith]};
  // the arc is (u -> v) and it was kth in v's IN  neighbors

  if (kth < neighbor[1 - l][v].size()) {
    // we need to replace it at rank kth (via a swap)

    auto swap_arc{neighbor[1 - l][v][kth]};
    int w{swap_arc};
    auto jth{neighbor_rank[1 - l][v][kth]};
    // back arc is (w -> v) and is jth in w's neighbors

    assert(neighbor_rank[l][w][jth] == kth);
    assert(static_cast<int>(neighbor[l][w][jth]) == v);

    neighbor_rank[l][w][jth] = neighbor[1 - l][v].size();

    neighbor[1 - l][v].push_back(swap_arc);
    neighbor_rank[1 - l][v].push_back(jth);

    neighbor[1 - l][v][kth] = reverse(a, u);
    neighbor_rank[1 - l][v][kth] = ith;
  } else {

    assert(kth == neighbor[1 - l][v].size());

    neighbor[1 - l][v].push_back(reverse(a, u));
    neighbor_rank[1 - l][v].push_back(ith);
  }

//
//    // it's going to end up at rank neighbor_rank[1-l][v].size()
//    neighbor_rank[l][u][ith] = neighbor_rank[1-l][v].size();
//    
//    // push it back there
//    neighbor_rank[1-l][v].emplace_back(ith);
//    neighbor[1-l][v].emplace_back(reverse(a,u));
}

template<class Arc>
void DirectedGraph<Arc>::recall(const int u) {
  for (auto l{IN}; l >= OUT; --l) {
    int n{static_cast<int>(neighbor[l][u].size())};
    for (int ith{0}; ith < n; ++ith) {
      addIthNeighbor(u, l, ith);
      ++numArc;
    }
  }

  vertices.add(u);

#ifdef DEBUG_SG
  verify("recall");
#endif
}

template<class Arc>
void DirectedGraph<Arc>::recall(const int u, const bool l) {

  int n{static_cast<int>(neighbor[l][u].size())};
  for (int ith{0}; ith < n; ++ith) {
    addIthNeighbor(u, l, ith);
    ++numArc;
  }

  active[l][u] = true;

#ifdef DEBUG_SG
  verify("recall");
#endif
}

template <class Arc> void DirectedGraph<Arc>::undo_(const int u, const int v) {

  assert(static_cast<int>(neighbor[OUT][u].back()) == v);
  neighbor[OUT][u].pop_back();
  neighbor_rank[OUT][u].pop_back();

  assert(static_cast<int>(neighbor[IN][v].back()) == u);
  neighbor[IN][v].pop_back();
  neighbor_rank[IN][v].pop_back();

  --numArc;
}

template<class Arc>
void DirectedGraph<Arc>::resize(const int n) {

  if (n > 0) {

    buffer.resize(n, 0);

    neighbor[IN].resize(n);
    neighbor[OUT].resize(n);
    neighbor_rank[IN].resize(n);
    neighbor_rank[OUT].resize(n);

    vertices.reserve(n);
    vertices.fill_back();

    active[IN].resize(n, true);
    active[OUT].resize(n, true);
  }

#ifdef DEBUG_SG
  // cout << *this << std::endl;
  verify("init");
#endif
}


template<class Arc>
void DirectedGraph<Arc>::newVertex(const int x) {

//    std::cout << "new vertex " << x << std::endl;

assert(static_cast<size_t>(x) >= size());
auto n{size()};

resize(static_cast<size_t>(x + 1));

for (auto i{n}; i < size() - 1; ++i) {
  vertices.remove_front(static_cast<int>(i));
}

#ifdef DEBUG_SG
// cout << *this << std::endl;
verify("new vertex");
#endif
}

template<class Arc>
size_t DirectedGraph<Arc>::arcCount() const { return numArc; }

template<class Arc>
size_t DirectedGraph<Arc>::vertexCount() const { return vertices.size(); }

template<class Arc>
size_t DirectedGraph<Arc>::size() const { return vertices.capacity(); }

template<class Arc>
size_t DirectedGraph<Arc>::outdegree(const int u) const { return neighbor[OUT][u].size(); }

template<class Arc>
size_t DirectedGraph<Arc>::indegree(const int u) const { return neighbor[IN][u].size(); }

//
template<class Arc>
void DirectedGraph<Arc>::merge(const int u, const int v) {

#ifdef DEBUG_SG
  verify("before merge");
#endif

  assert(vertices.has(v));
  assert(vertices.has(u));

  remove(v);

  buffer.reset();
  buffer.set(u);
  for (auto a : neighbor[OUT][u]) {
    assert(static_cast<int>(a) >= 0);
    assert(static_cast<size_t>(static_cast<int>(a)) < buffer.size());
    buffer.set(a);
  }

  for (auto a : neighbor[OUT][v]) {
    if (not buffer[a]) {
      add(u, a);
    }
  }

  buffer.reset();
  buffer.set(u);
  for (auto a : neighbor[IN][u]) {
    assert(static_cast<int>(a) >= 0);
    assert(static_cast<size_t>(static_cast<int>(a)) < buffer.size());
    buffer.set(a);
  }

  for (auto a : neighbor[IN][v]) {
    if (not buffer[a]) {
      add(a, reverse(a, u));
    }
  }

#ifdef DEBUG_SG
  verify("after merge");
#endif
}

// restore to the last saved state
template<class Arc>
void DirectedGraph<Arc>::undo() {

  auto v{trail.back()};
  trail.pop_back();

  // cout << "trail.back(): " << v << " " << vertices.has(v) << " " << vertices
  // << std::endl;

  if (vertices.has(v)) {

    if (not active[IN][v])
      recall(v, IN);
    else if (not active[OUT][v])
      recall(v, OUT);
    else {
      // undo an addArc(u, v)
      auto u{trail.back()};
      trail.pop_back();
      undo_(u, v);
    }

  } else {

    // cout << "undo rm " << v << std::endl ; //<< *this << std::endl;

    // undo a remove
    recall(v);
  }

  // cout << Arcs << std::endl;

#ifdef DEBUG_SG
  // cout << *this << std::endl;
  verify("undo");
#endif
}

template<class Arc>
template<class vprinter, class eprinter>
std::ostream &DirectedGraph<Arc>::display(std::ostream &os, vprinter p_vertex, eprinter p_edge) const {

  for (unsigned i{0}; i < size(); ++i) {
    if (not has(i))
      continue;
    os << std::setw(3) << p_vertex(i) << ":";

#ifdef DEBUG_SG
    int k = 0;
#endif
    for (auto a : neighbor[OUT][i]) {
      os << " " << p_edge(a);
//            a.display(os,p);
#ifdef DEBUG_SG
      << "@" << neighbor_rank[OUT][i][k++]
#endif
          ;
    }

#ifdef DEBUG_SG
    os << " |";
    k = 0;
    for (auto a : neighbor[IN][i]) {
      os << " " << p_edge(a);
      //            a.display(os,p);
      << "@" << neighbor_rank[IN][i][k++];
    }
#endif
    os << std::endl;
  }

  return os;
}

#ifdef DEBUG_SG

template<class Arc>
void DirectedGraph<Arc>::verify(const char *msg) {

  // cout << "beg verif " << msg << std::endl ;
  for (auto l{OUT}; l <= IN; ++l) {
    size_t count = 0;
    for (auto x : vertices) {

      // if(x==188)
      //     cout << "verif " << x << std::endl;
      //
      // int rank = 0;
      for (size_t ith{0}; ith < neighbor[l][x].size(); ++ith) {
        int y{neighbor[l][x][ith]};
        if (y == x) {
          std::cout << msg << ": error, there is a loop on " << x << std::endl;
          exit(1);
        }

        auto jth{neighbor_rank[l][x][ith]};

        // if(x==188)
        //     cout << ith << ": " << y << " (" << jth << " ==> " <<  neighbor[1
        //     - l][y][jth] << "/" << neighbor_rank[1 - l][y][jth] << ")" <<
        //     std::endl;
        //
        if (neighbor_rank[1 - l][y][jth] != ith or
            static_cast<int>(neighbor[1 - l][y][jth]) != x) {
          std::cout << msg << ": error on " << x << " " << y << "\n"
                    << *this << std::endl;
          exit(1);
        }

        // ++rank;
      }
      count += neighbor[l][x].size();
    }

    assert(count == arcCount());
  }

  // cout << "end verif " << msg << std::endl ;
}

#endif



template<class Arc>
std::ostream &operator<<(std::ostream &os, const DirectedGraph<Arc> &x) {
  return x.display(
      os, [](const int i) { return static_cast<int>(i); },
      [](const Arc &e) {
        std::stringstream ss;
        ss << e;
        return ss.str();
      });
  //        return static_cast<int>(e);});
}


template<typename T>
std::ostream &operator<<(std::ostream &os, const LabeledEdge<T> &x) {
  return x.display(os,
                   [](const LabeledEdge<T> &e) { return static_cast<int>(e); });
}

template<typename T, typename S>
std::ostream &operator<<(std::ostream &os, const StampedLabeledEdge<T,S> &x) {
  return x.display(os, [](const StampedLabeledEdge<T, S> &e) {
    return std::to_string(static_cast<int>(e)) + " (" +
           std::to_string(e.label()) + ")";
  });
}


//int concatenate(int, int b, DirectedGraph<int>&) {
//    return b;
//}
//
//template<typename T>
//LabeledEdge<T> concatenate(LabeledEdge<T>& a, LabeledEdge<T>& b, DirectedGraph<LabeledEdge<T>>&) {
//    return LabeledEdge<T>(static_cast<int>(b), a.label()+b.label());
//}
//
//template<typename T>
//StampedLabeledEdge<T> concatenate(StampedLabeledEdge<T>& a, StampedLabeledEdge<T>& b, DirectedGraph<StampedLabeledEdge<T>>& g) {
//    return StampedLabeledEdge<T>(static_cast<int>(b), a.label()+b.label(), g.arcCount());
//}


} // namespace tempo

#endif
