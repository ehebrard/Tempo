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

template <typename T>
class FullTransitivity : public Constraint<T>, public ReversibleObject {
private:
  Solver<T> &m_solver;

  // all the edges, in chronological order
  std::vector<DistanceConstraint<T>> edges;
  // the reason for each edge[: a literal, and the leading and following
  // sub-paths]
  std::vector<PathExplanation<T>> reason;
  // for each pair x,y, the index of the last literal
  std::vector<std::vector<index_t>> _index_;
  // for each edge, the index of the previous edge
  std::vector<index_t> previous;

  // the distance from
  std::vector<std::vector<T>> distance_from;
  // the distance to (transpose of 'distance_from')
  std::vector<std::vector<T>> distance_to;

  // for each pair x,y, the corresponding disjunct literal, or
  // Solver::Contradiction<T> if there is none
  std::vector<std::vector<Literal<T>>> literal;

  // to store bounds events (so that they are all propagated together)
  SparseSet<> changed_lbs;
  SparseSet<> changed_ubs;

  // util for Bellman-Ford
  SparseSet<> changed;
  SparseSet<> enqueued;

  bool direction_flag{false};
  int ref{-1};

public:
  FullTransitivity(Solver<T> &solver);
  virtual ~FullTransitivity();

  void clear();

  // add a set of disjuncts
  template <typename Iter>
  void addResource(const Iter beg_disjunct, const Iter end_disjunct);

  bool notify(const Literal<T>, const int rank) override;
  void post(const int idx) override;
  void propagate() override;

  T distance(const int x, const int y) const;
  bool addEdge(const int x, const int y, const T d,
               const PathExplanation<T> r = NoReason);
  void undo() override;

  void propagateForward(const DistanceConstraint<T> e,
                        const Literal<T> l = Solver<T>::Contradiction);
  void propagateBackward(const DistanceConstraint<T> e,
                         const Literal<T> l = Solver<T>::Contradiction);

  template <typename G>
  void update(const int x, const int y, const G &neighbors,
              std::vector<T> &distance,
              const std::vector<std::vector<T>> &current_distance);

  void xplain(const Literal<T> l, const hint h,
              std::vector<Literal<T>> &Cl) override;

  // helper for xplain
  void getPath(const index_t h, std::vector<Literal<T>> &Cl);

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

  void printMatrix() const;

  static PathExplanation<T> NoReason;

#ifdef DBG_FTRANS
  void verify(const char *msg) const;
#endif
};

template <typename T>
PathExplanation<T> FullTransitivity<T>::NoReason = PathExplanation<T>();


template <typename T>
T FullTransitivity<T>::distance(const int x, const int y) const {
  return edges[_index_[x][y]].distance;
}

template <typename T>
bool FullTransitivity<T>::addEdge(const int x, const int y, const T d, const PathExplanation<T> r) {

  if (distance(x, y) <= d)
    return false;

#ifdef DBG_FTRANS
  if (DBG_FTRANS) {
    std::cout << "  - x" << x << " -> x" << y << " (" << distance(x, y) << "/"
              << d << ")\n";
  }
#endif

  ReversibleObject::save();

  auto i{static_cast<index_t>(edges.size())};

  distance_from[x][y] = distance_to[y][x] = d;

  previous.push_back(_index_[x][y]);

  _index_[x][y] = i;
  edges.emplace_back(x, y, d);
  reason.push_back(r);

  if (literal[y][x] != Solver<T>::Contradiction and
      m_solver.boolean.getEdge(literal[y][x]).distance + distance_from[x][y] <
          0) {

#ifdef DBG_FTRANS
    if (DBG_FTRANS) {
      std::cout << "  --> infer literal " << m_solver.pretty(literal[x][y])
                << "\n";
    }
#endif

    try {
      m_solver.set(literal[x][y], {this, static_cast<hint>(i)});
    } catch (Failure<T> &f) {
      clear();
      throw f;
    }
  }

  return true;
}

template <typename T> void FullTransitivity<T>::undo() {

  auto e{edges.back()};
  auto i{previous.back()};

  assert(i == 0 or edges[i].from == e.from and edges[i].to == e.to);

  _index_[e.from][e.to] = i;
  distance_from[e.from][e.to] = distance_to[e.to][e.from] = edges[i].distance;

  reason.pop_back();
  edges.pop_back();
  previous.pop_back();

  changed_lbs.clear();
  changed_ubs.clear();
}

template <typename T>
FullTransitivity<T>::FullTransitivity(Solver<T> &solver)
    : ReversibleObject(&solver.getEnv()), m_solver(solver) {

  Constraint<T>::priority = Priority::Low;

  auto n{m_solver.numeric.size()};

  _index_.resize(n);
  distance_from.resize(n);
  distance_to.resize(n);
  literal.resize(n);

  changed.reserve(n);
  enqueued.reserve(n);
  changed_lbs.reserve(n);
  changed_ubs.reserve(n);

  edges.push_back({Constant::NoVar, Constant::NoVar, Constant::Infinity<T>});
  reason.push_back(NoReason);
  for (auto &row : _index_)
    while (row.size() < n)
      row.emplace_back(0);

  for (auto &row : distance_from)
    while (row.size() < n)
      row.emplace_back(Constant::Infinity<T>);
  for (auto &row : distance_to)
    while (row.size() < n)
      row.emplace_back(Constant::Infinity<T>);

  for (auto &row : literal)
    row.resize(n, Solver<T>::Contradiction);

  for (size_t x{0}; x < _index_.size(); ++x) {
    _index_[x][x] = Constant::NoIndex;
    distance_from[x][x] = 0;
    distance_to[x][x] = 0;
  }

  // HACK!!: add an edge from origin to end to allow propagation of
  // shirtest path
  m_solver.set({0, 1, m_solver.numeric.upper(1)});

#ifdef DBG_FTRANS
  if (DBG_FTRANS) {
    printMatrix();
  }
#endif

  for (auto x : m_solver.core) {

    if (x != 0) {
      if (m_solver.numeric.upper(x) != Constant::Infinity<T>) {

#ifdef DBG_FTRANS
        if (DBG_FTRANS) {
          std::cout << "\n**upper bound of " << x << std::endl;
        }
#endif

        if (addEdge(0, x, m_solver.numeric.upper(x))) {
          auto f{edges.back()};
          propagateForward(f);
          propagateBackward(f);
        }

#ifdef DBG_FTRANS
        if (DBG_FTRANS) {
          printMatrix();
        }
#endif
      }

      if (m_solver.numeric.lower(x) != -Constant::Infinity<T>) {

#ifdef DBG_FTRANS
        if (DBG_FTRANS) {
          std::cout << "\n**lower bound of " << x << std::endl;
        }
#endif

        if (addEdge(x, 0, -m_solver.numeric.lower(x))) {
          auto f{edges.back()};
          propagateForward(f);
          propagateBackward(f);
        }

#ifdef DBG_FTRANS
        if (DBG_FTRANS) {
          printMatrix();
        }
#endif
      }
    }

    for (auto e : m_solver.core[x]) {

      int y{e};

#ifdef DBG_FTRANS
      if (DBG_FTRANS) {
        std::cout << "\n**edge " << x << " -> " << y << std::endl;
      }
#endif

      if (addEdge(x, y, e.label())) {
        auto f{edges.back()};
        propagateForward(f);
        propagateBackward(f);
      }

#ifdef DBG_FTRANS
      if (DBG_FTRANS) {
        printMatrix();
      }
#endif
    }
  }

#ifdef DBG_FTRANS
  if (DBG_FTRANS) {
    printMatrix();
  }
#endif

#ifdef DBG_FTRANS
  verify("after constructor");
#endif
}

template <typename T> FullTransitivity<T>::~FullTransitivity() {}

template <typename T> void FullTransitivity<T>::clear() {

  if (direction_flag) {
    for (auto yp{changed.fbegin()}; yp != changed.fend(); ++yp) {
      auto y{*yp};
      assert(distance_from[ref][y] <= distance_to[y][ref]);
      assert(distance(ref, y) == distance_to[y][ref]);
      distance_from[ref][y] = distance_to[y][ref];
    }
  } else {
    for (auto xp{changed.fbegin()}; xp != changed.fend(); ++xp) {
      auto x{*xp};
      assert(distance_from[x][ref] >= distance_to[ref][x]);
      assert(distance(x, ref) == distance_from[x][ref]);
      distance_to[ref][x] = distance_from[x][ref];
    }
  }

  changed_lbs.clear();
  changed_ubs.clear();
}

template <typename T>
template <typename Iter>
void FullTransitivity<T>::addResource(const Iter beg_disjunct, const Iter end_disjunct) {
  for (auto disjunct{beg_disjunct}; disjunct != end_disjunct; ++disjunct) {
    auto l{m_solver.boolean.getLiteral(true, disjunct->id())};
    auto prec_true{m_solver.boolean.getEdge(l)};
    auto prec_false{m_solver.boolean.getEdge(~l)};
    literal[prec_true.from][prec_true.to] = l;
    literal[prec_false.from][prec_false.to] = ~l;
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

  int n{static_cast<int>(_index_.size())};
  for (int x{0}; x < n; ++x) {
    m_solver.wake_me_on(lb<T>(x), this->id());
    m_solver.wake_me_on(ub<T>(x), this->id());
    for (int y{0}; y < n; ++y) {
      if (literal[x][y] != Solver<T>::Contradiction)
        m_solver.wake_me_on(literal[x][y], this->id());
    }
  }
}



template <typename T>
bool FullTransitivity<T>::notify(const Literal<T> l, const int) {

#ifdef DBG_FTRANS
  if (DBG_FTRANS) {
//    if(m_solver.num_fails > 1000) {
        std::cout << "\nnotify " << m_solver.pretty(l) << std::endl;
//        exit(1);
//    }
  }
#endif

#ifdef DBG_FTRANS
  verify("before notify");
#endif

  if (l.isNumeric()) {
    if (l.sign() == bound::lower) {
      if (not changed_lbs.has(l.variable()))
        changed_lbs.add(l.variable());
    } else {
      if (not changed_ubs.has(l.variable()))
        changed_ubs.add(l.variable());
    }
    return true;
  }

  auto e{m_solver.boolean.getEdge(l)};

  if (addEdge(e.from, e.to, e.distance)) {
    auto f{edges.back()};
    propagateForward(f, l);
    propagateBackward(f, l);
  }

#ifdef DBG_FTRANS
  if (DBG_FTRANS) {
    printMatrix();
  }
#endif

#ifdef DBG_FTRANS
  verify("after notify ");
#endif

  return false;
}

template <typename T>
void FullTransitivity<T>::propagateForward(const DistanceConstraint<T> e,
                                           const Literal<T> l) {

#ifdef DBG_FTRANS
  verify("before prop forward");
#endif

  direction_flag = true;
  ref = e.from;

  update(e.from, e.to, m_solver.core, distance_from[e.from], distance_to);

#ifdef DBG_FTRANS
  if (DBG_FTRANS) {
    std::cout << "new distances (f) because of x" << e.from << " -> x" << e.to
              << ":\n";
  }
#endif

  int n{static_cast<int>(_index_.size())};

  assert(changed.empty());

  for (auto yp{changed.fbegin()}; yp != changed.fend(); ++yp) {

    auto y{*yp};
    for (int x{0}; x < n; ++x)
      if (x != y) {
        if ((distance_from[x][e.from] < Constant::Infinity<T>)and(
                distance_from[e.from][y] < Constant::Infinity<T>) and
            (distance_from[x][e.from] + distance_from[e.from][y] <
             distance_to[y][x])) {
          if (l == Solver<T>::Contradiction)
            addEdge(x, y, distance_from[x][e.from] + distance_from[e.from][y]);
          else
            addEdge(x, y, distance_from[x][e.from] + distance_from[e.from][y],
                    {l, _index_[x][e.from], _index_[e.to][y]});
        }
      }
  }

#ifdef DBG_FTRANS
  verify("after prop forward");
#endif
}

template <typename T>
void FullTransitivity<T>::propagateBackward(const DistanceConstraint<T> e,
                                            const Literal<T> l) {

#ifdef DBG_FTRANS
  verify("before prop backward");
#endif

  direction_flag = false;
  ref = e.to;

  update(e.to, e.from, m_solver.core.backward(), distance_to[e.to],
         distance_from);

#ifdef DBG_FTRANS
  if (DBG_FTRANS) {
    std::cout << "new distances (b) because of x" << e.from << " -> x" << e.to
              << ":\n";
  }
#endif

  int n{static_cast<int>(_index_.size())};
  for (auto xp{changed.fbegin()}; xp != changed.fend(); ++xp) {
    auto x{*xp};
    for (int y{0}; y < n; ++y)
      if (x != y) {
        if ((distance_to[e.to][x] < Constant::Infinity<T>)and(
                distance_to[y][e.to] < Constant::Infinity<T>) and
            (distance_to[e.to][x] + distance_to[y][e.to] <
             distance_from[x][y])) {
          if (l == Solver<T>::Contradiction)
            addEdge(x, y, distance_to[e.to][x] + distance_to[y][e.to]);
          else
            addEdge(x, y, distance_to[e.to][x] + distance_to[y][e.to],
                    {l, _index_[x][e.from], _index_[e.to][y]});
        }
      }
  }

#ifdef DBG_FTRANS
  if (DBG_FTRANS) {
    std::cout << "end (b)\n";
  }
#endif

#ifdef DBG_FTRANS
  verify("after prop backward");
#endif
}

// update the distance
template <typename T>
template <typename G>
void FullTransitivity<T>::update(
    const int x, const int y, const G &neighbors, std::vector<T> &shortest_path,
    const std::vector<std::vector<T>> &current_distance) {

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

      if (current_distance[v][u] != Constant::Infinity<T> and
          shortest_path[u] + w < shortest_path[v]) {

#ifdef DBG_BELLMAN_FT
        if (DBG_BELLMAN_FT) {
          std::cout << " * shorter path -> " << v << "("
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

          clear();
          throw Failure<T>({this, static_cast<hint>(_index_[y][x])});
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
//    if(m_solver.num_fails > 1000) {
        std::cout << "\npropagate " << *this << std::endl;
//        exit(1);
//    }
  }
#endif

#ifdef DBG_FTRANS
  verify("before propagate");
#endif

  while (not changed_lbs.empty()) {
    auto x{changed_lbs.back()};

    if (addEdge(x, 0, -m_solver.numeric.lower(x))) {
      auto f{edges.back()};
      auto lbl{m_solver.numeric.getLiteral(bound::lower, x)};
      propagateForward(f, lbl);
      propagateBackward(f, lbl);
    }

    changed_lbs.pop_back();
  }

  while (not changed_ubs.empty()) {
    auto x{changed_ubs.back()};

    if (addEdge(0, x, m_solver.numeric.upper(x))) {
      auto f{edges.back()};
      auto ubl{m_solver.numeric.getLiteral(bound::upper, x)};
      propagateForward(f, ubl);
      propagateBackward(f, ubl);
    }

    changed_ubs.pop_back();
  }

#ifdef DBG_FTRANS
  if (DBG_FTRANS) {
    printMatrix();
  }
#endif

#ifdef DBG_FTRANS
  verify("after propagate");
#endif
}

template <typename T>
void FullTransitivity<T>::getPath(const index_t h,
                                  std::vector<Literal<T>> &Cl) {
  if (h == Constant::NoIndex)
    return;

  auto pr{reason[h]};

  if (pr.literal == Solver<T>::Contradiction) {

#ifdef DBG_EXPL_FTRANS
    std::cout << edges[h];
#endif

    auto e{edges[h]};
    auto l{literal[e.from][e.to]};
    if (l != Solver<T>::Contradiction) {
      Cl.push_back(l);

#ifdef DBG_EXPL_FTRANS
      std::cout << " (external)\n";
#endif

    }

#ifdef DBG_EXPL_FTRANS
    else
      std::cout << " (ground)\n";
#endif

    //        return;
  } else {

    Cl.push_back(pr.literal);

#ifdef DBG_EXPL_FTRANS
    std::cout << m_solver.pretty(pr.literal) << " (external*)\n";
#endif
  }

  getPath(pr.prefix, Cl);
  getPath(pr.suffix, Cl);
}

template <typename T>
void FullTransitivity<T>::xplain(const Literal<T> l, const hint h,
                                 std::vector<Literal<T>> &Cl) {

  if (l == Solver<T>::Contradiction) {

#ifdef DBG_EXPL_FTRANS
    std::cout << "explain failure (negative cycle) with FullTransitivity "
                 "constraint (hint="
              << edges[h] << ")\n";
    exit(1);
#endif

  } else {

#ifdef DBG_EXPL_FTRANS
    std::cout << "explain " << m_solver.pretty(l)
              << " with FullTransitivity constraint (hint=" << edges[h]
              << ")\n";
#endif

    getPath(h, Cl);

    //#ifdef DBG_EXPL_FTRANS
    //        exit(1);
    //#endif
  }
}

template <typename T>
std::ostream &FullTransitivity<T>::display(std::ostream &os) const {
  os << "FullTransitivity";

#ifdef DEBUG_CONSTRAINT
  os << "[" << cons_id << "]";
#endif

  return os;
}

#ifdef DBG_FTRANS
template <typename T> void FullTransitivity<T>::verify(const char *msg) const {
  var_t endv{static_cast<var_t>(_index_.size())};
  for (var_t x{0}; x < endv; ++x) {
    for (var_t y{0}; y < endv; ++y) {
      if (distance_from[x][y] != distance_to[y][x]) {
        printMatrix();
        std::cout << " BUG " << msg << " (#cp=" << m_solver.num_choicepoints
                  << "): discrepancy to/from on " << x << "," << y << " ("
                  << distance_from[x][y] << " != " << distance_to[y][x]
                  << ") \n";
        exit(1);
      }
      if (x != y and (distance_from[x][y] != distance(x, y))) {
        printMatrix();
        std::cout << " BUG " << msg << " (#cp=" << m_solver.num_choicepoints
                  << "): discrepancy edges/matrix on " << x << "," << y << " ("
                  << distance_from[x][y] << " != " << distance(x, y) << ") \n";
        exit(1);
      }
    }
  }
}
#endif

template <typename T>
void FullTransitivity<T>::printMatrix() const {
  var_t endv{static_cast<var_t>(_index_.size())};
  std::cout << "     ";
  for (var_t x{0}; x < endv; ++x) {
    std::cout << "x" << std::setw(4) << std::left << x;
  }
  std::cout << std::endl;
  for (var_t x{0}; x < endv; ++x) {
    std::cout << "x" << std::setw(4) << std::left << x;
    for (var_t y{0}; y < endv; ++y) {
      if (distance_from[x][y] == Constant::Infinity<T>)
        std::cout << "    .";
      else if (distance_from[x][y] > Constant::Infinity<T> / 2)
        std::cout << "    *";
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
