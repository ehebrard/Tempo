/************************************************
 * Tempo Transitivity.hpp
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

#ifndef TEMPO_TRANSITIVITY_HPP
#define TEMPO_TRANSITIVITY_HPP

#include <cassert>
#include <map>
#include <vector>

#include "Explanation.hpp"
#include "Global.hpp"
#include "ReversibleObject.hpp"
#include "Solver.hpp"
#include "constraints/Constraint.hpp"
#include "util/DisjointSet.hpp"
#include "util/SparseSet.hpp"

//#define DBG_LTRANS

namespace tempo {

template <typename T> class Transitivity : public Constraint<T> {
private:
  Solver<T> &m_solver;
  Interval<T> schedule;
  std::vector<Interval<T>> the_tasks;

  // encoding: (y \in front(DAG[x]) && x \in back(y)) <=> edge (x,y)
  std::vector<SparseSet<int, Reversible<size_t>>> DAG;
  // encoding: y \in transitive_reduction[x] <=> edge (x,y)
  SparseSet<int, Reversible<size_t>> transitive_reduction;
  DisjointSet<> forest;
  std::vector<std::vector<T>> distance_matrix;

  std::vector<std::vector<Literal<T>>> disjunct;
  std::vector<int> scopex;
  std::vector<int> scopey;

  std::vector<int> sorted_tasks;
  std::vector<T> offset;

  std::vector<int> new_succ_of_x;
  std::vector<int> new_pred_of_y;

    // @TODO remove that
  std::vector<int> task_map;

  bool change_flag{false};

  bool transition_flag{false};

  T length(const size_t e) const {
      auto i{first(e)};
      auto j{second(e)};
      return transition_time(i,j);
  }
    size_t edge(const size_t i, const size_t j) const {
      return the_tasks.size() * i + j;
    }
  size_t first(const size_t e) const { return e / the_tasks.size(); }
  size_t second(const size_t e) const { return e % the_tasks.size(); }
  T transition_time(const int i, const int j) const;

public:
  template <typename ItTask, typename ItVar>
  Transitivity(Solver<T> &solver, Interval<T> &sched, const ItTask beg_task,
               const ItTask end_task, const ItVar beg_var);
  virtual ~Transitivity();

  void add_edge(const int x, const int y, const int r);

  bool notify(const Literal<T>, const int rank) override;
  void post(const int idx) override;
  void propagate() override;
  void min_spanning_tree();

  void xplain(const Literal<T> l, const hint h, std::vector<Literal<T>> &Cl) override;
  int getType() const override;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

  std::string prettyTask(const int i) const;

#ifdef DEBUG_CONSTRAINT
  int debug_flag{2};
#endif
};

template <typename T>
template <typename ItTask, typename ItVar>
Transitivity<T>::Transitivity(Solver<T> &solver, Interval<T> &sched,
                              const ItTask beg_task, const ItTask end_task,
                              const ItVar beg_var)
    : m_solver(solver), schedule(sched),
      transitive_reduction(std::distance(beg_task, end_task) *
                               std::distance(beg_task, end_task),
                           &m_solver.getEnv()) {

  Constraint<T>::priority = Priority::Low;

          task_map.resize(m_solver.numeric.size(), -1);

          for (auto jp{beg_task}; jp != end_task; ++jp) {
            int t{static_cast<int>(the_tasks.size())};
            task_map[jp->start.id()] = t;
            task_map[jp->end.id()] = t;
            the_tasks.push_back(*jp);
          }

  disjunct.resize(the_tasks.size());
  sorted_tasks.resize(the_tasks.size());
  offset.resize(the_tasks.size());

  transitive_reduction.fill();
  forest.resize(the_tasks.size());

  for (size_t i{0}; i < the_tasks.size(); ++i) {

    disjunct[i].resize(the_tasks.size());

    DAG.emplace_back(the_tasks.size(), &m_solver.getEnv());
    DAG.back().fill();

    transitive_reduction.remove_back(edge(i, i));
  }

  auto ep{beg_var};
  for (auto ip{beg_task}; ip != end_task; ++ip) {
    for (auto jp{ip + 1}; jp != end_task; ++jp) {
      auto x{*ep};

      auto i{std::distance(beg_task, ip)};
      auto j{std::distance(beg_task, jp)};
      disjunct[i][j] = m_solver.boolean.getLiteral(false, x);
      disjunct[j][i] = m_solver.boolean.getLiteral(true, x);

      //      if (not transition_flag) {
      //        transition_flag =
      //            ((transition_time(i, j) > 0) or (transition_time(j, i) >
      //            0));
      //      }
      ++ep;
    }
  }
}

template <typename T> Transitivity<T>::~Transitivity() {}

template <typename T>
T Transitivity<T>::transition_time(const int i, const int j) const {
  return -m_solver.boolean.getEdge(disjunct[i][j]).distance;
}

template <typename T> void Transitivity<T>::post(const int idx) {

    Constraint<T>::cons_id = idx;
    Constraint<T>::idempotent = true;

#ifdef DEBUG_CONSTRAINT
  if (debug_flag > 0) {
    std::cout << "post " << *this << std::endl;
  }
#endif

  for (size_t i{0}; i < the_tasks.size(); ++i) {
    m_solver.wake_me_on(lb<T>(the_tasks[i].getStart()), this->id());
    m_solver.wake_me_on(ub<T>(the_tasks[i].getEnd()), this->id());

    for (size_t j{0}; j < the_tasks.size(); ++j)
      if (i != j) {

        m_solver.wake_me_on(disjunct[i][j], this->id());
        scopex.push_back(i);
        scopey.push_back(j);
      }
  }
}

template <typename T>
void Transitivity<T>::add_edge(const int x, const int y, const int r) {

#ifdef DBG_TRANSITIVITY
  if (DBG_TRANSITIVITY)
    std::cout << " ==> add edge t" << the_tasks[x].id() << " -> t"
              << the_tasks[y].id() << " (" << m_solver.pretty(disjunct[x][y]) << ")" << std::endl;
#endif

  if (DAG[x].frontsize() > 0 or DAG[y].backsize() > 0)
    change_flag = true;

#ifdef DBG_LTRANS
  if (m_solver.boolean.isUndefined(disjunct[x][y].variable())) {
    std::cout << "edge pruning\n";
  }
#endif

  m_solver.set(disjunct[x][y], {this, r});

  assert(DAG[x].has(y));
  assert(DAG[y].has(x));

  DAG[x].remove_front(y);
  DAG[y].remove_back(x);

  if (transition_flag) {
    transitive_reduction.remove_back(edge(y, x));
  }
}

template <typename T>
bool Transitivity<T>::notify(const Literal<T> l, const int r) {
  if (l.isNumeric())
    return true;

  auto x{scopex[r]};
  auto y{scopey[r]};

#ifdef DBG_TRANSITIVITY
  if (DBG_TRANSITIVITY) {
    std::cout << std::endl;
    for (size_t i{0}; i < the_tasks.size(); ++i) {
      std::cout << "t" << the_tasks[i].id() << ":";
      for (auto j{DAG[i].fbegin()}; j != DAG[i].fend(); ++j) {
        std::cout << " -> t" << the_tasks[*j].id();
      }
      std::cout << std::endl;
    }

    for (size_t i{0}; i < the_tasks.size(); ++i) {
      std::cout << "t" << the_tasks[i].id() << ":";
      for (auto j{DAG[i].bbegin()}; j != DAG[i].bend(); ++j) {
        std::cout << " <- t" << the_tasks[*j].id();
      }
      std::cout << std::endl;
    }

      std::cout << std::endl;
 
      
      for (size_t i{0}; i < the_tasks.size(); ++i) {
        for (size_t j{0}; j < the_tasks.size(); ++j) if (i != j) {
                std::cout << " " << m_solver.pretty(disjunct[i][j])  ;
            std::cout.flush();
            if(m_solver.boolean.satisfied(disjunct[i][j]))
                std::cout << "[true]";
            if(m_solver.boolean.falsified(disjunct[i][j]))
                std::cout << "[false]";
        }
        std::cout << std::endl;
      }
      
    std::cout << "\nnotify edge " << the_tasks[x] << " -> " << the_tasks[y]
              << " / " << m_solver.pretty(l) << std::endl;
  }
#endif

  if (DAG[x].isfront(y)) {
    assert(DAG[y].isback(x));
    return false;
  }

  new_succ_of_x.clear();
  new_pred_of_y.clear();
  for (auto zp{DAG[y].frbegin()}; zp != DAG[y].frend(); ++zp) {
    auto z{*zp};
    if (not DAG[x].isfront(z)) {
      new_succ_of_x.push_back(z);
    }
  }
  for (auto tp{DAG[x].brbegin()}; tp != DAG[x].brend(); ++tp) {
    auto t{*tp};
    if (not DAG[t].isback(y)) {
      new_pred_of_y.push_back(t);
    }
  }
  for (auto z : new_succ_of_x) {
    for (auto tp{DAG[x].brbegin()}; tp != DAG[x].brend(); ++tp) {
      auto t{*tp};
      if (not DAG[z].isback(t)) {
        add_edge(t, z, r);
      }
    }
    if (not DAG[x].isfront(z))
      add_edge(x, z, r);
  }

  for (auto t : new_pred_of_y) {
    for (auto zp{DAG[y].frbegin()}; zp != DAG[y].frend(); ++zp) {
      auto z{*zp};
      if (not DAG[t].isfront(z)) {
        add_edge(t, z, r);
      }
    }
    if (not DAG[y].isback(t))
      add_edge(t, y, r);
  }

  assert(not DAG[x].isfront(y));

  DAG[x].remove_front(y);
  DAG[y].remove_back(x);

#ifdef DBG_TRANSITIVITY
  if (DBG_TRANSITIVITY) {

      std::cout << std::endl;
      
    for (size_t i{0}; i < the_tasks.size(); ++i) {
      std::cout << "t" << the_tasks[i].id() << ":";
      for (auto j{DAG[i].fbegin()}; j != DAG[i].fend(); ++j) {
        std::cout << " -> t" << the_tasks[*j].id();
      }
      std::cout << std::endl;
    }

    for (size_t i{0}; i < the_tasks.size(); ++i) {
      std::cout << "t" << the_tasks[i].id() << ":";
      for (auto j{DAG[i].bbegin()}; j != DAG[i].bend(); ++j) {
        std::cout << " <- t" << the_tasks[*j].id();
      }
      std::cout << std::endl;
    }

    for (size_t i{0}; i < the_tasks.size(); ++i) {
      for (auto j{DAG[i].fbegin()}; j != DAG[i].fend(); ++j) {
        assert(DAG[*j].isback(i));
      }
    }
  }
#endif

  if (transition_flag) {
    transitive_reduction.remove_back(edge(y, x));

    // for all successors of y
    for (auto zp{DAG[y].frbegin()}; zp != DAG[y].frend(); ++zp) {
      auto e{edge(x, *zp)};
      if (transitive_reduction.has(e)) {
        transitive_reduction.remove_back(e);
      }
    }

    // for all predecessors of x
    for (auto zp{DAG[x].brbegin()}; zp != DAG[x].brend(); ++zp) {
      auto e{edge(*zp, y)};
      if (transitive_reduction.has(e)) {
        transitive_reduction.remove_back(e);
      }
    }

    // for the cross-product:
    for (auto sy{DAG[y].frbegin()}; sy != DAG[y].frend(); ++sy) {
      for (auto px{DAG[x].brbegin()}; px != DAG[x].brend(); ++px) {
        auto e{edge(*px, *sy)};
        if (transitive_reduction.has(e)) {
          transitive_reduction.remove_back(e);
        }
      }
    }

    std::cout << "\nTRANS:\n";
    for (size_t i{0}; i < the_tasks.size(); ++i) {
      std::cout << "t" << the_tasks[i].id() << ":";
      for (auto j{DAG[i].fbegin()}; j != DAG[i].fend(); ++j) {
        std::cout << " -> t" << the_tasks[*j].id();
      }
      std::cout << std::endl;
    }

    std::cout << "TRED:\n";
    for (size_t i{0}; i < the_tasks.size(); ++i) {
      std::cout << "t" << the_tasks[i].id() << ":";
      for (size_t j{0}; j < the_tasks.size(); ++j) {
        if (transitive_reduction.has(edge(i, j)))
          std::cout << " -> t" << the_tasks[j].id();
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

  }

  return true;
}

template <typename T> void Transitivity<T>::min_spanning_tree() {

  std::cout << "\nmin spanning tree\n";

  std::sort(transitive_reduction.begin(), transitive_reduction.end(),
            [&](const int a, const int b) { return length(a) < length(b); });
  transitive_reduction.re_index(transitive_reduction.begin(),
                                transitive_reduction.end());

  T path_length{0};
  size_t path_count{0};
  for (auto e : transitive_reduction) {

    auto sx{forest.find(first(e))};
    auto sy{forest.find(second(e))};

    if (sx != sy) {
      forest.merge_roots(sx, sy);
      path_length += length(e);
      ++path_count;

      std::cout << "add " << first(e) << " -> " << second(e) << " ("
                << length(e) << ") ==> " << path_length << "\n";

      if (path_count == the_tasks.size() - 1) {
        std::cout << "the tree has " << (the_tasks.size() - 1) << " edges\n";
        break;
      }
    } else {
      std::cout << "skip " << first(e) << " -> " << second(e) << "\n";
    }
  }

  T min_start{Constant::Infinity<T>};
  //  for (auto t : m_tasks) {
  //    path_length += m_solver.minDuration(t);
  //    min_start = std::min(min_start, m_solver.lower(START(t)));
  //  }
  for (auto& t : the_tasks) {
    path_length += t.minDuration(m_solver);
    min_start = std::min(min_start, t.getEarliestStart(m_solver));
  }
  path_length += min_start;

    std::cout << "lb = " << path_length << " / " << schedule.getEarliestEnd(m_solver) << ".." << schedule.getLatestEnd(m_solver) << std::endl;

  forest.clear();
}

template <typename T> void Transitivity<T>::propagate() {

#ifdef DBG_LTRANS
  bool pruning{false};
#endif

    /// Update the upper bounds
  for (size_t x{0}; x < the_tasks.size(); ++x) {
    sorted_tasks[x] = x;
    offset[x] = 0;
  }
    // sort by non-decreasing latest end
  std::sort(sorted_tasks.begin(), sorted_tasks.end(),
            [&](const int x, const int y) -> bool {
              return the_tasks[x].getLatestEnd(m_solver) <
                     the_tasks[y].getLatestEnd(m_solver);
            });

#ifdef DBG_TRANSITIVITY
    if (DBG_TRANSITIVITY) {
        std::cout << "\npropagate w.r.t. subsequent tasks\n";
    }
#endif
    
    // starting with the earliest task
  for (auto x : sorted_tasks) {

#ifdef DBG_TRANSITIVITY
    if (DBG_TRANSITIVITY) {

      std::cout << prettyTask(x);
      for (auto yp{DAG[x].bbegin()}; yp != DAG[x].bend(); ++yp) {
        std::cout << " <- t" << the_tasks[*yp].id();
      }
      std::cout << std::endl;
    }
#endif

    // get all the successors y of x in the graph, (so predecessors in the
    // schedule: y < x)
    for (auto yp{DAG[x].bbegin()}; yp != DAG[x].bend(); ++yp) {
      //      offset[*yp] += m_solver.minDuration(m_tasks[x]);
      offset[*yp] += the_tasks[x].minDuration(m_solver);

      // x must end before ex
      //      auto ex{m_solver.upper(END(m_tasks[x]))};
      auto ex{the_tasks[x].getLatestEnd(m_solver)};

      if (ex == Constant::Infinity<T>)
        continue;

      // y must end before ey
      //      auto ey{m_solver.upper(END(m_tasks[*yp]))};
      auto ey{the_tasks[*yp].getLatestEnd(m_solver)};

      // y must end before ex -
      if ((ex - offset[*yp]) < ey) {
#ifdef DBG_LTRANS
        pruning = true;
#endif
          
          
#ifdef DBG_TRANSITIVITY
        if (DBG_TRANSITIVITY) {
          std::cout << " new bound (" << the_tasks[*yp] << " must end before "
                    << ex << " - " << offset[*yp] << ") "
                    << the_tasks[*yp].end.before(ex - offset[*yp]) << "/" << ey
                    << std::endl;
        }
#endif
        auto bc{the_tasks[*yp].end.before(ex - offset[*yp])};
        m_solver.set(bc, {this, ex});
      }
    }
  }

    
    /// Update the lower bounds
#ifdef DBG_TRANSITIVITY
    if (DBG_TRANSITIVITY) {
        std::cout << "\npropagate w.r.t. preceding tasks\n";
    }
#endif
  for (size_t x{0}; x < the_tasks.size(); ++x) {
    offset[x] = 0;
  }
    // sort by non-increasing earliest start
  std::sort(sorted_tasks.begin(), sorted_tasks.end(),
            [&](const int x, const int y) -> bool {
              return the_tasks[x].getEarliestStart(m_solver) >
                     the_tasks[y].getEarliestStart(m_solver);
            });

// starting with the latest task
  for (auto y : sorted_tasks) {
#ifdef DBG_TRANSITIVITY
    if (DBG_TRANSITIVITY) {
      std::cout << prettyTask(y);
      for (auto xp{DAG[y].fbegin()}; xp != DAG[y].fend(); ++xp) {
        std::cout << " -> t" << the_tasks[*xp].id();
      }
    }
    std::cout << std::endl;
#endif

      // the successors of y are pushed by y's duration
    for (auto xp{DAG[y].fbegin()}; xp != DAG[y].fend(); ++xp) {
      offset[*xp] += the_tasks[y].minDuration(m_solver);

      auto sy{the_tasks[y].getEarliestStart(m_solver)};
      auto sz{the_tasks[*xp].getEarliestStart(m_solver)};

      if (sy == -Constant::Infinity<T>)
        continue;

      if ((sy + offset[*xp]) > sz) {
#ifdef DBG_LTRANS
        pruning = true;
#endif

#ifdef DBG_TRANSITIVITY
      if (DBG_TRANSITIVITY)
        std::cout << " new bound (" << the_tasks[*xp] << " must start after "
                  << sy << " + " << offset[*xp] << ") "
                  << the_tasks[*xp].start.after(sy + offset[*xp]) << "/" << sz
                  << std::endl;
#endif

      m_solver.set(the_tasks[*xp].start.after(sy + offset[*xp]),
                   {this, sy});
      }
    }
  }

  if (transition_flag) {

    assert(false);
    min_spanning_tree();
  }

#ifdef DBG_LTRANS
  if (pruning)
    std::cout << "bound pruning\n";
#endif
}

template <typename T> int Transitivity<T>::getType() const {
  return TRANSITIVITYEXPL;
}

template <typename T>
void Transitivity<T>::xplain(const Literal<T> l, const hint h,
                             std::vector<Literal<T>> &Cl) {

    auto l_lvl{m_solver.propagationLevel(l)};
    
  if (l.isNumeric()) {

#ifdef DBG_EXPL_TRANS
    std::cout << "explain " << l << " with transitivity constraint (hint=" << h
              << ")\n";
#endif

    
#ifdef DBG_EXPL_TRANS
    std::cout << " reason was: the "
              << (l.sign() == bound::lower ? "predecessors" : "successors")
              << " of " << the_tasks[task_map[l.variable()]] << " that were "
              << (l.sign() == bound::lower ? "above " : "below ") << h
              << " @lvl " << l_lvl << std::endl;
#endif
      
    auto x{task_map[l.variable()]};
    if (l.sign() == bound::lower) {
#ifdef DBG_EXPL_TRANS
      std::cout << "lb of " << the_tasks[x]
                << " dur=" << the_tasks[x].minDuration(m_solver) << std::endl;
#endif
      for (auto yp{DAG[x].bbegin()}; yp != DAG[x].bend(); ++yp) {
        auto p{disjunct[*yp][x]};
          auto b{the_tasks[*yp].start.after(h)};
        if (m_solver.propagationLevel(p) < l_lvl and m_solver.numeric.satisfied(b)) {
          
            if(m_solver.propagationLevel(b) < l_lvl) {
                Cl.push_back(p);
                Cl.push_back(b);
                
#ifdef DBG_EXPL_TRANS
                std::cout << " - " << m_solver.pretty(p) << " & " << b << " ("
                          << the_tasks[*yp].minDuration(m_solver) << ") @"
                          << m_solver.propagationLevel(p) << " & "
                          << m_solver.propagationLevel(b) << std::endl;
#endif
            }
        }
      }
    } else {
#ifdef DBG_EXPL_TRANS
      std::cout << "ub of " << the_tasks[x]
                << " dur=" << the_tasks[x].minDuration(m_solver) << std::endl;
#endif
      for (auto yp{DAG[x].fbegin()}; yp != DAG[x].fend(); ++yp) {
        auto p{disjunct[x][*yp]};
          auto b{the_tasks[*yp].end.before(h)};
        if (m_solver.propagationLevel(p) < l_lvl and m_solver.numeric.satisfied(b)) {
          
            if(m_solver.propagationLevel(b) < l_lvl) {
                Cl.push_back(p);
                Cl.push_back(b);
                
                
                //            assert(m_solver.propagationLevel(b) < l_lvl);
                
#ifdef DBG_EXPL_TRANS
                std::cout << " - " << m_solver.pretty(p) << " & " << b << " ("
                          << the_tasks[*yp].minDuration(m_solver) << ") @"
                          << m_solver.propagationLevel(p) << " & "
                          << m_solver.propagationLevel(b) << std::endl;
#endif
            }
        }
      }
    }

  } else {

#ifdef DBG_EXPL_TRANS
    std::cout << "explain " << m_solver.pretty(l)
              << " with transitivity constraint (hint=" << h << ")\n";
#endif

    int i{scopex[h]};
    int j{scopey[h]};
    auto r{disjunct[i][j]};

#ifdef DBG_EXPL_TRANS
    std::cout << " reason was d[" << i << "][" << j
              << "] = " << m_solver.pretty(r) << std::endl;
#endif

    auto lc{m_solver.boolean.getEdge(l)};
    auto rc{m_solver.boolean.getEdge(r)};
      Cl.push_back(r);
      
      assert(m_solver.propagationLevel(r) < l_lvl);

#ifdef DBG_EXPL_TRANS
    std::cout << rc.from << " -> " << rc.to;
#endif

    if (lc.from != rc.from) {

      auto x{task_map[lc.from]};
      auto y{task_map[rc.from]};
      Cl.push_back(disjunct[y][x]);
        
        assert(m_solver.propagationLevel(disjunct[y][x]) < l_lvl);

#ifdef DBG_EXPL_TRANS
      std::cout << " & " << lc.from << " -> " << rc.from << " (" << the_tasks[x]
                << " -> " << the_tasks[y] << "/"
                << m_solver.pretty(disjunct[y][x]) << ")";
#endif
    }

    if (lc.to != rc.to) {

      auto x{task_map[lc.to]};
      auto y{task_map[rc.to]};
      Cl.push_back(disjunct[x][y]);
        
        assert(m_solver.propagationLevel(disjunct[x][y]) < l_lvl);

#ifdef DBG_EXPL_TRANS
      std::cout << " & " << rc.to << " -> " << lc.to << " (" << the_tasks[x]
                << " -> " << the_tasks[y] << "/"
                << m_solver.pretty(disjunct[x][y]) << ")";
#endif
    }

#ifdef DBG_EXPL_TRANS
    std::cout << " ===> " << lc.from << " -> " << lc.to << std::endl;
#endif

    //        if(lc.to != rc.to)
    //        exit(1);
  }
}

template <typename T>
std::ostream &Transitivity<T>::display(std::ostream &os) const {
  os << "Transitivity";

#ifdef DEBUG_CONSTRAINT
  os << "[" << cons_id << "]";
#endif

  os << "(";
  for (auto &t : the_tasks) {
    std::cout << " t" << t.id();
  }
  std::cout << " )";
  return os;
}

template <typename T>
std::ostream &Transitivity<T>::print_reason(std::ostream &os,
                                            const hint) const {
  //  display(os);
  os << "transitivity";
  //
  //  if (not explanations[h].empty()) {
  //
  //    auto l{explanations[h].begin()};
  //    m_solver.displayLiteral(os, *l);
  //    ++l;
  //    while (l != explanations[h].end()) {
  //      os << ", ";
  //      m_solver.displayLiteral(os, *l);
  //      ++l;
  //    }
  //  }
  //
  //  os << ")";
  return os;
}

template <typename T>
std::string Transitivity<T>::prettyTask(const int i) const {
  std::stringstream ss;
  //  ss << "t" << m_tasks[i] << ": [" << est(i) << ".." << lct(i) << "] ("
  //     << minduration(i) << ")";
  ss << "t" << the_tasks[i].id() << ": ["
     << the_tasks[i].getEarliestStart(m_solver) << ".."
     << the_tasks[i].getLatestEnd(m_solver) << "] ("
     << the_tasks[i].minDuration(m_solver) << ")";
  return ss.str();
}

// template <typename T> std::vector<int> Transitivity<T>::task_map;

} // namespace tempo

#endif
