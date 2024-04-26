#ifndef TEMPO_TRANSITIVITY_HPP
#define TEMPO_TRANSITIVITY_HPP

#include <cassert>
#include <map>
#include <vector>

#include "Explanation.hpp"
#include "Global.hpp"
#include "ReversibleObject.hpp"
#include "Scheduler.hpp"
#include "constraints/Constraint.hpp"
#include "util/DisjointSet.hpp"
#include "util/SparseSet.hpp"

//#define DBG_LTRANS

namespace tempo {

template <typename T> class Transitivity : public Constraint {
private:
  Scheduler<T> &m_schedule;
  std::vector<task> m_tasks;
  std::vector<int> task_map;

  //  std::vector<SparseSet<int, Reversible<size_t>>> forward;
  //  std::vector<SparseSet<int, Reversible<size_t>>> backward;

  // encoding: (y \in front(DAG[x]) && x \in back(y)) <=> edge (x,y)
  std::vector<SparseSet<int, Reversible<size_t>>> DAG;
  // encoding: y \in RED[x] <=> edge (x,y)
//  std::vector<SparseSet<int, Reversible<size_t>>> RED;
  SparseSet<int, Reversible<size_t>> transitive_reduction;
  DisjointSet<> forest;
  std::vector<std::vector<T>> distance_matrix;

  //    std::vector<var>
  std::vector<std::vector<lit>> disjunct;
  std::vector<int> scopex;
  std::vector<int> scopey;

  std::vector<int> sorted_tasks;
  std::vector<T> offset;

  //  std::vector<int> from;
  //  std::vector<int> to;

  //  SparseSet<> changed_pred;
  //  SparseSet<> changed_succ;

  std::vector<int> new_succ_of_x;
  std::vector<int> new_pred_of_y;

  bool change_flag{false};

  bool transition_flag{false};

  //  std::vector<lit> cur_explanation;
  //  std::vector<std::vector<lit>> explanations;

  //    int lvl;

  T length(const size_t e) const {
      auto i{first(e)};
      auto j{second(e)};
      return transition_time(i,j);
  }
    size_t edge(const size_t i, const size_t j) const {
      return m_tasks.size() * i + j;
    }
  size_t first(const size_t e) const { return e / m_tasks.size(); }
  size_t second(const size_t e) const { return e % m_tasks.size(); }
  T transition_time(const int i, const int j) const;

public:
  template <typename ItTask, typename ItVar>
  Transitivity(Scheduler<T> &scheduler, const ItTask beg_task,
               const ItTask end_task, const ItVar beg_var, const ItVar end_var);
  virtual ~Transitivity();

  void add_edge(const int x, const int y, const int r);

  bool notify_bound(const int lit, const int rank) override;
  bool notify_edge(const int lit, const int rank) override;
  void post(const int idx) override;
  void propagate() override;
  void min_spanning_tree();

  void xplain(const lit l, const hint h, std::vector<lit> &Cl) override;
  int getType() const override;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

#ifdef DEBUG_CONSTRAINT
  int debug_flag{2};
#endif
};

template <typename T>
template <typename ItTask, typename ItVar>
Transitivity<T>::Transitivity(Scheduler<T> &scheduler, const ItTask beg_task,
                              const ItTask end_task, const ItVar beg_var,
                              const ItVar end_var)
    : m_schedule(scheduler),
      transitive_reduction(std::distance(beg_task, end_task) *
                               std::distance(beg_task, end_task),
                           &m_schedule.getEnv())
//    forward(Reversible<size_t>(0, m_schedule.getEnv()), Reversible<size_t>(0,
//    m_schedule.getEnv())),
// forward(Reversible<size_t>(0, m_schedule.getEnv()), Reversible<size_t>(0,
// m_schedule.getEnv())),
{

  priority = LOW;

  task_map.resize(m_schedule.numTask());

  // get all tasks with non-zero duration
  size_t i{0}, j{0};
  for (auto jp{beg_task}; jp != end_task; ++jp) {

    task t{*jp};
    task_map[t] = i++;
    m_tasks.push_back(t);
  }

  //        forward.resize(m_tasks.size());
  //        for(auto &row : forward) {
  //            row.reserve(m_tasks.size());
  //        }
  //
  //        backward.resize(m_tasks.size());
  //        for(auto &col : backward) {
  //            col.reserve(m_tasks.size());
  //        }

  //  changed.reserve(m_tasks.size());

  disjunct.resize(m_tasks.size());
  sorted_tasks.resize(m_tasks.size());
  offset.resize(m_tasks.size());
  //  changed_pred.reserve(m_tasks.size());
  //  changed_succ.reserve(m_tasks.size());

  transitive_reduction.fill(); // resize(m_tasks.size() * m_tasks.size());
  forest.resize(m_tasks.size());

  for (i = 0; i < m_tasks.size(); ++i) {

    disjunct[i].resize(m_tasks.size());

    DAG.emplace_back(m_tasks.size(), &m_schedule.getEnv());
    DAG.back().fill();

//    RED.emplace_back(m_tasks.size(), &m_schedule.getEnv());
//    RED.back().fill();
//    RED.back().remove_back(i);

    transitive_reduction.remove_back(edge(i, i));
    //        SparseSet<int, Reversible<size_t>> row(Reversible<size_t>(0,
    //        &m_schedule.getEnv()), Reversible<size_t>(0,
    //        &m_schedule.getEnv()), m_tasks.size()); forward.push_back(row);
    //    forward.emplace_back(m_tasks.size(), &m_schedule.getEnv());
    //    backward.emplace_back(m_tasks.size(), &m_schedule.getEnv());
    //        backward.emplace_back(Reversible<size_t>(0, &m_schedule.getEnv()),
    //        Reversible<size_t>(0, &m_schedule.getEnv()), m_tasks.size());
  }

  //    int i,j;
  for (auto v{beg_var}; v != end_var; ++v) {
    auto ef{m_schedule.getEdge(POS(*v))};
    i = task_map[TASK(ef.from)];
    j = task_map[TASK(ef.to)];
    disjunct[i][j] = POS(*v);

    if (not transition_flag) {
      transition_flag = (transition_time(i, j) > 0);
    }

    auto eb{m_schedule.getEdge(NEG(*v))};
    i = task_map[TASK(eb.from)];
    j = task_map[TASK(eb.to)];
    disjunct[i][j] = NEG(*v);

    if (not transition_flag) {
      transition_flag = (transition_time(i, j) > 0);
    }
  }

  //    distance_matrix.resize(m_tasks.size());
  //    for(size_t i{0}; i<m_tasks.size(); ++i) {
  //        distance_matrix[i].resize(m_tasks.size());
  //        for(size_t j{0}; j<m_tasks.size(); ++j) if(i != j) {
  //            distance_matrix[i][j] =
  //            -m_schedule.getEdge(disjunct[i][j]).distance;
  //        }
  //    }
}

template <typename T> Transitivity<T>::~Transitivity() {}

template <typename T>
T Transitivity<T>::transition_time(const int i, const int j) const {
  return -m_schedule.getEdge(disjunct[i][j]).distance;
}

template <typename T> void Transitivity<T>::post(const int idx) {

  cons_id = idx;
  idempotent = true;

#ifdef DEBUG_CONSTRAINT
  if (debug_flag > 0) {
    std::cout << "post " << *this << std::endl;
  }
#endif

  //  for (auto e : edges) {
  //  for (size_t i{0}; i < m_tasks.size(); ++i) {
  //    for (size_t j{i + 1}; j < m_tasks.size(); ++j) {
  //      lit l{m_schedule.getEdgeLit({START(m_tasks[i]), END(m_tasks[j])})};
  //      m_schedule.wake_me_on_edge(l, cons_id);
  //      m_schedule.wake_me_on_edge(NOT(l), cons_id);
  //    }
  //  }

  for (size_t i{0}; i < m_tasks.size(); ++i) {
    m_schedule.wake_me_on_event(LOWERBOUND(START(m_tasks[i])), cons_id);
    m_schedule.wake_me_on_event(UPPERBOUND(END(m_tasks[i])), cons_id);

    for (size_t j{0}; j < m_tasks.size(); ++j)
      if (i != j) {
        //        lit l{m_schedule.getEdgeLit({START(m_tasks[i]),
        //        END(m_tasks[j])})};
        m_schedule.wake_me_on_edge(disjunct[i][j], cons_id);
        scopex.push_back(i);
        scopey.push_back(j);
        //        m_schedule.wake_me_on_edge(disjunct[j][i], cons_id);
      }
  }

  //    m_schedule.wake_me_on_edge(POS(e), cons_id);
  //    m_schedule.wake_me_on_edge(NEG(e), cons_id);
  //  }
}

template <typename T>
void Transitivity<T>::add_edge(const int x, const int y, const int r) {

#ifdef DBG_TRANSITIVITY
  if (DBG_TRANSITIVITY)
    std::cout << " ==> add edge t" << m_tasks[x] << " -> t" << m_tasks[y]
              << std::endl;
#endif

  if (DAG[x].frontsize() > 0 or DAG[y].backsize() > 0)
    change_flag = true;

  //            << m_schedule.prettyLiteral(EDGE(disjunct[x][y])) << std::endl;

#ifdef DBG_LTRANS
  if (m_schedule.isUndefined(VAR(disjunct[x][y]))) {
    std::cout << "edge pruning\n";
  }
#endif

  m_schedule.set(disjunct[x][y], {this, r});
  //    }

  //    else std::cout << "already here\n";

  assert(DAG[x].has(y));
  assert(DAG[y].has(x));

  DAG[x].remove_front(y);
  DAG[y].remove_back(x);

  if (transition_flag) {
    //    RED[y].remove_back(x);
    transitive_reduction.remove_back(edge(y, x));
  }
}

template <typename T> bool Transitivity<T>::notify_bound(const lit, const int) {
  return true;
}

template <typename T>
bool Transitivity<T>::notify_edge(const lit l, const int r) {
  //  change_flag = false;
  auto e{m_schedule.getEdge(l)};

#ifdef DBG_TRANSITIVITY
  if (DBG_TRANSITIVITY) {
    std::cout << std::endl;
    //      std::cout << "before notify\nnew succs:";
    //      for(auto x : changed_succ)
    //          std::cout << " t" << m_tasks[x];
    //      std::cout << std::endl;
    //      std::cout << "new preds:";
    //      for(auto x : changed_pred)
    //          std::cout << " t" << m_tasks[x];
    //      std::cout << std::endl;
    for (size_t i{0}; i < m_tasks.size(); ++i) {
      std::cout << "t" << m_tasks[i] << ":";
      for (auto j{DAG[i].fbegin()}; j != DAG[i].fend(); ++j) {
        std::cout << " -> t" << m_tasks[*j];
      }
      std::cout << std::endl;
    }

    for (size_t i{0}; i < m_tasks.size(); ++i) {
      std::cout << "t" << m_tasks[i] << ":";
      for (auto j{DAG[i].bbegin()}; j != DAG[i].bend(); ++j) {
        std::cout << " <- t" << m_tasks[*j];
      }
      std::cout << std::endl;
    }

    std::cout << "notify edge t" << TASK(e.from) << " -> t" << TASK(e.to)
              << std::endl;
  }
#endif

  auto x{task_map[TASK(e.from)]};
  auto y{task_map[TASK(e.to)]};

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

    for (size_t i{0}; i < m_tasks.size(); ++i) {
      std::cout << "t" << m_tasks[i] << ":";
      for (auto j{DAG[i].fbegin()}; j != DAG[i].fend(); ++j) {
        std::cout << " -> t" << m_tasks[*j];
      }
      std::cout << std::endl;
    }

    for (size_t i{0}; i < m_tasks.size(); ++i) {
      std::cout << "t" << m_tasks[i] << ":";
      for (auto j{DAG[i].bbegin()}; j != DAG[i].bend(); ++j) {
        std::cout << " <- t" << m_tasks[*j];
      }
      std::cout << std::endl;
    }

    for (size_t i{0}; i < m_tasks.size(); ++i) {
      for (auto j{DAG[i].fbegin()}; j != DAG[i].fend(); ++j) {
        assert(DAG[*j].isback(i));
      }
    }
  }
#endif

  if (transition_flag) {
    //    std::cout << "rm* " << y << " -> " << x << std::endl;
    //    RED[y].remove_back(x);
    transitive_reduction.remove_back(edge(y, x));

    // new edge (x,y)

    // for all successors of y
    for (auto zp{DAG[y].frbegin()}; zp != DAG[y].frend(); ++zp) {
      // z cannot be a direct successor of x
      //      if (RED[x].has(*zp)) {
      //        RED[x].remove_back(*zp);
      //      }
      auto e{edge(x, *zp)};
      if (transitive_reduction.has(e)) {
        transitive_reduction.remove_back(e);
      }
    }

    // for all predecessors of x
    for (auto zp{DAG[x].brbegin()}; zp != DAG[x].brend(); ++zp) {
      // z cannot be a direct successor of x
      //      if (RED[*zp].has(y)) {
      //        RED[*zp].remove_back(y);
      //      }
      auto e{edge(*zp, y)};
      if (transitive_reduction.has(e)) {
        transitive_reduction.remove_back(e);
      }
    }

    // for the cross-product:
    for (auto sy{DAG[y].frbegin()}; sy != DAG[y].frend(); ++sy) {
      for (auto px{DAG[x].brbegin()}; px != DAG[x].brend(); ++px) {
        //        if (RED[*px].has(*sy)) {
        //          RED[*px].remove_back(*sy);
        //        }
        auto e{edge(*px, *sy)};
        if (transitive_reduction.has(e)) {
          transitive_reduction.remove_back(e);
        }
      }
    }

    std::cout << "\nTRANS:\n";
    for (size_t i{0}; i < m_tasks.size(); ++i) {
      std::cout << "t" << m_tasks[i] << ":";
      for (auto j{DAG[i].fbegin()}; j != DAG[i].fend(); ++j) {
        std::cout << " -> t" << m_tasks[*j];
      }
      std::cout << std::endl;
    }
    //
    //    std::cout << "RED:\n";
    //    for (size_t i{0}; i < m_tasks.size(); ++i) {
    //      std::cout << "t" << m_tasks[i] << ":";
    //      for (auto j : RED[i]) {
    //        std::cout << " -> t" << m_tasks[j];
    //      }
    //      std::cout << std::endl;
    //    }
    //    std::cout << std::endl;

    std::cout << "TRED:\n";
    for (size_t i{0}; i < m_tasks.size(); ++i) {
      std::cout << "t" << m_tasks[i] << ":";
      for (size_t j{0}; j < m_tasks.size(); ++j) {
        if (transitive_reduction.has(edge(i, j)))
          std::cout << " -> t" << m_tasks[j];
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    //    for (size_t i{0}; i < m_tasks.size(); ++i) {
    //      for (size_t j{0}; j < m_tasks.size(); ++j) {
    //        if (RED[i].has(j) != transitive_reduction.has(edge(i, j))) {
    //          std::cout << "ERROR (" << m_tasks[i] << "," << m_tasks[j] <<
    //          ")\n"; exit(1);
    //        }
    //      }
    //    }
  }

  return true;
}

// template <typename T>
// bool Transitivity<T>::udpate_transitive_reduction(const lit l, const int r) {
//   //  change_flag = false;
//   auto e{m_schedule.getEdge(l)};
//
//#ifdef DBG_TRANSITIVITY
//   if (DBG_TRANSITIVITY) {
//     std::cout << std::endl;
//     //      std::cout << "before notify\nnew succs:";
//     //      for(auto x : changed_succ)
//     //          std::cout << " t" << m_tasks[x];
//     //      std::cout << std::endl;
//     //      std::cout << "new preds:";
//     //      for(auto x : changed_pred)
//     //          std::cout << " t" << m_tasks[x];
//     //      std::cout << std::endl;
//     for (size_t i{0}; i < m_tasks.size(); ++i) {
//       std::cout << "t" << m_tasks[i] << ":";
//       for (auto j{RED[i].fbegin()}; j != RED[i].fend(); ++j) {
//         std::cout << " -> t" << m_tasks[*j];
//       }
//       std::cout << std::endl;
//     }
//
//     for (size_t i{0}; i < m_tasks.size(); ++i) {
//       std::cout << "t" << m_tasks[i] << ":";
//       for (auto j{RED[i].bbegin()}; j != RED[i].bend(); ++j) {
//         std::cout << " <- t" << m_tasks[*j];
//       }
//       std::cout << std::endl;
//     }
//
//     std::cout << "notify edge t" << TASK(e.from) << " -> t" << TASK(e.to)
//               << std::endl;
//   }
//#endif
//
//   auto x{task_map[TASK(e.from)]};
//   auto y{task_map[TASK(e.to)]};
//
//   if (DAG[x].isfront(y)) {
//     assert(DAG[y].isback(x));
//     return false;
//   }
//
//   new_succ_of_x.clear();
//   new_pred_of_y.clear();
//   for (auto zp{DAG[y].frbegin()}; zp != DAG[y].frend(); ++zp) {
//     auto z{*zp};
//     if (not DAG[x].isfront(z)) {
//       new_succ_of_x.push_back(z);
//     }
//   }
//   for (auto tp{DAG[x].brbegin()}; tp != DAG[x].brend(); ++tp) {
//     auto t{*tp};
//     if (not DAG[t].isback(y)) {
//       new_pred_of_y.push_back(t);
//     }
//   }
//   for (auto z : new_succ_of_x) {
//     for (auto tp{DAG[x].brbegin()}; tp != DAG[x].brend(); ++tp) {
//       auto t{*tp};
//       if (not DAG[z].isback(t)) {
//         add_edge(t, z, r);
//       }
//     }
//     if (not DAG[x].isfront(z))
//       add_edge(x, z, r);
//   }
//
//   for (auto t : new_pred_of_y) {
//     for (auto zp{DAG[y].frbegin()}; zp != DAG[y].frend(); ++zp) {
//       auto z{*zp};
//       if (not DAG[t].isfront(z)) {
//         add_edge(t, z, r);
//       }
//     }
//     if (not DAG[y].isback(t))
//       add_edge(t, y, r);
//   }
//
//   assert(not DAG[x].isfront(y));
//
//   DAG[x].remove_front(y);
//   DAG[y].remove_back(x);
//
//#ifdef DBG_TRANSITIVITY
//   if (DBG_TRANSITIVITY) {
//
//     for (size_t i{0}; i < m_tasks.size(); ++i) {
//       std::cout << "t" << m_tasks[i] << ":";
//       for (auto j{DAG[i].fbegin()}; j != DAG[i].fend(); ++j) {
//         std::cout << " -> t" << m_tasks[*j];
//       }
//       std::cout << std::endl;
//     }
//
//     for (size_t i{0}; i < m_tasks.size(); ++i) {
//       std::cout << "t" << m_tasks[i] << ":";
//       for (auto j{DAG[i].bbegin()}; j != DAG[i].bend(); ++j) {
//         std::cout << " <- t" << m_tasks[*j];
//       }
//       std::cout << std::endl;
//     }
//
//     for (size_t i{0}; i < m_tasks.size(); ++i) {
//       for (auto j{DAG[i].fbegin()}; j != DAG[i].fend(); ++j) {
//         assert(DAG[*j].isback(i));
//       }
//     }
//   }
//#endif
//
//   return true;
// }

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

      if (path_count == m_tasks.size() - 1) {
        std::cout << "the tree has " << (m_tasks.size() - 1) << " edges\n";
        break;
      }
    } else {
      std::cout << "skip " << first(e) << " -> " << second(e) << "\n";
    }
  }

  T min_start{INFTY};
  for (auto t : m_tasks) {
    path_length += m_schedule.minDuration(t);
    min_start = std::min(min_start, m_schedule.lower(START(t)));
  }
  path_length += min_start;

  std::cout << "lb = " << path_length << " / " << m_schedule.lower(HORIZON)
            << ".." << m_schedule.upper(HORIZON) << std::endl;

  forest.clear();
}

template <typename T> void Transitivity<T>::propagate() {

#ifdef DBG_LTRANS
  bool pruning{false};
#endif

  for (size_t x{0}; x < m_tasks.size(); ++x) {
    sorted_tasks[x] = x;
    offset[x] = 0;
  }
  std::sort(sorted_tasks.begin(), sorted_tasks.end(),
            [&](const int x, const int y) -> bool {
              return m_schedule.upper(END(m_tasks[x])) <
                     m_schedule.upper(END(m_tasks[y]));
            });

  for (auto x : sorted_tasks) {

#ifdef DBG_TRANSITIVITY
    if (DBG_TRANSITIVITY) {
      std::cout << "t" << m_tasks[x] << " ["
                << m_schedule.lower(START(m_tasks[x])) << ".."
                << m_schedule.upper(END(m_tasks[x])) << "] ("
                << m_schedule.minDuration(m_tasks[x]) << ")";
      for (auto yp{DAG[x].fbegin()}; yp != DAG[x].fend(); ++yp) {
        std::cout << " -> t" << m_tasks[*yp];
      }
      std::cout << std::endl;
    }
#endif

    // get all the successors of x in the graph, (so predecessors in the
    // schedule)
    for (auto yp{DAG[x].fbegin()}; yp != DAG[x].fend(); ++yp) {
      offset[*yp] += m_schedule.minDuration(m_tasks[x]);

      // x must end before ex
      auto ex{m_schedule.upper(END(m_tasks[x]))};

      auto ey{m_schedule.upper(END(m_tasks[*yp]))};

#ifdef DBG_TRANSITIVITY
      if (DBG_TRANSITIVITY)
        std::cout << " new bound " << prettyEvent(END(m_tasks[*yp]))
                  << " <= " << ex - offset[*yp] << "/" << ey << std::endl;
#endif

      if ((ex - offset[*yp]) < ey) {
#ifdef DBG_LTRANS
        pruning = true;
#endif
        m_schedule.set({UPPERBOUND(END(m_tasks[*yp])), ex - offset[*yp]},
                       {this, offset[*yp]});
      }
    }
//        }
#ifdef DBG_TRANSITIVITY
    if (DBG_TRANSITIVITY)
      std::cout << std::endl;
#endif
  }

#ifdef DBG_TRANSITIVITY
  if (DBG_TRANSITIVITY)
    std::cout << std::endl;
#endif

  for (size_t x{0}; x < m_tasks.size(); ++x) {
    //        sorted_tasks[x] = x;
    offset[x] = 0;
  }
  std::sort(sorted_tasks.begin(), sorted_tasks.end(),
            [&](const int x, const int y) -> bool {
              return m_schedule.lower(START(m_tasks[x])) >
                     m_schedule.lower(START(m_tasks[y]));
            });

  //    for(auto x : sorted_tasks) {
  //        std::cout << "[" << m_schedule.lower(START(m_tasks[x])) << ".." <<
  //        m_schedule.upper(END(m_tasks[x])) << "]\n";
  //    }

  for (auto y : sorted_tasks) {
#ifdef DBG_TRANSITIVITY
    if (DBG_TRANSITIVITY) {
      std::cout << "t" << m_tasks[y] << " ["
                << m_schedule.lower(START(m_tasks[y])) << ".."
                << m_schedule.upper(END(m_tasks[y])) << "] ("
                << m_schedule.minDuration(m_tasks[y]) << ")";
      //        if(DAG[x].frontsize() > 0) {
      for (auto xp{DAG[y].bbegin()}; xp != DAG[y].bend(); ++xp) {
        std::cout << " <- t" << m_tasks[*xp];
      }
    }
    std::cout << std::endl;
#endif
    //        if(DAG[x].frontsize() > 0) {
    for (auto xp{DAG[y].bbegin()}; xp != DAG[y].bend(); ++xp) {
      offset[*xp] += m_schedule.minDuration(m_tasks[y]);

      auto sy{m_schedule.lower(START(m_tasks[y]))};
      auto sz{m_schedule.lower(START(m_tasks[*xp]))};

#ifdef DBG_TRANSITIVITY
      if (DBG_TRANSITIVITY)
        std::cout << " new bound " << prettyEvent(START(m_tasks[*xp]))
                  << " >= " << sy + offset[*xp] << "/" << sz << std::endl;
#endif

      if ((sy + offset[*xp]) > sz) {
#ifdef DBG_LTRANS
        pruning = true;
#endif
        m_schedule.set({LOWERBOUND(START(m_tasks[*xp])), -sy - offset[*xp]},
                       {this, offset[*xp]});
      }
    }
//        }
#ifdef DBG_TRANSITIVITY
    if (DBG_TRANSITIVITY)
      std::cout << std::endl;
#endif
  }

#ifdef DBG_LTRANS
  if (pruning)
    std::cout << "bound pruning\n";
#endif

  if (transition_flag)
    min_spanning_tree();
}

template <typename T> int Transitivity<T>::getType() const {
  return TRANSITIVITYEXPL;
}

template <typename T>
void Transitivity<T>::xplain(const lit l, const hint h, std::vector<lit> &Cl) {

#ifdef DBG_EXPL_TRANS
  std::cout << "explain " << m_schedule.prettyLiteral(l)
            << " with transitivity constraint (h=" << h << ")\n";
#endif

  if (LTYPE(l) == EDGE_LIT) {

    //        int n{static_cast<int>(m_tasks.size())};
    //        int i{h/n};
    //        int j{h%n};
    int i{scopex[h]};
    int j{scopey[h]};
    auto r{disjunct[i][j]};

#ifdef DBG_EXPL_TRANS
    std::cout << " reason was d[" << i << "][" << j << "] = ";
    std::cout << m_schedule.prettyLiteral(EDGE(r));
#endif

    Cl.push_back(EDGE(r));
    auto el{m_schedule.getEdge(FROM_GEN(l))};
    auto er{m_schedule.getEdge(r)};
    if (el.from != er.from) {

      //            std::cout << " {" << prettyEvent(el.from) << "->" <<
      //            prettyEvent(er.from) << "}"; std::cout.flush();

      lit p{m_schedule.getEdgeLit({el.from, NOT(er.from)})};

#ifdef DBG_EXPL_TRANS
      std::cout << " and " << m_schedule.prettyLiteral(EDGE(p));
      std::cout.flush();
#endif

      Cl.push_back(EDGE(p));
    }
    if (el.to != er.to) {

      //            std::cout << " {" << prettyEvent(er.to) << "->" <<
      //            prettyEvent(el.to) << "}"; std::cout.flush();

      lit p{m_schedule.getEdgeLit({NOT(er.to), el.to})};

#ifdef DBG_EXPL_TRANS
      std::cout << " and " << m_schedule.prettyLiteral(EDGE(p));
      std::cout.flush();
#endif

      Cl.push_back(EDGE(p));
    }

#ifdef DBG_EXPL_TRANS
    std::cout << std::endl;
#endif

  } else {

    auto bc{m_schedule.getBound(FROM_GEN(l))};

#ifdef DBG_EXPL_TRANS
    std::cout << " reason was all predecessors within " << h << " of "
              << prettyEvent(EVENT(bc.l)) << std::endl;
#endif

    if (SIGN(bc.l) == LOWER) {
      auto t{TASK(EVENT(bc.l))};
      task x{-1};
      for (unsigned i{0}; i < m_tasks.size(); ++i) {
        if (m_tasks[i] == t) {
          x = i;
          break;
        }
      }

      for (auto yp{DAG[x].fbegin()}; yp != DAG[x].fend(); ++yp) {
        BoundConstraint<T> yc{LIT(START(m_tasks[*yp]), LOWER), bc.distance + h};

#ifdef DBG_EXPL_TRANS
        std::cout << " implicant of " << yc << ": ";
#endif

        auto p{m_schedule.getImplicant(yc)};

        if (p != NoLit) {
          if (p < FROM_GEN(l) and
              m_schedule.getBound(p).distance <= yc.distance) {

#ifdef DBG_EXPL_TRANS
            std::cout << m_schedule.prettyLiteral(BOUND(p)) << " ("
                      << m_schedule.minDuration(
                             TASK(EVENT(m_schedule.getBound(p).l)))
                      << ") & "
                      << m_schedule.prettyLiteral(EDGE(disjunct[x][*yp]))
                      << std::endl;
#endif

            assert(m_schedule.getBound(p).distance <= yc.distance);
            Cl.push_back(BOUND(p));
            Cl.push_back(EDGE(disjunct[x][*yp]));
          }
#ifdef DBG_EXPL_TRANS
          else {
            std::cout << "none\n";
          }
#endif
        }
      }
    } else {
      auto t{TASK(EVENT(bc.l))};
      task y{-1};
      for (unsigned i{0}; i < m_tasks.size(); ++i) {
        if (m_tasks[i] == t) {
          y = i;
          break;
        }
      }
      for (auto xp{DAG[y].bbegin()}; xp != DAG[y].bend(); ++xp) {
        BoundConstraint<T> xc{LIT(END(m_tasks[*xp]), UPPER), bc.distance + h};

#ifdef DBG_EXPL_TRANS
        std::cout << " implicant of " << xc << ": ";
#endif

        auto p{m_schedule.getImplicant(xc)};

        if (p != NoLit) {
          if (p < FROM_GEN(l) and
              m_schedule.getBound(p).distance <= xc.distance) {

#ifdef DBG_EXPL_TRANS
            std::cout << m_schedule.prettyLiteral(BOUND(p)) << " ("
                      << m_schedule.minDuration(
                             TASK(EVENT(m_schedule.getBound(p).l)))
                      << ") & "
                      << m_schedule.prettyLiteral(EDGE(disjunct[*xp][y]))
                      << std::endl;
#endif

            assert(m_schedule.getBound(p).distance <= xc.distance);
            Cl.push_back(BOUND(p));
            Cl.push_back(EDGE(disjunct[*xp][y]));
          }

#ifdef DBG_EXPL_TRANS
          else {
            std::cout << "none\n";
          }
#endif
        }
      }
    }
  }
}

template <typename T>
std::ostream &Transitivity<T>::display(std::ostream &os) const {
  os << "Transitivity";

#ifdef DEBUG_CONSTRAINT
  os << "[" << cons_id << "]";
#endif

  os << "(";
  for (auto t : m_tasks) {
    std::cout << " t" << t;
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
  //    m_schedule.displayLiteral(os, *l);
  //    ++l;
  //    while (l != explanations[h].end()) {
  //      os << ", ";
  //      m_schedule.displayLiteral(os, *l);
  //      ++l;
  //    }
  //  }
  //
  //  os << ")";
  return os;
}

// template <typename T> std::vector<int> Transitivity<T>::task_map;

} // namespace tempo

#endif
