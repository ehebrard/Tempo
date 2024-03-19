#ifndef TEMPO_DISJUNCTIVEPRECEDENCE_HPP
#define TEMPO_DISJUNCTIVEPRECEDENCE_HPP

#include <cassert>
#include <map>
#include <vector>

#include "Constraint.hpp"
#include "Explanation.hpp"
#include "Global.hpp"
#include "Scheduler.hpp"
#include "util/SparseSet.hpp"

#define DBG_GT

namespace tempo {

template <typename T>
class DisjunctivePrecedence : public Constraint, public Explainer {
private:
  Scheduler<T> &m_schedule;
  std::vector<task> m_tasks;
  //    std::vector<int> from;// for each var, the rank in m_tasks of the "from"
  //    event's task
  //        std::vector<int> to;// for each var, the rank in m_tasks of the "to"
  //        event's task
  //
  //        SparseSet<int, Reversible<size_t>> sources; // dynamic: the current
  //        sources in the precedence graph
  //    SparseSet<int, Reversible<size_t>> sinks; // dynamic: the current sinks
  //    in the precedence graph
  //

  SparseSet<> visited;
  std::vector<bool> concerned; // static: tasks requiring that resource
  //    SparseSet<int, Reversible<size_t>> sources; // dynamic: the current
  //    sources in the precedence graph SparseSet<int, Reversible<size_t>>
  //    sinks; // dynamic: the current sinks in the precedence graph
  //
  std::vector<var> edges;

  std::vector<lit> cur_explanation;
  std::vector<std::vector<lit>> explanations;

  //    int lvl;

public:
  template <typename ItTask, typename ItVar>
  DisjunctivePrecedence(Scheduler<T> &scheduler, const ItTask beg_task,
                        const ItTask end_task, const ItVar beg_var,
                        const ItVar end_var);
  virtual ~DisjunctivePrecedence();

  bool notify_edge(const int lit, const int rank) override;
  void post(const int idx) override;
  void propagate() override;

  template <class G> T DFS(const task s, const G &neighbors);

  template <class G>
  T getTrail(const task s, const G &neighbors
#ifdef DBG_GT
             ,
             int lvl
#endif
  );
  template <class G>
  T getHead(const task s, const G &neighbors
#ifdef DBG_GT
            ,
            int lvl
#endif
  );

  // void xplain(const lit l, const hint h, std::vector<lit> &Cl) const
  // override;
  void xplain(const lit l, const hint h, std::vector<lit> &Cl) override;
  int getType() const override;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

  static std::vector<T> task_intrail;
  static std::vector<T> task_outtrail;
  static std::vector<T> task_inhead;
  static std::vector<T> task_outhead;
  static std::vector<int> task_map;
  static SparseSet<> dfs_stack;
#ifdef DEBUG_CONSTRAINT
  int debug_flag{2};
#endif
};

template <typename T>
template <typename ItTask, typename ItVar>
DisjunctivePrecedence<T>::DisjunctivePrecedence(Scheduler<T> &scheduler,
                                                const ItTask beg_task,
                                                const ItTask end_task,
                                                const ItVar beg_var,
                                                const ItVar end_var)
    : m_schedule(scheduler) {

  dfs_stack.reserve(
      m_schedule.numEvent()); // static_cast<size_t>(end_task - beg_task));

  task_map.resize(m_schedule.numTask());
  concerned.resize(m_schedule.numTask(), false);
  visited.reserve(m_schedule.numTask());
  task_intrail.resize(m_schedule.numTask(), 0);
  task_inhead.resize(m_schedule.numTask(), 0);
  task_outtrail.resize(m_schedule.numTask(), 0);
  task_outhead.resize(m_schedule.numTask(), 0);

  // get all tasks with non-zero duration
  auto i{0};
  for (auto j{beg_task}; j != end_task; ++j) {

    task t{*j};
    task_map[t] = i++;
    m_tasks.push_back(t);

    if (m_schedule.maxDuration(t) > 0) {
      concerned[t] = true;
    }
  }

  for (auto v{beg_var}; v != end_var; ++v) {
    edges.push_back(*v);

    auto e{m_schedule.getEdge(POS(*v))};

    from.push_back(task_map[e.from]);
    to.push_back(task_map[e.to]);
  }
  //
  //#ifdef DEBUG_CONSTRAINT
  //  debug_flag = 2;
  //#endif
}

template <typename T> DisjunctivePrecedence<T>::~DisjunctivePrecedence() {}

template <typename T> void DisjunctivePrecedence<T>::post(const int idx) {

#ifdef DEBUG_CONSTRAINT
  if (debug_flag > 0) {
    std::cout << "post " << *this << std::endl;
  }
#endif

  cons_id = idx;

  for (auto e : edges) {
    m_schedule.wake_me_on_edge(POS(e), cons_id);
    m_schedule.wake_me_on_edge(NEG(e), cons_id);

    //      std::cout << "post on " << m_schedule.prettyLiteral(EDGE(POS(e))) <<
    //      " and " << m_schedule.prettyLiteral(EDGE(NEG(e))) << std::endl;
    //
    // might be useful to be called on bounds as well (but not on jobshop/open
    // shop)
  }
}

template <typename T>
bool DisjunctivePrecedence<T>::notify_edge(const lit, const int) {
  return true;
}

/*
 TO EXPLORE : IN
 STACK : FRONT
 EXPLORED : BACK

 s|....|
 u,v,w|...|s
 u,v|...|w,s

 */

template <typename T>
template <class G>
T DisjunctivePrecedence<T>::DFS(const G &neighbors) {
  while (dfs_stack.frontsize() > 0) {
    auto x = *(dfs_stack.frbegin());
    dfs_stack.front_to_back(x);
    for (auto e : neighbors[x]) {
      if (e.label() <= 0 and dfs_stack.has(e)) {
        dfs_stack.remove_front(e);
      }
    }
  }
}

template <typename T>
template <class G>
bool DisjunctivePrecedence<T>::isSink(const event s, const G &neighbors) {
  for (auto e : neighbors[s]) {
    if (e.label() <= 0 and dfs_stack.has(e)) {
      return false;
    }
  }
  return true;
}

template <typename T>
template <class G>
T DisjunctivePrecedence<T>::topologicalForwardSort() {
  auto &F{m_schedule.domain.getForwardGraph()};
  auto &B{m_schedule.domain.getBackwardGraph()};
  dfs_stack.clear();
  for (auto t : m_tasks) {
    dfs_stack.add(START(t));
    dfs_stack.add(END(t));
  }
  for (auto t : m_tasks) {
    if (isSink(END(t), B))
      dfs_stack.remove_front(END(t));
  }
  DFS(N);
  std::cout << dfs_stack << std::endl;
}

template <typename T> T DisjunctivePrecedence<T>::ForwardSort() {
  auto &N{m_schedule.domain.getForwardGraph()};
  dfs_stack.clear();
  for (auto t : m_tasks) {
    dfs_stack.add(t);
    auto source{true};
    for (auto e : N[END(t)]) {
      if (e.label() <= 0 and concerned[e]) {
        source = false;
        break;
      }
    }
    if (source)
      dfs_stack.remove(t);
  }
  DFS(N);
  std::cout << dfs_stack << std::endl;
}

template <typename T>
template <class G>
T DisjunctivePrecedence<T>::getTrail(const task s, const G &neighbors
#ifdef DBG_GT
                                     ,
                                     int lvl
#endif
) {

#ifdef DEBUG_CONSTRAINT
  //    if(++debug_flag > 100)
  //        exit(1);
  if (debug_flag > 1) {
#ifdef DBG_GT
    for (auto i{0}; i < lvl; ++i)
      std::cout << " ";
#endif

    std::cout << "trail of t" << s << "=\n";
  }
#endif

  if (visited.has(s)) {
#ifdef DEBUG_CONSTRAINT
    if (debug_flag > 0) {

#ifdef DBG_GT
      for (auto i{0}; i < lvl; ++i)
        std::cout << " ";
#endif
      std::cout << "already counted (0)\n";
    }
#endif
    return 0;
  }

  visited.add(s);

  if (task_intrail[s] != 0) {
#ifdef DEBUG_CONSTRAINT
    if (debug_flag > 1) {
#ifdef DBG_GT
      for (auto i{0}; i < lvl; ++i)
        std::cout << " ";
#endif
      std::cout << "already computed (" << task_intrail[s] << ")\n";
    }
#endif
    return task_intrail[s];
  }

  T trail{0};
  T outtrail{-1}; // m_schedule.upper(HORIZON)};
  for (auto e : neighbors[END(s)]) {
    if (e.label() <= 0) {

      task t{TASK(e)};

      if (t >= 0 and concerned[t]) {

        trail += m_schedule.minDuration(t);

#ifdef DEBUG_CONSTRAINT
        if (debug_flag > 1) {
#ifdef DBG_GT
          for (auto i{0}; i < lvl; ++i)
            std::cout << " ";
#endif
          std::cout << "explore predecessor t" << t << " (" << trail << "+"
                    << m_schedule.minDuration(t) << ")\n";
        }
#endif

        trail += getTrail(t, neighbors
#ifdef DBG_GT
                          ,
                          lvl + 1
#endif
        );
        outtrail = std::max(outtrail, task_outtrail[t]);
      }
    }
  }

  if (outtrail == -1) {
#ifdef DEBUG_CONSTRAINT
    if (debug_flag > 1) {
#ifdef DBG_GT
      for (auto i{0}; i < lvl; ++i)
        std::cout << " ";
#endif
      std::cout << "sink! (0)\n";
    }
#endif
    task_outtrail[s] = m_schedule.upper(END(s));
    return 0;
  }

#ifdef DEBUG_CONSTRAINT
  if (debug_flag > 1) {
#ifdef DBG_GT
    for (auto i{0}; i < lvl; ++i)
      std::cout << " ";
#endif
    std::cout << "trail(" << s << ") = " << trail << "\n";
  }
#endif

  task_outtrail[s] = outtrail;
  return task_intrail[s] = trail;
}

template <typename T>
template <class G>
T DisjunctivePrecedence<T>::getHead(const task s, const G &neighbors
#ifdef DBG_GT
                                    ,
                                    int lvl
#endif
) {

#ifdef DEBUG_CONSTRAINT
  //    if(++debug_flag > 100)
  //        exit(1);
  if (debug_flag > 1) {
#ifdef DBG_GT
    for (auto i{0}; i < lvl; ++i)
      std::cout << " ";
#endif
    std::cout << "head of t" << s << "=\n";
  }
#endif

  if (visited.has(s)) {
#ifdef DEBUG_CONSTRAINT
    if (debug_flag > 0) {
#ifdef DBG_GT
      for (auto i{0}; i < lvl; ++i)
        std::cout << " ";
#endif
      std::cout << "already counted\n";
    }
#endif
    return 0;
  }

  visited.add(s);

  //  if (neighbors[END(s)].empty()) {
  //#ifdef DEBUG_CONSTRAINT
  //    if (debug_flag > 0) {
  //      std::cout << " sink! (0)\n";
  //    }
  //#endif
  //    task_outtrail[s] = m_schedule.upper(END(s));
  //    return 0;
  //  }
  if (task_inhead[s] != 0) {
#ifdef DEBUG_CONSTRAINT
    if (debug_flag > 1) {
#ifdef DBG_GT
      for (auto i{0}; i < lvl; ++i)
        std::cout << " ";
#endif
      std::cout << "already computed (" << task_inhead[s] << ")\n";
    }
#endif
    return task_inhead[s];
  }
  T head{0};
  T outhead{m_schedule.upper(HORIZON)};
  for (auto e : neighbors[START(s)]) {
    if (e.label() <= 0) {

      task t{TASK(e)};

      //#ifdef DEBUG_CONSTRAINT
      //    if (debug_flag > 1) {
      //      std::cout << " explore edge " << prettyEvent(START(s)) << " -> "
      //      << prettyEvent(int(e)) << " (t" << t << ")\n";
      //    }
      //#endif

      //            task t{TASK(e)};
      if (t >= 0 and concerned[t]) {

        head += m_schedule.minDuration(t);

#ifdef DEBUG_CONSTRAINT
        if (debug_flag > 1) {
#ifdef DBG_GT
          for (auto i{0}; i < lvl; ++i)
            std::cout << " ";
#endif
          std::cout << "explore successor t" << t << " (" << head << "+"
                    << m_schedule.minDuration(t) << ")\n";
        }
#endif

        head += getHead(t, neighbors
#ifdef DBG_GT
                        ,
                        lvl + 1
#endif
        );
        outhead = std::min(outhead, task_outhead[t]);
      }
    }
  }

  if (outhead == m_schedule.upper(HORIZON)) {
#ifdef DEBUG_CONSTRAINT
    if (debug_flag > 1) {
#ifdef DBG_GT
      for (auto i{0}; i < lvl; ++i)
        std::cout << " ";
#endif
      std::cout << "sink! (0)\n";
    }
#endif
    task_outhead[s] = m_schedule.lower(START(s));
    return 0;
  }

#ifdef DEBUG_CONSTRAINT
  if (debug_flag > 1) {
#ifdef DBG_GT
    for (auto i{0}; i < lvl; ++i)
      std::cout << " ";
#endif
    std::cout << "head(" << s << ") = " << head << "\n";
  }
#endif

  task_outhead[s] = outhead;
  return task_inhead[s] = head;
}

template <typename T> void DisjunctivePrecedence<T>::propagate() {

#ifdef DEBUG_CONSTRAINT
  if (debug_flag > 0) {
    std::cout << std::endl << "propagate precedences()\n";
  }
#endif

  for (auto t : m_tasks) {

    //      std::cout << t << "/" << task_intrail.size() << std::endl;

    task_intrail[t] = 0;
    task_inhead[t] = 0;
  }

  //    for(auto t : sources) {
  //        task_trail[t] = getTrail(END(t), m_schedule.getForwardGraph());
  //    }
  //
  //  for (auto t : m_tasks) {
  //    task_intrail[t] = 0;
  //    //        task_head[t] = 0;
  //  }

  for (auto t : m_tasks) {
    //      debug_flag = 1;

    visited.clear();
    //      lvl = 0;

#ifdef DEBUG_CONSTRAINT
    if (debug_flag > 1) {
      std::cout << std::endl;
    }
#endif
    task_intrail[t] = getTrail(t, m_schedule.domain.getBackwardGraph()
#ifdef DBG_GT
                                      ,
                               0
#endif
    );
    visited.clear();
    //      lvl = 0;

#ifdef DEBUG_CONSTRAINT
    if (debug_flag > 1) {
      std::cout << std::endl;
    }
#endif
    task_inhead[t] = getHead(t, m_schedule.domain.getForwardGraph()
#ifdef DBG_GT
                                    ,
                             0
#endif
    );

    auto new_ub{(m_schedule.upper(HORIZON) - task_intrail[t])};
    auto cur_ub{m_schedule.upper(END(t))};
    auto better_ub{(task_outtrail[t] - task_intrail[t])};

    std::cout << "task t" << t << " " << cur_ub << "/" << new_ub << "/"
              << better_ub << (cur_ub > better_ub ? "**\n" : "\n");

    auto new_lb{task_inhead[t]};
    auto cur_lb{m_schedule.lower(START(t))};
    auto better_lb{task_outhead[t] + task_inhead[t]};

    std::cout << "task t" << t << " " << cur_lb << "/" << new_lb << "/"
              << better_lb << (cur_lb < better_lb ? "**\n" : "\n");

    assert(better_ub <= new_ub);
    assert(better_lb >= new_lb);

    if (better_ub < 0) {
      std::cout << m_schedule << std::endl;

      exit(1);
    }
    assert(better_lb >= 0);

    //      task_inhead[t] = getHead(t, m_schedule.domain.getForwardGraph());
  }

  //
  //#ifdef DEBUG_CONSTRAINT
  //  if (debug_flag > 0) {
  //    std::cout << std::endl << "ST\n";
  //    for (auto i : est_order) {
  //      std::cout << prettyEvent(START(m_tasks[i])) << ": [" << est(i)
  //                << ".." << lst(i) << "]\n";
  //    }
  //
  //    std::cout << "\nCT\n";
  //    for (auto i : lct_order) {
  //      std::cout << prettyEvent(END(m_tasks[i])) << ": [" << ect(i) << ".."
  //                << lct(i) << "]\n";
  //    }
  //  }
  //#endif
}

template <typename T> int DisjunctivePrecedence<T>::getType() const {
  return EDGEFINDINGEXPL;
}

template <typename T>
void DisjunctivePrecedence<T>::xplain(const lit, const hint h,
                                      std::vector<lit> &Cl) {
  for (auto r : explanations[h]) {
    Cl.push_back(r);
  }
}

template <typename T>
std::ostream &DisjunctivePrecedence<T>::display(std::ostream &os) const {
  os << "Disjunctive Precedence";

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
std::ostream &DisjunctivePrecedence<T>::print_reason(std::ostream &os,
                                                     const hint) const {
  //  display(os);
  os << "disj-prec";
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

template <typename T> std::vector<int> DisjunctivePrecedence<T>::task_map;
template <typename T> std::vector<T> DisjunctivePrecedence<T>::task_intrail;
template <typename T> std::vector<T> DisjunctivePrecedence<T>::task_outtrail;
template <typename T> std::vector<T> DisjunctivePrecedence<T>::task_inhead;
template <typename T> std::vector<T> DisjunctivePrecedence<T>::task_outhead;

} // namespace tempo

#endif
