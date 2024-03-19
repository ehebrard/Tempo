#ifndef TEMPO_TRANSITIVITY_HPP
#define TEMPO_TRANSITIVITY_HPP

#include <cassert>
#include <map>
#include <vector>

#include "Constraint.hpp"
#include "Explanation.hpp"
#include "Global.hpp"
#include "ReversibleObject.hpp"
#include "Scheduler.hpp"
#include "util/SparseSet.hpp"

//#define DBG_TRANSITIVITY

namespace tempo {

template <typename T> class Transitivity : public Constraint {
private:
  Scheduler<T> &m_schedule;
  std::vector<task> m_tasks;
  std::vector<int> task_map;

  //  std::vector<SparseSet<int, Reversible<size_t>>> forward;
  //  std::vector<SparseSet<int, Reversible<size_t>>> backward;

  std::vector<SparseSet<int, Reversible<size_t>>> DAG;

  //    std::vector<var>
  std::vector<std::vector<lit>> disjunct;

  //  std::vector<int> from;
  //  std::vector<int> to;

  SparseSet<> changed_pred;
  SparseSet<> changed_succ;
    
    std::vector<int> new_succ_of_x;
    std::vector<int> new_pred_of_y;
    
    bool change_flag{false};

  //  std::vector<lit> cur_explanation;
  //  std::vector<std::vector<lit>> explanations;

  //    int lvl;

public:
  template <typename ItTask, typename ItVar>
  Transitivity(Scheduler<T> &scheduler, const ItTask beg_task,
               const ItTask end_task, const ItVar beg_var, const ItVar end_var);
  virtual ~Transitivity();

  void add_edge(const int x, const int y, const lit r);

  bool notify_edge(const int lit, const int rank) override;
  void post(const int idx) override;
  void propagate() override;

  // void xplain(const lit l, const hint h, std::vector<lit> &Cl) const
  // override;
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
    : m_schedule(scheduler)
//    forward(Reversible<size_t>(0, m_schedule.getEnv()), Reversible<size_t>(0,
//    m_schedule.getEnv())),
// forward(Reversible<size_t>(0, m_schedule.getEnv()), Reversible<size_t>(0,
// m_schedule.getEnv())),
{

  task_map.resize(m_schedule.numTask());

  // get all tasks with non-zero duration
  auto i{0};
  for (auto j{beg_task}; j != end_task; ++j) {

    task t{*j};
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
  changed_pred.reserve(m_tasks.size());
  changed_succ.reserve(m_tasks.size());

  for (size_t i{0}; i < m_tasks.size(); ++i) {

    disjunct[i].resize(m_tasks.size());

    DAG.emplace_back(m_tasks.size(), &m_schedule.getEnv());
    DAG.back().fill();
    //        SparseSet<int, Reversible<size_t>> row(Reversible<size_t>(0,
    //        &m_schedule.getEnv()), Reversible<size_t>(0,
    //        &m_schedule.getEnv()), m_tasks.size()); forward.push_back(row);
    //    forward.emplace_back(m_tasks.size(), &m_schedule.getEnv());
    //    backward.emplace_back(m_tasks.size(), &m_schedule.getEnv());
    //        backward.emplace_back(Reversible<size_t>(0, &m_schedule.getEnv()),
    //        Reversible<size_t>(0, &m_schedule.getEnv()), m_tasks.size());
  }

  for (auto v{beg_var}; v != end_var; ++v) {
    auto ef{m_schedule.getEdge(POS(*v))};
    disjunct[task_map[TASK(ef.from)]][task_map[TASK(ef.to)]] = POS(*v);
    auto eb{m_schedule.getEdge(NEG(*v))};
    disjunct[task_map[TASK(eb.from)]][task_map[TASK(eb.to)]] = NEG(*v);
  }
}

template <typename T> Transitivity<T>::~Transitivity() {}

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
    for (size_t j{0}; j < m_tasks.size(); ++j)
      if (i != j) {
        //        lit l{m_schedule.getEdgeLit({START(m_tasks[i]),
        //        END(m_tasks[j])})};
        m_schedule.wake_me_on_edge(disjunct[i][j], cons_id);
        //        m_schedule.wake_me_on_edge(disjunct[j][i], cons_id);
      }
  }

  //    m_schedule.wake_me_on_edge(POS(e), cons_id);
  //    m_schedule.wake_me_on_edge(NEG(e), cons_id);
  //  }
}

template <typename T>
void Transitivity<T>::add_edge(const int x, const int y, const lit r) {

#ifdef DBG_TRANSITIVITY
  std::cout << " ==> add edge t" << m_tasks[x] << " -> t" << m_tasks[y]
            << std::endl;
#endif
    
change_flag=true;
    
  //            << m_schedule.prettyLiteral(EDGE(disjunct[x][y])) << std::endl;

  //    if(m_schedule.isUndefined(VAR(disjunct[x][y]))) {
  m_schedule.set(disjunct[x][y], {this, EDGE(r)});
  //    }

  //    else std::cout << "already here\n";

  assert(DAG[x].has(y));
  assert(DAG[y].has(x));

  DAG[x].remove_front(y);
  DAG[y].remove_back(x);

  if (not changed_pred.has(y))
    changed_pred.add(y);

  if (not changed_succ.has(x))
    changed_succ.add(x);
}

template <typename T>
bool Transitivity<T>::notify_edge(const lit l, const int) {
    change_flag=false;
  auto e{m_schedule.getEdge(l)};

#ifdef DBG_TRANSITIVITY
  std::cout << std::endl;
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
#endif

  auto x{task_map[TASK(e.from)]};
  auto y{task_map[TASK(e.to)]};

  if (DAG[x].isfront(y)) {
    assert(DAG[y].isback(x));
    return false;
  }

  // new edge x -> y
    // add egde between x and every successor of y
    // add egde between every predecessor of x and y

////  auto change_flag{false};
//    for (auto zp{DAG[y].frbegin()}; zp != DAG[y].frend(); ++zp) {
//        // there's a path x -> y -> z
//        auto z{*zp};
//        if (not DAG[x].isfront(z)) {
//            // there was no path x -> z
//            add_edge(x, z, l);
////            change_flag = true;
//            for (auto tp{DAG[x].brbegin()}; tp != DAG[x].brend(); ++tp) {
//                // there's a path t -> x -> y -> z
//                auto t{*tp};
//                if (not DAG[t].isfront(z)) {
//                    add_edge(t, z, l);
//                }
//            }
//        }
//    }
    
    for (auto zp{DAG[y].frbegin()}; zp != DAG[y].frend(); ++zp) {
        auto z{*zp};
        if(not DAG[x].isfront(z)) {
            add_edge(x, z, l);
        }
        for (auto tp{DAG[x].brbegin()}; tp != DAG[x].brend(); ++tp) {
            auto t{*tp};
            // there's a path t->x->y->z
            if(not DAG[t].isfront(y)) {
                add_edge(t, y, l);
            }
            if(not DAG[t].isfront(z)) {
                add_edge(t, z, l);
            }
        }
    }
 
  //  std::cout << forward[x] << std::endl;
  //
  //
  //    for(auto z : forward[y]) {
  //        if(not forward[x].has(v)) {
  //            forward[x].add(v);
  //            backward[v].add(x);
  //        }
  //    }

  //  // search for vertices reachable from y and not reachable from x
  //  changed.add(y);
  //  while (not changed.empty()) {
  //    auto u{changed.front()};
  //    changed.pop_front();
  //    for (auto v : forward[u])
  //      if (not changed.isfront(v) /*not visited*/ and
  //          not forward[x].has(v) /*not reachable from x*/) {
  //        changed.add(v);
  //      }
  //  }
  //    // changed contains all those vertives in "front"
  //    assert(not changed.isfront(x));
  //    changed.add(x);
  //    while (not changed.empty()) {
  //      auto u{changed.back()};
  //      changed.pop_back();
  //      for (auto v : backward[u])
  //        if (not changed.isfront(v) /*not visited*/) {
  //          changed.add(v);
  //        }
  //    }

  assert(not DAG[x].isfront(y));

  //  forward[x].add(y);
  //  backward[y].add(x);

  DAG[x].remove_front(y);
  DAG[y].remove_back(x);

#ifdef DBG_TRANSITIVITY
  if (change_flag)
    std::cout << "***\n";

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
#endif
if(change_flag)
{
    change_flag = false;
    return true;
}
    
  return false;
}

template <typename T> void Transitivity<T>::propagate() {
  auto cmax{m_schedule.upper(HORIZON)};

#ifdef DBG_TRANSITIVITY
    std::cout << "propagate transitivity:\n" << m_schedule << std::endl;
#endif
    
    try{
        for (auto x : changed_succ) {
            
#ifdef DBG_TRANSITIVITY
    std::cout << "t" << m_tasks[x] << " has new successors:\n";
#endif
            
            T est{cmax};
            T total_duration{0};
            for (auto yp{DAG[x].fbegin()}; yp != DAG[x].fend(); ++yp) {
                total_duration += m_schedule.minDuration(m_tasks[*yp]);
                est = std::min(est, m_schedule.lower(START(m_tasks[*yp])));
                
#ifdef DBG_TRANSITIVITY
    std::cout << " - t" << m_tasks[*yp] << " " 
                << m_schedule.minDuration(m_tasks[*yp])
                << "/" << m_schedule.lower(START(m_tasks[*yp]))
                << " -> " << total_duration << "/" << est << "\n";
#endif
            }
            
#ifdef DBG_TRANSITIVITY
    std::cout << "deduce " << prettyEvent(START(m_tasks[x])) << " >= "
            << (est + total_duration) << " (was "
            << m_schedule.lower(START(m_tasks[x])) << ")"<< std::endl;
#endif
            
            m_schedule.set({LOWERBOUND(START(m_tasks[x])), -est - total_duration}, {this, NoHint});
        }
        
        for (auto y : changed_pred) {
            
#ifdef DBG_TRANSITIVITY
    std::cout << "t" << m_tasks[y] << " has new predecessors:\n";
#endif
            
            T lct{0};
            T total_duration{0};
            for (auto xp{DAG[y].bbegin()}; xp != DAG[y].bend(); ++xp) {
                total_duration += m_schedule.minDuration(m_tasks[*xp]);
                lct = std::max(lct, m_schedule.upper(END(m_tasks[*xp])));
                
#ifdef DBG_TRANSITIVITY
    std::cout << " - t" << m_tasks[*xp] << " "
                << m_schedule.minDuration(m_tasks[*xp])
                << "/" << m_schedule.upper(END(m_tasks[*xp]))
                << " -> " << total_duration << "/" << lct << "\n";
#endif
                
            }
            
#ifdef DBG_TRANSITIVITY
    std::cout << "deduce " << prettyEvent(END(m_tasks[y])) << " <= "
            << (lct - total_duration) << " (was "
            << m_schedule.upper(END(m_tasks[y])) << ")" << std::endl;
#endif
            
            m_schedule.set({UPPERBOUND(END(m_tasks[y])), lct - total_duration},
                           {this, NoHint});
        }
    } catch(Failure &f) {
        changed_succ.clear();
        changed_pred.clear();
        throw f;
    }
    changed_succ.clear();
    changed_pred.clear();
}

template <typename T> int Transitivity<T>::getType() const {
  return TRANSITIVITYEXPL;
}

template <typename T>
void Transitivity<T>::xplain(const lit, const hint h, std::vector<lit> &Cl) {
    if(h != NoHint)
        Cl.push_back(h);
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
