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
//  std::vector<task> m_tasks;
    std::vector<Task<T>*> the_tasks;
  //    std::vector<int> task_from; // for each disjunct (in the order they are
  //    declared as triggers, the 'from' event) std::vector<int> task_to; // for
  //    each disjunct (in the order they are declared as triggers, the 'to'
  //    event)

  // encoding: (y \in front(DAG[x]) && x \in back(y)) <=> edge (x,y)
  std::vector<SparseSet<int, Reversible<size_t>>> DAG;
  // encoding: y \in transitive_reduction[x] <=> edge (x,y)
  SparseSet<int, Reversible<size_t>> transitive_reduction;
  DisjointSet<> forest;
  std::vector<std::vector<T>> distance_matrix;

  std::vector<std::vector<lit>> disjunct;
  std::vector<int> scopex;
  std::vector<int> scopey;

  std::vector<int> sorted_tasks;
  std::vector<T> offset;

  std::vector<int> new_succ_of_x;
  std::vector<int> new_pred_of_y;

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
  Transitivity(Scheduler<T> &scheduler,
               const ItTask beg_task,
               const ItTask end_task,
               const ItVar beg_var);
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
Transitivity<T>::Transitivity(Scheduler<T> &scheduler, 
                              const ItTask beg_task,
                              const ItTask end_task,
                              const ItVar beg_var)
    : m_schedule(scheduler),
      transitive_reduction(std::distance(beg_task, end_task) *
                               std::distance(beg_task, end_task),
                           &m_schedule.getEnv()) {

  priority = Priority::Low;

//  auto n{std::distance(beg_task, end_task)};
//  auto m{std::distance(beg_var, end_var)};
//  assert(m = n * (n - 1) / 2);

  // get all tasks with non-zero duration
  //  size_t i{0}, j{0};
//  for (auto jp{beg_task}; jp != end_task; ++jp) {
//
//    task t{*jp};
//    m_tasks.push_back(t);
//  }
          
          for (auto jp{beg_task}; jp != end_task; ++jp) {
            the_tasks.push_back(*jp);
          }

  disjunct.resize(the_tasks.size());
  sorted_tasks.resize(the_tasks.size());
  offset.resize(the_tasks.size());

  transitive_reduction.fill(); // resize(m_tasks.size() * m_tasks.size());
  forest.resize(the_tasks.size());

  for (size_t i{0}; i < the_tasks.size(); ++i) {

    disjunct[i].resize(the_tasks.size());

    DAG.emplace_back(the_tasks.size(), &m_schedule.getEnv());
    DAG.back().fill();

    transitive_reduction.remove_back(edge(i, i));
  }

  auto ep{beg_var};
  for (auto ip{beg_task}; ip != end_task; ++ip) {
    for (auto jp{ip + 1}; jp != end_task; ++jp) {
      auto x{*ep};

      auto i{std::distance(beg_task, ip)};
      auto j{std::distance(beg_task, jp)};
      disjunct[i][j] = NEG(x);
      disjunct[j][i] = POS(x);

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

  for (size_t i{0}; i < the_tasks.size(); ++i) {
//    m_schedule.wake_me_on_event(LOWERBOUND(START(m_tasks[i])), cons_id);
//    m_schedule.wake_me_on_event(UPPERBOUND(END(m_tasks[i])), cons_id);
      m_schedule.wake_me_on_event(LOWERBOUND(the_tasks[i]->getStart()), cons_id);
      m_schedule.wake_me_on_event(UPPERBOUND(the_tasks[i]->getEnd()), cons_id);

      
    for (size_t j{0}; j < the_tasks.size(); ++j)
      if (i != j) {

        m_schedule.wake_me_on_edge(disjunct[i][j], cons_id);
        scopex.push_back(i);
        scopey.push_back(j);
      }
  }
}

template <typename T>
void Transitivity<T>::add_edge(const int x, const int y, const int r) {

#ifdef DBG_TRANSITIVITY
  if (DBG_TRANSITIVITY)
    std::cout << " ==> add edge t" << the_tasks[x]->id() << " -> t" << the_tasks[y]->id()
              << std::endl;
#endif

  if (DAG[x].frontsize() > 0 or DAG[y].backsize() > 0)
    change_flag = true;

#ifdef DBG_LTRANS
  if (m_schedule.isUndefined(VAR(disjunct[x][y]))) {
    std::cout << "edge pruning\n";
  }
#endif

  m_schedule.set(disjunct[x][y], {this, r});

  assert(DAG[x].has(y));
  assert(DAG[y].has(x));

  DAG[x].remove_front(y);
  DAG[y].remove_back(x);

  if (transition_flag) {
    transitive_reduction.remove_back(edge(y, x));
  }
}

template <typename T> bool Transitivity<T>::notify_bound(const lit, const int) {
  return true;
}

template <typename T>
bool Transitivity<T>::notify_edge(const lit, const int r) {

  auto x{scopex[r]};
  auto y{scopey[r]};

#ifdef DBG_TRANSITIVITY
  if (DBG_TRANSITIVITY) {
    std::cout << std::endl;
    for (size_t i{0}; i < the_tasks.size(); ++i) {
      std::cout << the_tasks[i]->id() << ":";
      for (auto j{DAG[i].fbegin()}; j != DAG[i].fend(); ++j) {
        std::cout << " -> t" << the_tasks[*j]->id();
      }
      std::cout << std::endl;
    }

    for (size_t i{0}; i < the_tasks.size(); ++i) {
      std::cout << the_tasks[i]->id() << ":";
      for (auto j{DAG[i].bbegin()}; j != DAG[i].bend(); ++j) {
        std::cout << " <- " << the_tasks[*j]->id();
      }
      std::cout << std::endl;
    }

    std::cout << "notify edge " << *the_tasks[x] << " -> " << *the_tasks[y]
      << " / " << m_schedule.prettyLiteral(EDGE(l))
              << std::endl;
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

    for (size_t i{0}; i < the_tasks.size(); ++i) {
      std::cout << "t" << the_tasks[i]->id() << ":";
      for (auto j{DAG[i].fbegin()}; j != DAG[i].fend(); ++j) {
        std::cout << " -> t" << the_tasks[*j]->id();
      }
      std::cout << std::endl;
    }

    for (size_t i{0}; i < the_tasks.size(); ++i) {
      std::cout << "t" << the_tasks[i]->id() << ":";
      for (auto j{DAG[i].bbegin()}; j != DAG[i].bend(); ++j) {
        std::cout << " <- t" << the_tasks[*j]->id();
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
      std::cout << "t" << the_tasks[i]->id() << ":";
      for (auto j{DAG[i].fbegin()}; j != DAG[i].fend(); ++j) {
        std::cout << " -> t" << the_tasks[*j]->id();
      }
      std::cout << std::endl;
    }

    std::cout << "TRED:\n";
    for (size_t i{0}; i < the_tasks.size(); ++i) {
      std::cout << "t" << the_tasks[i]->id() << ":";
      for (size_t j{0}; j < the_tasks.size(); ++j) {
        if (transitive_reduction.has(edge(i, j)))
          std::cout << " -> t" << the_tasks[j]->id();
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

  T min_start{INFTY};
//  for (auto t : m_tasks) {
//    path_length += m_schedule.minDuration(t);
//    min_start = std::min(min_start, m_schedule.lower(START(t)));
//  }
    for (auto t : the_tasks) {
      path_length += t->minDuration();
      min_start = std::min(min_start, t->getEarliestStart());
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

  for (size_t x{0}; x < the_tasks.size(); ++x) {
    sorted_tasks[x] = x;
    offset[x] = 0;
  }
  std::sort(sorted_tasks.begin(), sorted_tasks.end(),
            [&](const int x, const int y) -> bool {
      //              return m_schedule.upper(END(m_tasks[x])) <
      //                     m_schedule.upper(END(m_tasks[y]));
      return the_tasks[x]->getLatestEnd() <
      the_tasks[y]->getLatestEnd();});

#ifdef DBG_TRANSITIVITY
    if (DBG_TRANSITIVITY) {
        std::cout << "\npropagate w.r.t. subsquent tasks\n";
    }
#endif
    

  for (auto x : sorted_tasks) {

#ifdef DBG_TRANSITIVITY
    if (DBG_TRANSITIVITY) {
      std::cout << "t" << the_tasks[x]->id() << " ["
        << the_tasks[x]->getEarliestStart() << ".."
        << the_tasks[x]->getLatestEnd() << "] ("
        << the_tasks[x]->minDuration() << ")";
//                << m_schedule.lower(START(m_tasks[x])) << ".."
//                << m_schedule.upper(END(m_tasks[x])) << "] ("
//                << m_schedule.minDuration(m_tasks[x]) << ")";
      for (auto yp{DAG[x].fbegin()}; yp != DAG[x].fend(); ++yp) {
        std::cout << " -> t" << the_tasks[*yp]->id();
      }
      std::cout << std::endl;
    }
#endif

    // get all the successors y of x in the graph, (so predecessors in the
    // schedule: y < x)
    for (auto yp{DAG[x].fbegin()}; yp != DAG[x].fend(); ++yp) {
//      offset[*yp] += m_schedule.minDuration(m_tasks[x]);
        offset[*yp] += the_tasks[x]->minDuration();

      // x must end before ex
//      auto ex{m_schedule.upper(END(m_tasks[x]))};
        auto ex{the_tasks[x]->getLatestEnd()};

      // y must end before ey
//      auto ey{m_schedule.upper(END(m_tasks[*yp]))};
        auto ey{the_tasks[*yp]->getLatestEnd()};

      // y must end before ex -
      if ((ex - offset[*yp]) < ey) {
#ifdef DBG_LTRANS
        pruning = true;
#endif
          
          
#ifdef DBG_TRANSITIVITY
        if (DBG_TRANSITIVITY) {
            std::cout << " new bound (" << *the_tasks[*yp] << " must end before " <<
            ex << " - " << offset[*yp] << ") " ;
            m_schedule.displayBound(std::cout, the_tasks[*yp]->end.before(ex - offset[*yp]));
            //                  << " <= " << ex - offset[*yp]
            std::cout << "/" << ey << std::endl;
        }
#endif
          
//        m_schedule.set({UPPERBOUND(END(m_tasks[*yp])), ex - offset[*yp]},
//                       {this, offset[*yp]});
          
//          std::cout << "hello\n";
          
          auto bc{the_tasks[*yp]->end.before(ex - offset[*yp])};
          
//          m_schedule.set(bc, {this, offset[*yp]});
          m_schedule.set(bc, {this, ex});

      }
    }
//        }
//#ifdef DBG_TRANSITIVITY
//    if (DBG_TRANSITIVITY)
//      std::cout << std::endl;
//#endif
  }

#ifdef DBG_TRANSITIVITY
    if (DBG_TRANSITIVITY) {
        std::cout << "\npropagate w.r.t. preceding tasks\n";
    }
#endif

  for (size_t x{0}; x < the_tasks.size(); ++x) {
    //        sorted_tasks[x] = x;
    offset[x] = 0;
  }
  std::sort(sorted_tasks.begin(), sorted_tasks.end(),
            [&](const int x, const int y) -> bool {
//              return m_schedule.lower(START(m_tasks[x])) >
//                     m_schedule.lower(START(m_tasks[y]));
      return the_tasks[x]->getEarliestStart() >
      the_tasks[y]->getEarliestStart();
            });

  //    for(auto x : sorted_tasks) {
  //        std::cout << "[" << m_schedule.lower(START(m_tasks[x])) << ".." <<
  //        m_schedule.upper(END(m_tasks[x])) << "]\n";
  //    }

  for (auto y : sorted_tasks) {
#ifdef DBG_TRANSITIVITY
    if (DBG_TRANSITIVITY) {
      std::cout << "t" << the_tasks[y]->id() << " ["
                << the_tasks[y]->getEarliestStart() << ".."
                << the_tasks[y]->getLatestEnd() << "] ("
                << the_tasks[y]->minDuration() << ")";
      //        if(DAG[x].frontsize() > 0) {
      for (auto xp{DAG[y].bbegin()}; xp != DAG[y].bend(); ++xp) {
        std::cout << " <- t" << the_tasks[*xp]->id();
      }
    }
    std::cout << std::endl;
#endif
    //        if(DAG[x].frontsize() > 0) {
    for (auto xp{DAG[y].bbegin()}; xp != DAG[y].bend(); ++xp) {
        offset[*xp] += the_tasks[y]->minDuration(); //m_schedule.minDuration(m_tasks[y]);

//      auto sy{m_schedule.lower(START(m_tasks[y]))};
//      auto sz{m_schedule.lower(START(m_tasks[*xp]))};
        auto sy{the_tasks[y]->getEarliestStart()};
        auto sz{the_tasks[*xp]->getEarliestStart()};


      if ((sy + offset[*xp]) > sz) {
#ifdef DBG_LTRANS
        pruning = true;
#endif
          
          
//          (" << *the_tasks[*yp] << " must end before " <<
//          ex << " - " << offset[*yp] << ") " ;
//          m_schedule.displayBound(std::cout, the_tasks[*yp]->end.before(ex - offset[*yp]));
          
#ifdef DBG_TRANSITIVITY
      if (DBG_TRANSITIVITY)
           std::cout << " new bound (" << *the_tasks[*xp] << " must start after " << sy << " + " << offset[*xp] << ") " ;
           m_schedule.displayBound(std::cout, the_tasks[*xp]->start.after(sy + offset[*xp]));
//          << prettyEvent(START(m_tasks[*xp]))
//                  << " >= " << sy + offset[*xp]
          std::cout << "/" << sz << std::endl;
#endif
          
//        m_schedule.set({LOWERBOUND(START(m_tasks[*xp])), -sy - offset[*xp]},
//                       {this, offset[*xp]});
          m_schedule.set(the_tasks[*xp]->start.after(sy + offset[*xp]),
//                         {this, offset[*xp]});
                         {this, sy});
      }
    }
//        }
//#ifdef DBG_TRANSITIVITY
//    if (DBG_TRANSITIVITY)
//      std::cout << std::endl;
//#endif
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

//      lit p{m_schedule.getEdgeLit({el.from, NOT(er.from)})};
        lit p{m_schedule.getEdgeLit({el.from, er.from})};

#ifdef DBG_EXPL_TRANS
      std::cout << " and " << m_schedule.prettyLiteral(EDGE(p));
      std::cout.flush();
#endif

      Cl.push_back(EDGE(p));
    }
    if (el.to != er.to) {

      //lit p{m_schedule.getEdgeLit({NOT(er.to), el.to})};
        lit p{m_schedule.getEdgeLit({er.to, el.to})};

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



    if (SIGN(bc.l) == LOWER) {
        
#ifdef DBG_EXPL_TRANS
      std::cout << " reason was all predecessors before " << h << std::endl;
#endif
        
//        std::cout << "hello\n";
        
//      auto t{TASK(EVENT(bc.l))};
        Task<T>& t{m_schedule.getTask(EVENT(bc.l))};
      task x{-1};
      for (unsigned i{0}; i < the_tasks.size(); ++i) {
        if (the_tasks[i]->id() == t.id()) {
          x = i;
          break;
        }
      }

      for (auto yp{DAG[x].fbegin()}; yp != DAG[x].fend(); ++yp) {
//        BoundConstraint<T> yc_old{LIT(START(m_tasks[*yp]), LOWER), bc.distance + h};
//          BoundConstraint<T> yc{the_tasks[*yp]->start.after(-bc.distance - h)};
          BoundConstraint<T> yc{the_tasks[*yp]->start.after(h)};

//          assert(yc_old == yc);
          
#ifdef DBG_EXPL_TRANS
        std::cout << " implicant of " << yc << ": ";
#endif

        auto p{m_schedule.getImplicant(yc)};

        if (p != NoLit) {
          if (p < FROM_GEN(l) and
              m_schedule.getBound(p).distance <= yc.distance) {

#ifdef DBG_EXPL_TRANS
            std::cout << m_schedule.prettyLiteral(BOUND(p)) << " ("
//                      << m_schedule.minDuration(
//                             TASK(EVENT(m_schedule.getBound(p).l)))
              << m_schedule.getTask(EVENT(m_schedule.getBound(p).l)).minDuration()
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
        
#ifdef DBG_EXPL_TRANS
      std::cout << " reason was all predecessors after " << h << std::endl;
#endif
        
        Task<T>& t{m_schedule.getTask(EVENT(bc.l))};
//      auto t{TASK(EVENT(bc.l))};
      task y{-1};
      for (unsigned i{0}; i < the_tasks.size(); ++i) {
        if (the_tasks[i]->id() == t.id()) {
          y = i;
          break;
        }
      }
      for (auto xp{DAG[y].bbegin()}; xp != DAG[y].bend(); ++xp) {
//        BoundConstraint<T> xc_old{LIT(END(m_tasks[*xp]), UPPER), bc.distance + h};
//          BoundConstraint<T> xc{the_tasks[*xp]->end.before(bc.distance + h)};
          BoundConstraint<T> xc{the_tasks[*xp]->end.before(h)};
          
//          assert(xc_old == xc);

#ifdef DBG_EXPL_TRANS
        std::cout << " implicant of " << xc << ": ";
#endif

        auto p{m_schedule.getImplicant(xc)};

        if (p != NoLit) {
          if (p < FROM_GEN(l) and
              m_schedule.getBound(p).distance <= xc.distance) {

#ifdef DBG_EXPL_TRANS
            std::cout << m_schedule.prettyLiteral(BOUND(p)) << " ("
//                      << m_schedule.minDuration(
//                             TASK(EVENT(m_schedule.getBound(p).l)))
              << m_schedule.getTask(EVENT(m_schedule.getBound(p).l)).minDuration()
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
  for (auto t : the_tasks) {
    std::cout << " t" << t->id();
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
