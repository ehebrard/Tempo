
#ifndef _TEMPO_SCHEDULER_HPP
#define _TEMPO_SCHEDULER_HPP


#include <sstream>
#include <fstream>

#include "ClauseBase.hpp"
#include "Constant.hpp"
#include "ConstraintQueue.hpp"
#include "DistanceConstraint.hpp"
#include "Global.hpp"
#include "Objective.hpp"
#include "Restart.hpp"
#include "Task.hpp"
#include "TemporalNetwork.hpp"
#include "constraints/DisjunctiveEdgeFinding.hpp"
#include "constraints/EdgeConstraint.hpp"
#include "constraints/Transitivity.hpp"
#include "heuristics/HeuristicManager.hpp"
#include "heuristics/ValueHeuristicsManager.hpp"
#include "util/Heap.hpp"
#include "util/KillHandler.hpp"
#include "util/Options.hpp"
#include "util/SubscribableEvent.hpp"
#include "util/parsing/format.hpp"

//#define DBG_MINIMIZATION

namespace tempo {

template<typename T>
class Scheduler : public Explainer, public ReversibleObject
{

public:
  /**
   * @name constructors
   */
  ///@{
  Scheduler(Options opt);
  ~Scheduler();
  ///@}

  /**
   * @name accessors
   */
  ///@{
  /// Number of tasks created
  size_t numTask() const;
  /// Number of event
  size_t numEvent() const;
  /// Number of Boolean variables
  size_t numVariable() const;
  /// Number of learnt clauses
  size_t numClause() const;
  /// Number of constraints
  size_t numConstraint() const;
  /// Number of bound changes
  size_t numBoundLiteral() const;
  /// Number of Boolean changes
  size_t numEdgeLiteral() const;
  ///@}

  /**
   * @name modelling methods
   */
  ///@{
  void load(ProblemInstance &p);
//  event newEvent();
  task newTask(const T min_dur, const T max_dur);
  task newTask(const T release, const T deadline, const T min_dur,
               const T max_dur);
  //    Task<T> getTask(const task i) {return tasks[i];}
  var newVariable(const DistanceConstraint<T> &if_true = Constant::NoEdge<T>,
                  const DistanceConstraint<T> &if_false = Constant::NoEdge<T>);

  void newPrecedence(event x, event y, T delay);
  void newMaximumLag(event x, event y, T maxlag);
    void addConstraint(const DistanceConstraint<T>& c);

  template <typename Container> void addClause(Container &c);

  template <typename ItTask, typename ItVar>
  void postEdgeFinding(const ItTask beg_task, const ItTask end_task,
                       const ItVar beg_var
    );

  template <typename ItTask, typename ItVar>
  void postTransitivityReasoning(const ItTask beg_task, const ItTask end_task,
                                 const ItVar beg_var);
    
    

  bool value(const var x) const;
  bool isTrue(const var x) const;
  bool isFalse(const var x) const;
  bool isUndefined(const var x) const;
  bool falsified(const lit l) const;
  bool satisfied(const lit l) const;

  bool falsified(const BoundConstraint<T> &c) const;
  bool satisfied(const BoundConstraint<T> &c) const;

  T upper(const event) const;
  T lower(const event) const;
  T bound(const bool s, const event) const;
//  T minDuration(const task) const;
//  T maxDuration(const task) const;
  T distance(const event x, const event y) const;
  //  T getMakespan() const;

  auto getSolution() const noexcept -> const std::vector<bool> &;

  // get edge literal from an edge propagation index
  lit getEdgeLiteral(const lit s) const;
  // get bound literal from a bound propagation index
  lit getBoundLiteral(const lit s) const;
  // get explaantion lit from bound lit
  //    lit getLiteral(const lit l) const;
  // get literal index for explanation
  lit getBoundIndex(const lit b) const;
  // get the oldest literal that implies constraint c
  lit getImplicant(const BoundConstraint<T> &c) const;
  size_t getIndex(const var x) const;

  lit stamp(const lit l) const;

  const SparseSet<var> &getBranch() const;

  void post(Constraint *);
  void relax(Constraint *);
  BoundConstraint<T> getBound(const lit) const;
  DistanceConstraint<T> getEdge(const lit) const;
  lit getEdgeLit(std::pair<event, event>) const;
  bool hasEdge(std::pair<event, event>) const;
  void wake_me_on_event(const lit, const int);
  void wake_me_on_edge(const lit, const int);
  //    lit getTriggerLiteral(const int) const;

  void set(const var x, const boolean_state val,
           Explanation e = Constant::NoReason);
  void set(const lit p, Explanation e = Constant::NoReason);
  void set(const BoundConstraint<T> &c, Explanation e = Constant::NoReason);

  void setprimalBound(const T b);

  void trigger(const lit l);
  void propagate();

  void saveState();
  void restoreState(const int);
  void undo() override;
  template <typename S> void updatedualBound(S &objective);
  boolean_state search();
  void restart(const bool on_solution = false);
  //  void notifySolution();
  void backtrack(Explanation e);
  void branchRight();
  void learnConflict(Explanation e);
  bool hasSolution() const;

  void initializeSearch();
  boolean_state satisfiable();
  template <typename S> void optimize(S &objective);
  template <typename S> void optimize_dichotomy(S &objective);

  int decisionLevel(const genlit l) const;
  //    void incrementActivity(const var x);

  void getCriticalPath(std::vector<genlit> &clause);

  // explanation stuff
  void xplain(const lit l, const hint h, std::vector<lit> &Cl) override;
  std::ostream &print_reason(std::ostream &, const hint) const override;
  int getType() const override;

  std::ostream &display(std::ostream &os, const bool dom = true,
                        const bool bra = true, const bool sva = true,
                        const bool cla = false, const bool egr = false,
                        const bool vgr = false, const bool con = false) const;
  std::ostream &displayBranch(std::ostream &os) const;
  //  std::ostream &displayStats(std::ostream &os, const char *msg) const;
  std::ostream &displayVariables(std::ostream &os) const;
  std::ostream &displayConstraints(std::ostream &os) const;
    std::ostream &displayBound(std::ostream &os, BoundConstraint<T> c) const {return domain.bounds.displayBound(os, c);}
  std::ostream &displayProgress(std::ostream &os) const;
  std::ostream &displayHeader(std::ostream &os, const int width = 59) const;
  std::ostream &displaySummary(std::ostream &os, std::string msg) const;

    std::string prettyEvent(const event e) const;
  std::string prettyLiteral(const genlit el) const;

#ifdef DBG_TRACE
  void printTrace() const;
#endif

  BacktrackEnvironment &getEnv() { return env; }
  const Options &getOptions() const { return options; }

  /**
   * @name subscribable events
   */
  ///@{
  mutable SubscribableEvent<const std::vector<genlit> &>
      ClauseAdded; ///< triggered when a new clause is learned
  mutable SubscribableEvent<Explanation &>
      ConflictEncountered; ///< triggered when a conflict is encountered
  mutable SubscribableEvent<> SearchRestarted; ///< triggered on restart
  //    mutable SubscribableEvent<T, T, std::function<T(event, event)>,
  //    std::size_t> SolutionFound; ///< triggered when a solution is found
  //    mutable SubscribableEvent<std::function<T(event, event)>, std::size_t>
  //    SignificantSubstepFound; ///< triggered when an interesting partial
  //    solution has been found
  ///@}
  ///

  //    const std::vector<std::vector<Arc>> & getForwardGraph() const;
  //    const std::vector<std::vector<Arc>> & getBackwardGraph() const;

  void setActivityMap(heuristics::impl::EventActivityMap<T> *map) {
    activityMap = map;
  }
  heuristics::impl::EventActivityMap<T> *getActivityMap() {
    return activityMap;
  }
  double looseness(const DistanceConstraint<T> &c) const {
    // to - from <= distance
    return static_cast<double>(upper(c.from) - lower(c.to) + c.distance) /
           static_cast<double>(upper(HORIZON));
  }
  double looseness(const BoundConstraint<T> &c) const {
    auto x{EVENT(c.l)};
    if (SIGN(c.l) == LOWER) {
      // x >= distance
      return static_cast<double>(upper(x) + c.distance) /
             static_cast<double>(upper(HORIZON));
    } else {
      // x <= distance
      return static_cast<double>(c.distance - lower(x)) /
             static_cast<double>(upper(HORIZON));
    }
  }

  void clearLearnedClauses() { clauses.forgetAll(); }

  Task<T> &getTask(const event e);
    
    void resize(const size_t n);

private:
  Options options;

  BacktrackEnvironment env;

  Task<T> schedule;

  /*
  // the current temporal graph (manages its own reversibility), used to
  compute the bounds, i.e., the shortest paths from an to the origin
  DirectedGraph<StampedLabeledEdge<T>> core;

  // the following two vectors are indexed with dualBound(e)/primalBound(e)
  // the current bounds (sortest path)
  std::vector<T> bound;
  std::vector<lit> bound_index;

  // all literals 1 bit for type (bound/edge),
  //  * if edge: 1 bit for pos/neg and the rest for the variable
  //  * if bound: 1 bit
  std::vector<anylit> trail;
  std::vector<lit> prev;
  std::vector<Explanation> reason;

  */

public:
  TemporalNetwork<T> domain;

private:
  ClauseBase<T> clauses;

  // edges[l] with l=2*i+t is a distance constraint x-y <= k (for t=0) and y-x
  // <= k' this->satisfied(l) this->falsified(l)
  std::vector<DistanceConstraint<T>> edges;

  std::map<std::pair<event, event>, lit> edge_map;

//  std::vector<T> min_duration;
//  std::vector<T> max_duration;

  std::vector<Constraint *> constraints;
  //    std::vector<int> task_cons_id; // the "id" of the constraint (among
  //    constraints on tasks)

  DirectedGraph<int> evt_constraint_network;
  DirectedGraph<int> var_constraint_network;

  ConstraintQueue<3> propagation_queue;

  SparseSet<var> search_vars;
  Reversible<size_t> edge_propag_pointer;
  std::vector<bool> polarity;

  std::vector<bool> best_solution;

  Reversible<size_t> bound_propag_pointer;
  std::vector<Explanation> reason;

  std::vector<int> var_level;
  //    std::vector<double> activity;

  std::vector<Task<T>> tasks;
  std::vector<size_t> event2task_map; // for each event, the task it belongs to
                                      // (and '-1' if it's independent)
  //    std::vector<std::vector<size_t>> task2rank_map; // for each task j, and
  //    for the constraint of task_id i

  void initialize_baseline();

  //

  //    T dualBound(const event x) const;
  //    T primalBound(const event x) const;

  std::optional<heuristics::HeuristicManager<T>> heuristic;
  std::optional<heuristics::ValueHeuristicsManager> valueHeuristic;
  RestartPolicy *restart_policy = nullptr;
  unsigned int restart_limit{static_cast<unsigned int>(-1)};

  // parameters
  double weight_unit{1};
  double decay{1};
  double guard{1.0e+25};
  int init_level{0};

  // buffers
  std::vector<lit> conflict_edges;
  std::vector<lit> conflict_bounds;
  std::vector<genlit> conflict_set;
  std::vector<genlit> conflict;
  std::vector<bool> visited_edge;

  //    T lb{0};
  //    T ub{INFTY};

  double start_time;

  Scheduler<T> *baseline;
  heuristics::impl::EventActivityMap<T> *activityMap{NULL};

public:
  long unsigned int num_fails{0};
  long unsigned int num_choicepoints{0};
  long unsigned int num_backtracks{0};
  long unsigned int num_literals{0};
  long unsigned int num_updates{0};
  long unsigned int num_clause_propagations{0};
  long unsigned int num_cons_propagations{0};
  long unsigned int num_edge_prunings{0};
  long unsigned int num_cons_prunings{0};
  long unsigned int learnt_size{0};
  //    double gap_ratio{0.0};
  //    long unsigned int num_tight{0};

  lit endEdgeLit() const;
  lit endBoundLit() const;

  bool is_younger(const genlit a, const genlit b) const {
    if (LTYPE(a) == EDGE_LIT) {
      if (LTYPE(b) == EDGE_LIT) {
        return getIndex(VAR(FROM_GEN(a))) > getIndex(VAR(FROM_GEN(b)));
      } else {
        return stamp(FROM_GEN(b)) <
               static_cast<lit>(getIndex(VAR(FROM_GEN(a))));
      }
    } else {
      if (LTYPE(b) == EDGE_LIT) {
        return stamp(FROM_GEN(a)) >=
               static_cast<lit>(getIndex(VAR(FROM_GEN(b))));
      } else {
        return a > b;
      }
    }
  }

  std::function<bool(const genlit, const genlit)> younger =
      [&](const genlit a, const genlit b) { return this->is_younger(a, b); };

private:
  bool stoppingCondition() const;
  bool needExplanation(const genlit l) const;
  bool groundFact(const genlit l) const;
  bool inExplanation(const genlit l) const;
  void markVisited(const genlit l);
  Explanation getExplanation(const genlit l) const;
  bool isResponsibleForEdge(const lit l, const Constraint *c) const;
  bool isResponsibleForBound(const lit l, const Constraint *c) const;

  void analyze(Explanation e);
  //  void quickxplain();
  //  void greedy_minimization();
  void minimization(const size_t w);
  SparseSet<int> cons;
  std::vector<int> necessary;
  void resolve(const lit l, Explanation e);
  //    bool relevant(const lit l);
  bool relevant(const lit l, const size_t w);
  void clearVisited();

#ifdef DBG_CL

  //  TemporalNetwork<T> cut;

  std::ofstream *cl_file{NULL};
  int num_clauses{0};

  void writeLiteral(const lit l) const;
  void writeConstraint(const BoundConstraint<T> &c) const;
#endif

#ifdef DBG_SOL
  std::vector<bool> ref_solution;
  var is_on_track() {

    //      assert(ref_solution.empty());

    if (ref_solution.empty())
      return NoVar;
    auto n{static_cast<var>(numVariable())};
    for (var x{0}; x < n; ++x) {
      if (not isUndefined(x) and (value(x) != ref_solution[x]))
        return x;
    }
    return NoVar;
  }

public:
  void load(std::vector<bool> &sol) { ref_solution = sol; }
#endif
};

/// SCHEDULER
///
template <typename T>
Scheduler<T>::Scheduler(Options opt)
    : ReversibleObject(&env), options(std::move(opt)), schedule(*this),
      domain(*this), clauses(*this), evt_constraint_network(&env),
      var_constraint_network(&env), propagation_queue(constraints),
      edge_propag_pointer(0, &env), bound_propag_pointer(0, &env) {
  resize(2);
          domain.bounds.setLabel(0, "src");
          domain.bounds.setLabel(1, "cmax");
          
  //          tasks.push_back(schedule);
  //          std::cout << *this << std::endl;

  if (options.restart_policy == "luby") {
    restart_policy = new Luby(options.restart_base);
  } else if (options.restart_policy == "geom") {
    restart_policy =
        new Geometric(options.restart_base, options.restart_factor);
  } else {
    restart_policy = new NoRestart();
  }

#ifdef DBG_CL
  if (options.dbg_file != "")
    cl_file = new std::ofstream(options.dbg_file, std::ofstream::out);
#endif
}

template<typename T>
Scheduler<T>::~Scheduler() {
  for (auto c : constraints)
    delete c;
  delete restart_policy;
#ifdef DBG_CL
  if (cl_file != NULL)
    delete cl_file;
#endif
}

template <typename T> Task<T> &Scheduler<T>::getTask(const event e) {
  assert(event2task_map[e] < numTask());
  return tasks[event2task_map[e]];
}

template<typename T>
size_t Scheduler<T>::numTask() const {

  //    std::cout << min_duration.size() << " / " << tasks.size() << std::endl;

//  assert(min_duration.size() == tasks.size() or
//         min_duration.size() == (tasks.size() + 1));

  return tasks.size();
}

template<typename T>
size_t Scheduler<T>::numEvent() const {
  return domain.size();
}

template <typename T> size_t Scheduler<T>::numVariable() const {
  return search_vars.capacity();
}

template<typename T>
size_t Scheduler<T>::numClause() const {
  return clauses.size();
}

template<typename T>
size_t Scheduler<T>::numConstraint() const {
  return constraints.size();
}

template<typename T>
size_t Scheduler<T>::numBoundLiteral() const {
  return domain.bounds.numLiteral();
}

template <typename T> size_t Scheduler<T>::numEdgeLiteral() const {
  return search_vars.frontsize();
}

template<typename T>
lit Scheduler<T>::endEdgeLit() const {
  return static_cast<lit>(2 * numVariable());
}

template<typename T>
lit Scheduler<T>::endBoundLit() const {
  return static_cast<lit>(2 * numEvent());
}

template<typename T>
lit Scheduler<T>::stamp(const lit l) const {
  return domain.bounds.getStamp(l);
}

template<typename T>
const SparseSet<>& Scheduler<T>::getBranch() const {
//   return branch;
return search_vars;
}

// template<typename T>
// const std::vector<std::vector<Arc>> & Scheduler<T>::getForwardGraph() const {
//     return domain.getForwardGraph();
// }
//
// template<typename T>
// const std::vector<std::vector<Arc>> & Scheduler<T>::getBackwardGraph() const
// {
//     return domain.getBackwardGraph();
// }

template<typename T>
void Scheduler<T>::resize(const size_t n) {

//    assert(env.level() == 0);

//    std::cout << 11 << std::endl;

event2task_map.resize(n, static_cast<size_t>(-1));
domain.resize(n);
//    front_propag_pointer = 0;
edge_propag_pointer = 0;

//    visited_bound[LOWER].resize(n, INFTY);
//    visited_bound[UPPER].resize(n, INFTY);
}

template <typename T> void Scheduler<T>::load(ProblemInstance &data) {
  //    resize(data.durations.size());
  for (auto d : data.durations) {
    newTask(d, d);
    //        task ti{static_cast<task>(tasks.size())};
    //        tasks.emplace_back(*this, START(ti), END(ti), d, d);
    addConstraint(tasks.back().start.after(schedule.start));
    addConstraint(tasks.back().end.before(schedule.end));
  }
  //  for (auto [x, y, k] : data.constraints) {
  ////      // HACK
  ////      auto c{tasks[TASK(x)].end.before(tasks[TASK(y)].start, k)};
  ////
  ////      std::cout << c << std::endl;
  ////
  //////      addConstraint(c);
  ////
  ////
  ////      std::cout << prettyEvent(y) << " " << prettyEvent(x) << " " << -k <<
  /// std::endl;
  ////
  //
  //    newPrecedence(x, y, k);
  //  }

  std::vector<var> scope;
  std::vector<var> task_ids;
  std::vector<Task<int> *> the_tasks;
  for (auto &job : data.resources) {
    for (size_t i{0}; i < job.size(); ++i) {
      task_ids.push_back(job[i]);
      the_tasks.push_back(&tasks[job[i]]);
      for (size_t j{i + 1}; j < job.size(); ++j) {
          scope.push_back(newVariable(
                                      tasks[job[i]].start.after(tasks[job[j]].end),
                                      tasks[job[j]].start.after(tasks[job[i]].end)));
      }
    }
    if (options.edge_finding) {
      postEdgeFinding(
                      the_tasks.begin(),
                      the_tasks.end(), scope.begin()
                      );
    }
    if (options.transitivity) {
      postTransitivityReasoning(the_tasks.begin(),
                                the_tasks.end(), 
                                scope.begin());
    }
    the_tasks.clear();
    task_ids.clear();
    scope.clear();
  }

//    std::cout << *this << std::endl;
  //    exit(1);
}

template <typename T>
void Scheduler<T>::addConstraint(const DistanceConstraint<T> &c) {
  domain.newEdge(c.from, c.to, c.distance);
}

// [y - x <= d]
template<typename T>
void Scheduler<T>::newMaximumLag(event x, event y, T d) {

//    assert(env.level() == 0);

#ifdef DBG_TRACE
if (DBG_BOUND and (DBG_TRACE & PROPAGATION)) {

  if (x == ORIGIN) {
    BoundConstraint<T> c(LIT(y, UPPER), d);
    std::cout << " add " << c << std::endl;
  } else if (y == ORIGIN) {
    BoundConstraint<T> c(LIT(x, LOWER), d);
    std::cout << " add " << c << std::endl;
  } else {
    DistanceConstraint<T> e(x, y, d);
    std::cout << " add " << e << std::endl;
  }
}
#endif

if (x == ORIGIN)
  domain.newBound({LIT(y, UPPER), d});
else if (y == ORIGIN)
  domain.newBound({LIT(x, LOWER), d});
else if (edge_map.contains({x, y})) {
  auto l{edge_map.at({x, y})};
  assert(d == edges[l].distance);
  set(l);
} else {
  domain.newEdge(x, y, d);
}
}

// [x + d <= y] <-> [x - y <= -d-1]
template<typename T>
void Scheduler<T>::newPrecedence(event x, event y, T d) {
  newMaximumLag(y, x, -d);
}

//template <typename T> task Scheduler<T>::newEvent() {
//
//  auto ei{static_cast<event>(numEvent())};
//  resize(numEvent() + 1);
//
//  if ((ei % 2) == 0) {
//    min_duration.push_back(0);
//    max_duration.push_back(INFTY);
//  }
//
//  return ei;
//}

template <typename T>
task Scheduler<T>::newTask(const T min_dur, const T max_dur) {
  assert(env.level() == 0);
  assert(numEvent() >= 2);

//  if (2 * (numTask() + 1) != numEvent()) {
//    std::cout << "not ok\n";
//  }
  task ti{static_cast<task>(numTask())};
    tasks.emplace_back(*this, min_dur, max_dur);
    
    domain.bounds.setLabel(tasks.back().start.id(), "s"+std::to_string(tasks.back().id()));
    if(tasks.back().start.id() != tasks.back().getEnd())
        domain.bounds.setLabel(tasks.back().getStart(), "e"+std::to_string(tasks.back().id()));
    

//  newMaximumLag(START(ti), END(ti), max_dur);
//  newMaximumLag(END(ti), START(ti), -min_dur);
    
//    newMaximumLag(START(ti), END(ti), max_dur);
//    newMaximumLag(END(ti), START(ti), -min_dur);

//  min_duration.push_back(min_dur);
//  max_duration.push_back(max_dur);

  
  event2task_map[tasks.back().getStart()] =
      event2task_map[tasks.back().getEnd()] = (tasks.size() - 1);

  return ti;
}

template <typename T>
task Scheduler<T>::newTask(const T release, const T deadline, const T min_dur, const T max_dur) {
  assert(release + max_dur <= deadline);

  task k = newTask(min_dur, max_dur);

  // Minimal starting date
//  newPrecedence(ORIGIN, START(tk), release);
    set(tasks[k].start.after(release));
  // Maximal ending date
    set(tasks[k].end.before(deadline));
//  newMaximumLag(ORIGIN, END(tk), deadline);

  return k;
}

template <typename T>
var Scheduler<T>::newVariable(const DistanceConstraint<T> &if_true,
                              const DistanceConstraint<T> &if_false) {
  assert(env.level() == 0);

  var x{static_cast<var>(numVariable())};

  if (if_true != Constant::NoEdge<T>) {
    std::pair<event, event> ptrue{if_true.from, if_true.to};
    edge_map[ptrue] = static_cast<lit>(edges.size());
  }
  edges.emplace_back(if_true);

  if (if_false != Constant::NoEdge<T>) {
    std::pair<event, event> pfalse{if_false.from, if_false.to};
    edge_map[pfalse] = static_cast<lit>(edges.size());
    edges.emplace_back(if_false);
  } else if (if_true != Constant::NoEdge<T>) {
    edges.emplace_back(~if_true);
  } else {
    edges.emplace_back(if_false);
  }

  search_vars.resize(edges.size() / 2);
  polarity.resize(edges.size() / 2, false);

  //    branch.resize(edges.size()/2);
  //    back_propag_pointer = branch.size();

  clauses.resize(numEvent(), numVariable());

  //    var_level.resize(branch.size(),-1);
  //    activity.resize(branch.size(), 0);
  //
  //    reason.resize(branch.size(), Constant::NoReason);
  //    visited_edge.resize(branch.size(), false);

  var_level.resize(numVariable(), -1);
  //    activity.resize(numVariable(), 0);

  reason.resize(numVariable(), Constant::NoReason);
  visited_edge.resize(numVariable(), false);

  if (if_true != Constant::NoEdge<T>) {
    post(new EdgeConstraint<T>(*this, POS(x)));
    post(new EdgeConstraint<T>(*this, NEG(x)));
  }

  return x;
}

//template<typename T>
//void Scheduler<T>::newVariable(const event x, const event y, const T d) {
//    
//    assert(env.level() == 0);
//    
//    edges.emplace_back(x,y,d);
//    branch.resize(edges.size());
//    back_propag_pointer = branch.size();
//    
//    clauses.resize(edges.size());
//    
//    var_level.resize(edges.size(),-1);
//    activity.resize(edges.size(), 0);
//    
//    reason.resize(edges.size(), Constant::NoReason);
//
//}

//template<typename T>
//void Scheduler<T>::newDisjunct(const event x1, const event y1, const T d1, const event x2, const event y2, const T d2) {
//    var x{static_cast<var>(numVariable())};
//    newVariable(x1, y1, d1);
//    
//    var y{static_cast<var>(numVariable())};
//    newVariable(x2, y2, d2);
//    
//    lit_buffer.clear();
//    lit_buffer.push_back(POS(x));
//    lit_buffer.push_back(POS(y));
//    
//    clauses.add_to(lit_buffer.begin(), lit_buffer.end(), clauses.base, numClause());
//    
//    lit_buffer.clear();
//    lit_buffer.push_back(NEG(x));
//    lit_buffer.push_back(NEG(y));
//    
//    clauses.add_to(lit_buffer.begin(), lit_buffer.end(), clauses.base, numClause());
//    
//    post(new EdgeConstraint<T>(*this, POS(x)));
//    post(new EdgeConstraint<T>(*this, POS(y)));
//    post(new EdgeConstraint<T>(*this, NEG(x)));
//    post(new EdgeConstraint<T>(*this, NEG(y)));
//}

template<typename T> bool Scheduler<T>::value(const var x) const {
  assert(not isUndefined(x));

  //    assert(branch.isfront(x) == polarity[x]);

  //    return branch.isfront(x);

  return polarity[x];
}

template<typename T> bool Scheduler<T>::isTrue(const var x) const {

//    assert(branch.isfront(x) == (polarity[x] and not search_vars.has(x)));

//    return branch.isfront(x);
return (polarity[x] and not search_vars.has(x));
}

template<typename T> bool Scheduler<T>::isFalse(const var x) const {
//    assert(branch.isback(x) == not(polarity[x] or search_vars.has(x)));

//    return branch.isback(x);

return not(polarity[x] or search_vars.has(x));
}

template<typename T> bool Scheduler<T>::isUndefined(const var x) const {
//    assert(branch.has(x) == search_vars.has(x));

//    return branch.has(x);

return search_vars.has(x);
}

template<typename T>
bool Scheduler<T>::falsified(const lit l) const {
  if (SIGN(l))
    return isFalse(VAR(l));
  return isTrue(VAR(l));
}

template<typename T>
bool Scheduler<T>::satisfied(const lit l) const {
  if (SIGN(l))
    return isTrue(VAR(l));
  return isFalse(VAR(l));
  //    auto lt{LTYPE(gl)};
  //    lit l{FROM_GEN(gl)};
  //
  //    if(lt == EDGE_LIT) {
  //        if(SIGN(l))
  //            return isTrue(VAR(l));
  //        return isFalse(VAR(l));
  //    }
  //    return domain.bounds.satisfied(l);
}

template<typename T>
bool Scheduler<T>::falsified(const BoundConstraint<T>& c) const {
  return domain.bounds.falsified(c);
}

template<typename T>
bool Scheduler<T>::satisfied(const BoundConstraint<T>& c) const {
  return domain.bounds.satisfied(c);
}

template<typename T>
T Scheduler<T>::upper(const event e) const {
  return domain.bounds.upper()[e];
}

template<typename T>
T Scheduler<T>::lower(const event e) const {
  return -domain.bounds.lower()[e];
}

template<typename T>
T Scheduler<T>::bound(const bool s, const event e) const {
  return domain.bounds.get(s, e);
}

//template <typename T> T Scheduler<T>::minDuration(const task t) const {
//  return min_duration[t];
//}
//
//template <typename T> T Scheduler<T>::maxDuration(const task t) const {
//  return max_duration[t];
//}

template <typename T>
T Scheduler<T>::distance(const event x, const event y) const {
  T ub{upper(y) - lower(x)};
  std::pair<event, event> e{x, y};
  if (hasEdge(e)) {
    return std::min(getEdge(getEdgeLit(e)).distance, ub);
  }
  return ub;
}

// template <typename T> T Scheduler<T>::getMakespan() const { return ub; }

template <typename T>
auto Scheduler<T>::getSolution() const noexcept -> const std::vector<bool> & {
  return best_solution;
}

template<typename T>
lit Scheduler<T>::getEdgeLiteral(const lit s) const {

//    assert(propagation_stack[s] == LIT(search_vars[s], polarity[search_vars[s]]));

//    return propagation_stack[s];

auto x{search_vars[s]};
return LIT(x, polarity[x]);
}

template<typename T>
lit Scheduler<T>::getBoundLiteral(const lit s) const {
  return domain.bounds.getLiteral(s);
}

// template<typename T>
// lit Scheduler<T>::getExplanation(const lit l) const {
//   return BOUND(domain.bounds.getIndex(l));
// }

template<typename T>
lit Scheduler<T>::getBoundIndex(const lit l) const {
  return domain.bounds.getIndex(l);
}

template<typename T>
lit Scheduler<T>::getImplicant(const BoundConstraint<T>& c) const {
  return domain.bounds.getImplicant(c);
}

template<typename T>
size_t Scheduler<T>::getIndex(const var x) const {
  return search_vars.index(x);
}

template<typename T>
DistanceConstraint<T> Scheduler<T>::getEdge(const lit l) const {
  return edges[l];
}

template<typename T>
BoundConstraint<T> Scheduler<T>::getBound(const lit l) const {
  return domain.bounds.getConstraint(l);
}

template<typename T>
lit Scheduler<T>::getEdgeLit(std::pair<event,event> arc) const {
  return edge_map.at(arc);
}

template<typename T>
bool Scheduler<T>::hasEdge(std::pair<event,event> arc) const {
  return edge_map.contains(arc);
}

template<typename T>
void Scheduler<T>::wake_me_on_event(const lit l, const int c) {
  evt_constraint_network.add(l, c);
}

template<typename T>
void Scheduler<T>::wake_me_on_edge(const lit l, const int c) {
  var_constraint_network.add(l, c);
}

template <typename T>
template <typename Container>
void Scheduler<T>::addClause(Container &c) {
  clauses.add(c.begin(), c.end());
}

template <typename T>
template <typename ItTask, typename ItVar>
void Scheduler<T>::postEdgeFinding(const ItTask beg_task, const ItTask end_task,
                                   const ItVar beg_var
                                   ) {
  post(new DisjunctiveEdgeFinding(*this, beg_task, end_task,
                                  beg_var
                                  ));
}

template <typename T>
template <typename ItTask, typename ItVar>
void Scheduler<T>::postTransitivityReasoning(const ItTask beg_task,
                                             const ItTask end_task,
                                             const ItVar beg_var) {
  post(new Transitivity(*this, beg_task, end_task, beg_var));
}

template <typename T> void Scheduler<T>::post(Constraint *con) {

  constraints.push_back(con);
  propagation_queue.resize(constraints.size());

  var_constraint_network.resize(std::max(2 * numVariable(), numConstraint()));
  evt_constraint_network.resize(std::max(2 * numEvent(), numConstraint()));

  con->post(numConstraint() - 1);
}

template<typename T>
void Scheduler<T>::relax(Constraint *con) {

//    constraints.push_back(con);
//    propagation_queue.resize(constraints.size());
//    
//    var_constraint_network.resize(std::max(2*numVariable(), numConstraint()));
//    evt_constraint_network.resize(std::max(2*numEvent(), numConstraint()));
//    
//    con->post(numConstraint() - 1);

//    std::cout << "relax " << *con << std::endl;

if (var_constraint_network.indegree(con->id()) > 0)
  var_constraint_network.remove(con->id(), IN);
if (evt_constraint_network.indegree(con->id()) > 0)
  evt_constraint_network.remove(con->id(), IN);

//    displayConstraints(std::cout)
}

template <typename T> void Scheduler<T>::setprimalBound(const T b) {
  //    assert(env.level() == 0);
  newMaximumLag(ORIGIN, HORIZON, b);
  //    domain.newEdge(ORIGIN, HORIZON, b);
}

template<typename T>
void Scheduler<T>::set(const BoundConstraint<T>& c, Explanation e) {
  domain.newBound(c, e);
}

template<typename T>
void Scheduler<T>::set(const var x, const boolean_state val, Explanation e) {
  set(LIT(x, val), e);
}

template<typename T>
void Scheduler<T>::set(const lit l, Explanation e) {

#ifdef DBG_TRACE
  if (DBG_BOUND and (DBG_TRACE & PROPAGATION)) {
    std::cout << "new literal: " << prettyLiteral(EDGE(l));
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
      std::cout << " b/c " << e << " (" << e.expl->id() << ")";
    }
    std::cout << std::endl;
  }
#endif

  var x{VAR(l)};
  auto val{SIGN(l)};
  if (isUndefined(x)) {
    search_vars.remove_front(x);
    polarity[x] = val;

    reason[x] = e;
    var_level[x] = env.level();

    if (edges[l] != Constant::NoEdge<T>)
      domain.newEdge(edges[l].from, edges[l].to, edges[l].distance);

  } else if (value(x) != val) {
    throw Failure(e);
  }
}

//template<typename T>
//void Scheduler<T>::set(const lit p, Explanation e) {
//    set(VAR(p), SIGN(p), e);
//}

//template<typename T>
//void Scheduler<T>::newVariables(const std::vector<DistanceConstraint<T>>& vars) {
//    
//    assert(env.level() == 0);
//    
//    for(auto e : vars)
//        edges.emplace_back(e);
//    domain.resize(edges.size());
//    back_propag_pointer = edges.size();
//}

template<typename T>
void tempo::Scheduler<T>::trigger(const lit l) {

  auto r{reason[VAR(l)].expl};

#ifdef DBG_TRACE
  if (DBG_TRACE & QUEUE) {
    std::cout << "triggers for (" << l << ") " << edges[l] << " b/c " << r->id()
              << "/" << reason[VAR(l)] << std::endl;
  }
#endif

  clauses.unit_propagate(EDGE(NOT(l)));

  const std::vector<int> &cons = var_constraint_network[l];
  const std::vector<unsigned> &rank = var_constraint_network.rank(l);

  //    std::cout <<reason[VAR(l)] << std::endl;

  for (auto i{cons.size()}; i-- > 0;) {
    //    for (unsigned i{0}; i < cons.size(); ++i) {

#ifdef DBG_TRACE
    if (DBG_TRACE & QUEUE) {
      std::cout << " -" << *(constraints[cons[i]]) << " ("
                << constraints[cons[i]]->id() << ")" << std::endl;
    }
#endif

    propagation_queue.edge_triggers(l, rank[i], cons[i], r);
  }
}

// template<typename T>
// void tempo::Scheduler<T>::propagate() {
////#ifdef DBG_SOL
////    var x;
////
////#endif
//
//    size_t edge_pointer{static_cast<size_t>(edge_propag_pointer)};
//    size_t bound_pointer{static_cast<size_t>(bound_propag_pointer)};
//  while (not propagation_queue.empty() or
//         (search_vars.frontsize() > edge_propag_pointer) or
//         (numBoundLiteral() > static_cast<size_t>(bound_propag_pointer))) {
////         (search_vars.frontsize() > edge_pointer) or
////         (numBoundLiteral() > static_cast<size_t>(bound_pointer))) {
//
//    while (search_vars.frontsize() > edge_propag_pointer) {
////      while (search_vars.frontsize() > edge_pointer) {
//      ++num_literals;
//      lit l{LIT(search_vars[edge_propag_pointer],
//                polarity[search_vars[edge_propag_pointer]])};
//      //
//      //            std::cout << " prop: (" <<
//      //            LIT(search_vars[edge_propag_pointer],
//      //            polarity[edge_propag_pointer]) << ") [" <<
//      //            edges[LIT(search_vars[edge_propag_pointer],
//      //            polarity[search_vars[edge_propag_pointer]])] << "] <" <<
//      //            polarity[search_vars[edge_propag_pointer]] << ">\n";
//
////#ifdef DBG_SOL
////      x = is_on_track();
////      if (x != NoVar) {
////        std::cout << "bug before trigger " << edges[l] << std::endl;
////        exit(1);
////      }
////
////#endif
//
//      trigger(l);
//
////#ifdef DBG_SOL
////      x = is_on_track();
////      if (x != NoVar) {
////        std::cout << "bug after trigger " << edges[l] << std::endl;
////        exit(1);
////      }
////
////#endif
//
//      ++edge_propag_pointer;
//    }
//
//    while (numBoundLiteral() > bound_propag_pointer) {
//      ++num_literals;
//      lit l{getBoundLiteral(bound_propag_pointer)};
//
////#ifdef DBG_SOL
////      x = is_on_track();
////      if (x != NoVar) {
////        std::cout << "bug before trigger " <<
///prettyLiteral(BOUND(bound_propag_pointer)) << std::endl; /        exit(1); /
///}
////
////#endif
//
//      const std::vector<int> &cons = evt_constraint_network[l];
//      const std::vector<unsigned> &rank = evt_constraint_network.rank(l);
//
//      auto r{domain.bounds.getExplanation(bound_propag_pointer).expl};
//
//#ifdef DBG_TRACE
//      if (DBG_TRACE & QUEUE) {
//        std::cout << "triggers for " << prettyEventLit(l) << " b/c " <<
//        r->id()
//                  << std::endl;
//      }
//#endif
//
//      //        std::cout << domain.bounds.getExplanation(l) << std::endl;
//
//      // important to visit in reverse order to be robust to relax
//      for (auto i{cons.size()}; i-- > 0;) {
//
//#ifdef DBG_TRACE
//        if (DBG_TRACE & QUEUE) {
//          std::cout << " -" << *(constraints[cons[i]]) << " ("
//                    << constraints[cons[i]]->id() << ")" << std::endl;
//        }
//#endif
//
//        propagation_queue.bound_triggers(l, rank[i], cons[i], r);
//      }
//
//
////#ifdef DBG_SOL
////      x = is_on_track();
////      if (x != NoVar) {
////        std::cout << "bug after trigger " <<
///prettyLiteral(BOUND(bound_propag_pointer)) << std::endl; /        exit(1); /
///}
////
////#endif
//
//      ++bound_propag_pointer;
//    }
//
//    if (not propagation_queue.empty()) {
//      auto cons{propagation_queue.pop_front()};
//
//#ifdef DBG_TRACE
//      if (DBG_TRACE & QUEUE) {
//        std::cout << "propagate " << *cons << std::endl;
//      }
//#endif
//
//
//        ++num_cons_propagations;
//
////#ifdef DBG_SOL
////        x = is_on_track();
////      if (x != NoVar) {
////        std::cout << "bug before propagation " << num_cons_propagations
////          << " on var " << edges[POS(x)] << " <> " << edges[NEG(x)]
////                  << std::endl;
////        exit(1);
////      }
////#endif
//
//      cons->propagate();
//
//
////#ifdef DBG_SOL
////        x = is_on_track();
////      if (x != NoVar) {
////        std::cout << "bug after propagation " << num_cons_propagations
////          << " on var " << edges[POS(x)] << " <> " << edges[NEG(x)]
////                  << std::endl;
////        exit(1);
////      }
////
////#endif
//
////#ifdef DBG_SOL
////      //        auto z{is_on_track()};
////      //        std::cout << was_on_track << " -> " << z << std::endl;
////      if (was_on_track and not is_on_track()) {
////        std::cout << "bug at propagation " << num_cons_propagations
////                  << std::endl;
////        exit(1);
////      }
////#endif
//    }
//  }
//}

template <typename T> void tempo::Scheduler<T>::propagate() {

  size_t edge_pointer{static_cast<size_t>(edge_propag_pointer)};
  size_t bound_pointer{static_cast<size_t>(bound_propag_pointer)};
  //    edge_pointer = static_cast<size_t>(edge_propag_pointer);
  //    bound_pointer = static_cast<size_t>(bound_propag_pointer);
  while (not propagation_queue.empty() or
         (search_vars.frontsize() > edge_pointer) or
         (numBoundLiteral() > static_cast<size_t>(bound_pointer))) {

    while (search_vars.frontsize() > edge_pointer) {
      ++num_literals;
      lit l{
          LIT(search_vars[edge_pointer], polarity[search_vars[edge_pointer]])};
      trigger(l);

      //      auto e{getEdge(l)};
      //      std::cout << "prop edge " << e << std::endl;
      //      std::cout << "dist(" << e.from << "," << e.to
      //                << ") = " << distance(e.from, e.to) << std::endl;
      //      std::cout << distance(END(12), START(7)) << std::endl <<
      //      std::endl;

      ++edge_pointer;
    }

    while (numBoundLiteral() > bound_pointer) {
      ++num_literals;
      lit l{getBoundLiteral(bound_pointer)};

      const std::vector<int> &cons = evt_constraint_network[l];
      const std::vector<unsigned> &rank = evt_constraint_network.rank(l);

      auto r{domain.bounds.getExplanation(bound_pointer).expl};

#ifdef DBG_TRACE
      if (DBG_TRACE & QUEUE) {
        std::cout << "triggers for " << prettyEventLit(l) << " b/c " << r->id()
                  << std::endl;
      }
#endif

      //        std::cout << domain.bounds.getExplanation(l) << std::endl;

      // important to visit in reverse order to be robust to relax
      for (auto i{cons.size()}; i-- > 0;) {

#ifdef DBG_TRACE
        if (DBG_TRACE & QUEUE) {
          std::cout << " -" << *(constraints[cons[i]]) << " ("
                    << constraints[cons[i]]->id() << ")" << std::endl;
        }
#endif

        propagation_queue.bound_triggers(l, rank[i], cons[i], r);
      }

      ++bound_pointer;
    }

    if (not propagation_queue.empty()) {
      auto cons{propagation_queue.pop_front()};

#ifdef DBG_TRACE
      if (DBG_TRACE & QUEUE) {
        std::cout << "propagate " << *cons << std::endl;
      }
#endif

      ++num_cons_propagations;
      cons->propagate();
    }
  }

  edge_propag_pointer = edge_pointer;
  bound_propag_pointer = bound_pointer;
}

template<typename T>
int Scheduler<T>::decisionLevel(const genlit l) const {
  if (LTYPE(l) == EDGE_LIT)
    return var_level[VAR(FROM_GEN(l))];
  lit bl{FROM_GEN(l)};
  auto idx{stamp(bl)};
  if (idx < 0)
    return init_level;
  auto x{search_vars[idx]};
  return var_level[x];
}

//template<typename T>
//void Scheduler<T>::incrementActivity(const var x) {
//    activity[x] += weight_unit;
//}

template<typename T>
void Scheduler<T>::getCriticalPath(std::vector<genlit>& clause) {
  domain.findExplanationPath(HORIZON, ORIGIN, clause,
                             static_cast<lit>(numEdgeLiteral()) - 1);
}

template<typename T>
void Scheduler<T>::xplain(const lit l, const hint h, std::vector<lit> &Cl) {
  if (l == NoLit) {
    //        std::cout << "domain exhausted on var [" << edges[POS(h)] <<
    //        "]<>[" << edges[NEG(h)] << "]";

    Cl.push_back(EDGE(POS(h)));
    Cl.push_back(EDGE(NEG(h)));
  } else {

    assert(LTYPE(l) == BOUND_LIT);

    lit bl{getBoundLiteral(FROM_GEN(l))};
    auto x{EVENT(bl)};

    if (SIGN(bl) == LOWER) {

      //            std::cout << "hint: " << prettyLiteral(EDGE(h)) <<
      //            std::endl;

      if (x != edges[h].from) {
        std::pair<event, event> arc{x, edges[h].from};

        if (edge_map.contains(arc)) {

          //                    std::cout << " *edge*\n" ;

          auto q{edge_map.at(arc)};

          assert(is_younger(l, EDGE(q)));

          Cl.push_back(EDGE(q));
        } else {

          //                    std::cout << " *path* (<= " <<
          //                    domain.bounds.getStamp(FROM_GEN(l)) << ")\n" ;

          domain.findExplanationPath(x, edges[h].from, Cl,
                                     static_cast<lit>(stamp(FROM_GEN(l))));
        }
      }

      auto q{domain.bounds.getPastIndex(LOWERBOUND(edges[h].to), FROM_GEN(l))};
      //            auto
      //            q{BOUND(domain.bounds.getIndex(dualBound(edges[h].to)))};

      if (q != NoLit) {
        assert(is_younger(l, BOUND(q)));

        Cl.push_back(BOUND(q));
      }
#ifdef DBG_TRACE
      else if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
        std::cout << " no need to explain lower bound of "
                  << prettyEvent(edges[h].to) << "\n";
      }
#endif

      assert(SIGN(getBoundLiteral(FROM_GEN(BOUND(
                 domain.bounds.getIndex(LOWERBOUND(edges[h].to)))))) == LOWER);

    } else {
      if (x != edges[h].to) {
        std::pair<event, event> arc{edges[h].to, x};

        if (edge_map.contains(arc)) {

          auto q{edge_map.at(arc)};

          assert(is_younger(l, EDGE(q)));

          Cl.push_back(EDGE(q));
        } else {

          //                    std::cout << " *path* (<= " <<
          //                    domain.bounds.getStamp(FROM_GEN(l)) << ")\n" ;

          domain.findExplanationPath(edges[h].to, x, Cl, stamp(FROM_GEN(l)));
        }
      }
      auto q{
          domain.bounds.getPastIndex(UPPERBOUND(edges[h].from), FROM_GEN(l))};
      //            auto
      //            q{BOUND(domain.bounds.getIndex(primalBound(edges[h].from)))};

      if (q != NoLit) {
        assert(is_younger(l, BOUND(q)));

        Cl.push_back(BOUND(q));
      }
#ifdef DBG_TRACE
      else if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
        std::cout << " no need to explain upper bound of "
                  << prettyEvent(edges[h].from) << "\n";
      }
#endif
    }

    assert(is_younger(l, EDGE(h)));
    Cl.push_back(EDGE(h));
  }

  /*TODO*/
  // h is the literal that lead to the bound change l
  //    DistanceConstraint<T> edge{edges[h]}; // (x,y)
  //
  //    BoundConstraint<T> b{bounds.getConstraint(l)}; // z
  //
  //    if(b.dis)
  //
  //    explainBound()
}

template <typename T>
std::ostream &Scheduler<T>::print_reason(std::ostream &os, const hint) const {
  os << "transitivity";
  return os;
}

template <typename T> int Scheduler<T>::getType() const { return PATHEXPL; }

template<typename T>
void Scheduler<T>::saveState() {
  env.save();
  ReversibleObject::save();
}

template<typename T>
void Scheduler<T>::restoreState(const int l) {
  env.restore(l);
}

template<typename T>
void Scheduler<T>::undo() {
  search_vars.setStart(edge_propag_pointer);
}

template<typename T>
void Scheduler<T>::branchRight() {

  restoreState(env.level() - 1);
  lit deduction_var{search_vars[edge_propag_pointer]};
  bool sign{not polarity[deduction_var]};

#ifdef DBG_TRACE
  if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
    std::cout << "-- backtrack to lvl " << env.level() << " & deduce "
              << edges[LIT(deduction_var, sign)] << " [i=" << num_choicepoints
              << "] --\n";
  }
#endif

  set(deduction_var, sign);
}

#ifdef DBG_CL
template<typename T>
void Scheduler<T>::writeLiteral(const lit l) const {
  auto j{FROM_GEN(l)};
  if (LTYPE(l) == EDGE_LIT) {
    *cl_file << " " << edges[j].from << " " << edges[j].to << " "
             << edges[j].distance;
  } else {
    auto c{domain.bounds.getConstraint(j)};
    writeConstraint(c);
  }
}

template<typename T>
void Scheduler<T>::writeConstraint(const BoundConstraint<T>& c) const {
  auto s{SIGN(c.l)};
  if (s == LOWER) {
    *cl_file << " " << EVENT(c.l) << " 0 " << c.distance;
  } else {
    *cl_file << " 0 " << EVENT(c.l) << " " << c.distance;
  }
}
#endif

template <typename T>
bool Scheduler<T>::relevant(const lit l, const size_t max_width) {
  //    assert(conflict_set.empty());
  conflict_set.clear();

#ifdef DBG_MINIMIZATION
  std::cout << "is " << prettyLiteral(l) << " relevant?\n";
#endif

  bool redundant{true};
  conflict_set.push_back(l);
  size_t next{0}, i, s, ve{conflict_edges.size()}, vb{conflict_bounds.size()};
  do {
    lit p{conflict_set[next++]};
    Explanation e{getExplanation(p)};

#ifdef DBG_MINIMIZATION
    std::cout << " resolve " << prettyLiteral(p) << " (" << e << ")\n";
#endif

    if (e != Constant::NoReason) {
      i = s = conflict_set.size();
      e.explain(p, conflict_set);
      for (; i != conflict_set.size(); ++i) {

#ifdef DBG_MINIMIZATION
        std::cout << "  - " << prettyLiteral(conflict_set[i]);
#endif

        //                if(getExplanation(p) == Constant

        if (not(groundFact(conflict_set[i]) or
                inExplanation(conflict_set[i]))) {

#ifdef DBG_MINIMIZATION
          std::cout << " new\n";
#endif

          if (conflict_set.size() >= (max_width + next)) {
            redundant = false;
            break;
          }
          markVisited(conflict_set[i]);
          conflict_set[s] = conflict_set[i];
          ++s;
        }
#ifdef DBG_MINIMIZATION
        else {
          if (inExplanation(conflict_set[i]))
            std::cout << " already in cut\n";
          else
            std::cout << " ground fact\n";
        }
#endif
      }
      conflict_set.resize(s);

#ifdef DBG_MINIMIZATION
      std::cout << " remain to explain:";
      for (i = next; i < conflict_set.size(); ++i) {
        std::cout << " " << prettyLiteral(conflict_set[i]);
      }
      std::cout << std::endl;
#endif

    } else {
      redundant = false;
#ifdef DBG_MINIMIZATION
      std::cout << " (decision)\n";
#endif
    }
  } while (conflict_set.size() > next and redundant);

  if (not redundant) {

#ifdef DBG_MINIMIZATION
    std::cout << "RELEVANT!\n";
#endif

    while (conflict_edges.size() > ve) {
      visited_edge[conflict_edges.back()] = false;
      conflict_edges.pop_back();
    }
    while (conflict_bounds.size() > vb) {
      domain.bounds.visited[conflict_bounds.back()] = false;
      conflict_bounds.pop_back();
    }

    return true;
  }

#ifdef DBG_MINIMIZATION
  std::cout << "NOT RELEVANT!\n";
#endif

  return false;
}

// template <typename T>
// bool Scheduler<T>::relevant(const lit l) {
//     Explanation e{getExplanation(l)};
//     if(e == Constant::NoReason) {
//         return true;
//     }
//
//     conflict_set.clear();
//     e.explain(l, conflict_set);
//
//     for(auto p : conflict_set) {
//         if(not (groundFact(p) or inExplanation(p))) {
//             return true;
//         }
//     }
//
//     return false;
// }

// template <typename T> void Scheduler<T>::greedy_minimization() {
//     necessary.clear();
//     for(auto l : conflict) {
//         if(relevant(l)) {
//             necessary.push_back(l);
//             assert(relevant(l,0));
//         } else {
//             assert(not relevant(l,0));
//         }
//     }
//     if (conflict.size() > necessary.size()) {
//       conflict = necessary;
//     }
// }

template <typename T> void Scheduler<T>::minimization(const size_t max_width) {
  necessary.clear();
  for (auto l : conflict) {
    if (relevant(l, max_width))
      necessary.push_back(l);
  }
  if (conflict.size() > necessary.size()) {
    conflict = necessary;
  }
}

// template <typename T> void Scheduler<T>::quickxplain(const T ub) {
//
//   baseline->setprimalBound(ub - Gap<T>::epsilon());
//
//   cons.reserve(conflict.size());
//   cons.clear();
//   for (size_t i{0}; i < conflict.size(); ++i) {
//     cons.add(i);
//   }
//   necessary.clear();
//
//   while (cons.size() > necessary.size()) {
//
//#ifdef DBG_MINIMIZATION
//     std::cout << std::endl
//               << *baseline << std::endl
//               << "new round " << cons.size() << "/" << necessary.size() <<
//               "\n";
//     std::cout << cons << std::endl;
//#endif
//
//     baseline->saveState();
//
//     for (auto i : necessary) {
//       auto l{conflict[i]};
//       //            auto pl{baseline->num_literals};
//
//#ifdef DBG_MINIMIZATION
//       std::cout << "-add necessary cons (" << i << ") ";
//       if (LTYPE(l) == EDGE_LIT)
//         std::cout << getEdge(FROM_GEN(l)) << std::endl;
//       else
//         std::cout << getBound(FROM_GEN(l)) << std::endl;
//#endif
//
//       try {
//         if (LTYPE(l) == EDGE_LIT)
//           baseline->set(FROM_GEN(l));
//         else
//           baseline->set(getBound(FROM_GEN(l)));
//         baseline->propagate();
//       } catch (Failure &f) {
//         cons.setStart(cons.end_idx());
//         break;
//       }
//       cons.remove_back(i);
//       //            assert(pl != baseline->num_literals);
//     }
//
//#ifdef DBG_MINIMIZATION
//     std::cout << cons << std::endl << *baseline << std::endl;
//#endif
//
//     while (not cons.empty()) {
//       auto i{cons.back()};
//       auto l{conflict[i]};
//
//       //            auto [x,y,k] = C[i];
//       auto pl{baseline->num_literals};
//
//#ifdef DBG_MINIMIZATION
//       std::cout << " -try (" << i << ") ";
//       if (LTYPE(l) == EDGE_LIT)
//         std::cout << getEdge(FROM_GEN(l)) << std::endl;
//       else
//         std::cout << getBound(FROM_GEN(l)) << std::endl;
//#endif
//
//       bool fail{false};
//       try {
//         if (LTYPE(l) == EDGE_LIT)
//           baseline->set(FROM_GEN(l));
//         else
//           baseline->set(getBound(FROM_GEN(l)));
//         baseline->propagate();
//       } catch (Failure &f) {
//         fail = true;
//       }
//       if (fail) {
//#ifdef DBG_MINIMIZATION
//         std::cout << " --> all the rest is subsumed (fail)!\n";
//#endif
//         necessary.push_back(i);
//         cons.remove_back(i);
//         while (not cons.empty()) {
//           cons.remove_front(cons.back());
//         }
//         //                    cons.setEnd(cons.start_idx());
//       } else if (pl == baseline->num_literals) {
//#ifdef DBG_MINIMIZATION
//         std::cout << " --> subsumed!\n";
//#endif
//         cons.remove_front(i);
//       } else {
//         if (cons.size() == 1) {
//#ifdef DBG_MINIMIZATION
//           std::cout << " --> necessary!\n";
//#endif
//           necessary.push_back(i);
//         }
//         cons.pop_back();
//       }
//#ifdef DBG_MINIMIZATION
//       std::cout << cons << std::endl;
//#endif
//     }
//
//     cons.setEnd(conflict.size());
//     baseline->restoreState(1);
//   }
//
//#ifdef DBG_MINIMIZATION
//   std::cout << "original cl:";
//   for (auto i : conflict) {
//     std::cout << " " << prettyLiteral(i);
//   }
//   std::cout << std::endl;
//   std::cout << "minimal cl:";
//   for (auto i : necessary) {
//     std::cout << " " << prettyLiteral(conflict[i]);
//   }
//   std::cout << std::endl;
//#endif
//
//
//
////    if(conflict.size() > necessary.size()) {
////        cons.clear();
////        for(auto i : necessary)
////        {
////            cons.add(i);
////        }
////        std::cout << std::endl;
////        for (auto i : conflict) {
////          std::cout << " " << prettyLiteral(i);
////            if(not cons.has(i))
////                std::cout << " (deleted)";
////            std::cout << std::endl;
////        }
////        std::cout << std::endl;
////    }
//
//
//  //    if(num_fails==6)
//  //        exit(1);
//
//  //    std::cout << conflict.size() << "/" << necessary.size() << std::endl;
//
//  conflict_set.clear();
//  for (auto i : necessary) {
//    //        std::cout << i << "/" << conflict.size() << std::endl;
//    conflict_set.push_back(conflict[i]);
//  }
//  conflict = conflict_set;
//
//  //    std::cout << *this << std::endl;
//  //  std::cout << *baseline << std::endl;
//
//  //  exit(1);
//
//  //    for(auto l : conflict) {
//  //
//  //    }
//}

template<typename T>
void Scheduler<T>::learnConflict(Explanation e) {

  analyze(e);

  //  if (options.minimization >= 1000) {
  //    quickxplain();
  //  } else
  if (options.minimization >= 0) {
    minimization(static_cast<size_t>(options.minimization));
  }

  clearVisited();

  std::sort(conflict.begin(), conflict.end(), [&](const lit a, const lit b) {
    return decisionLevel(a) > decisionLevel(b);
  });

  ClauseAdded.trigger(conflict);

  assert(decisionLevel(conflict[0]) == env.level());

  int max_level{
      (conflict.size() > 1 ? decisionLevel(conflict[1]) : init_level)};

  //  assert(upper(HORIZON) == ub - Gap<T>::epsilon());

#ifdef DBG_CL
//  if (++num_clauses > DBG_CL)
//    exit(1);
  if (cl_file != NULL) {
    *cl_file << "0 " << (conflict.size() + 1) << " 0 1 " << upper(HORIZON);
  }
#endif

  for (size_t i{0}; i < conflict.size(); ++i) {

#ifdef DBG_CL
    if (cl_file != NULL) {
      writeLiteral(conflict[i]);
    }
#endif

    conflict[i] = clauses.newNegLiteral(conflict[i]);
  }

#ifdef DBG_CL
  if (cl_file != NULL)
    *cl_file << std::endl;
#endif

  restoreState(max_level);

#ifdef DBG_TRACE
  auto cl =
#endif
      clauses.add(conflict.begin(), conflict.end(), true);

#ifdef DBG_TRACE
  if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
    std::cout << "learn conflict";
    if (clauses.size() > 0 and cl != NULL)
      clauses.displayClause(std::cout, cl);
    //      else
    //          std::cout << "(" << << ")";
    std::cout << std::endl;
  }
#endif

#ifdef DBG_CL
  if (++num_clauses > DBG_CL)
    exit(1);
#endif
}



template<typename T>
void Scheduler<T>::backtrack(Explanation e) {

  ++num_fails;

#ifdef DBG_TRACE
  if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
    conflict_set.clear();
    std::cout << "failure @level " << env.level() << "/" << init_level
              << " b/c " << e << ":\n";
    e.explain(NoLit, conflict_set);
    for (auto gl : conflict_set) {
      std::cout << " " << prettyLiteral(gl);
    }
    std::cout << std::endl;
  }
#endif

  ConflictEncountered.trigger(e);
  propagation_queue.clear();

  ++num_backtracks;

  if (env.level() == init_level) {
    throw SearchExhausted();
  }

  try {
    if (options.learning)
      learnConflict(e);
    else
      branchRight();
  } catch (const Failure &f) {
    backtrack(f.reason);
  }
}

template<typename T>
void Scheduler<T>::restart(const bool on_solution) {
  env.restore(init_level);
  undo();

  if (on_solution) {
    restart_policy->initialize(restart_limit);
    restart_limit += num_fails;
    //      std::cout << num_fails << " / " << restart_limit << std::endl;
  } else {
    restart_policy->reset(restart_limit);
    //      std::cout << num_fails << " / " << restart_limit << std::endl;
    //    displayStats(std::cout, "             ");
  }

  SearchRestarted.trigger();
}

#ifdef DBG_TRACE
template<typename T>
void Scheduler<T>::printTrace() const {
  if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
    //      std::cout << "bounds = ";
    //
    //      std::cout << "\n";
    display(std::cout, (DBG_TRACE & DOMAINS), (DBG_TRACE & BRANCH),
            false, (DBG_TRACE & CLAUSES), false, false);
  }
}
#endif

template <typename T> bool Scheduler<T>::hasSolution() const {
  return not best_solution.empty();
}

template <typename T> void Scheduler<T>::initialize_baseline() {
  baseline = new Scheduler<T>(options);
  baseline->resize(numEvent());
  for (size_t i{0}; i < numTask(); ++i)
//    baseline->newTask(minDuration(i), maxDuration(i));
      baseline->newTask(tasks[i].minDuration(), tasks[i].maxDuration());
  for (var x{0}; x < static_cast<var>(numVariable()); ++x) {
    baseline->newVariable(getEdge(NEG(x)), getEdge(POS(x)));
  }
  for (lit l{0}; l < static_cast<lit>(numEdgeLiteral()); ++l) {
    assert(getEdge(l) == baseline->getEdge(l));
  }
  auto n{static_cast<event>(numEvent())};
  auto &arcs{domain.getForwardGraph()};
  for (event i{0}; i < n; ++i) {
    for (auto j : arcs[i]) {
      baseline->newMaximumLag(i, j, j.label());
    }
  }
  baseline->setprimalBound(upper(HORIZON));

  //    std::cout << *this << std::endl;
  //  std::cout << *baseline << std::endl;
}

//<<<<<<< HEAD
// template <typename T> lit Scheduler<T>::polarityChoice(const var cp) {
//  lit d{NoLit};
//  if ((random() % 10) == 0) {
//    d = (random() % 2 ? POS(cp) : NEG(cp));
//  } else {
//    auto prec_a{getEdge(POS(cp))};
//    auto prec_b{getEdge(NEG(cp))};
//
//    T gap_a{0};
//    if (prec_a != Constant::NoEdge<T>)
//      gap_a = upper(prec_a.from) - lower(prec_a.to);
//
//    T gap_b{0};
//    if (prec_b != Constant::NoEdge<T>)
//      gap_b = upper(prec_b.from) - lower(prec_b.to);
//
//    //      double g{1};
//    ////          std::cout << static_cast<double>(gap_a) << " <> " <<
//    /// static_cast<double>(gap_b) << ": " << gap_ratio << " -> ";
//    //          if(gap_a < gap_b) {
//    ////              std::cout << "(" <<
//    ///(static_cast<double>(gap_a)/static_cast<double>(gap_b)) << ") ";
//    //
//    //              g =
//    (static_cast<double>(gap_a)/static_cast<double>(gap_b));
//    //          } else {
//    ////              std::cout << "(" <<
//    ///(static_cast<double>(gap_b)/static_cast<double>(gap_a)) << ") ";
//    //
//    //              g = (gap_a > 0 ?
//    //              static_cast<double>(gap_b)/static_cast<double>(gap_a) :
//    1);
//    //          }
//    ////          std::cout << gap_ratio << "\n";
//    //
//    //      gap_ratio += g;
//    //      if(g == 1) {
//    //          ++num_tight;
//    //      }
//
//    d = (gap_a < gap_b ? NEG(cp) : POS(cp));
//  }
//
//#ifdef DBG_TRACE
//  if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
//    std::cout << "\n-- new decision: " << prettyLiteral(EDGE(d))
//              << std::endl;
//  }
//#endif
//
//  return d;
//}
//
//=======
//>>>>>>> 953a3b06c78681b75576b6e0f8ddc428e1cbc5e3
template <typename T> void Scheduler<T>::initializeSearch() {
  if (options.minimization >= 1000)
    initialize_baseline();

  heuristic.emplace(*this, options);
  valueHeuristic.emplace(*this);

  restart_policy->initialize(restart_limit);
  start_time = cpu_time();
  //  init_level = env.level();
}

template <typename T>
template <typename S>
void Scheduler<T>::optimize_dichotomy(S &objective) {
  initializeSearch();
  objective.initDual();
  bool sat;

  displayHeader(std::cout, 62);

  while (objective.gap()) {

    saveState();
    std::cout << std::setw(6) << std::right << objective.dualBound() << " .. "
              << std::setw(6) << std::left << objective.primalBound()
              << std::right;
    displayProgress(std::cout);

    auto target = (objective.primalBound() + objective.dualBound()) / 2;

    try {
      objective.apply(target);
      sat = search();
    } catch (Failure &f) {
      sat = false;
    }

    if (KillHandler::instance().signalReceived())
      break;

    if (sat) {
      auto best{objective.value()};
      restart(true);
      try {
        objective.setPrimal(best);
      } catch (Failure &f) {
        objective.setDual(objective.primalBound());
      }
    } else {
      clearLearnedClauses();
      restart(true);
      objective.setDual(target + Gap<T>::epsilon());
    }

    assert(env.level() == 1);
    restoreState(0);
  }
}

template <typename T>
template <typename S>
void Scheduler<T>::optimize(S &objective) {

  initializeSearch();
  displayHeader(std::cout);
  std::cout << std::right;

  while (objective.gap() and not KillHandler::instance().signalReceived()) {
    //      std::cout << objective.gap() << std::endl;
    auto satisfiability = search();
    if (satisfiability == True) {
      auto best{objective.value()};
      if (options.verbosity >= Options::NORMAL) {
        //        objective.display(std::cout);
        //        displayStats(std::cout, "");
        std::cout << std::setw(10) << best;
        displayProgress(std::cout);
      }
      best_solution = polarity;
      restart(true);
      try {
        objective.setPrimal(best);
      } catch (Failure &f) {
        objective.setDual(objective.primalBound());
      }
    } else if (satisfiability == False) {
      objective.setDual(objective.primalBound());
    }
    //    else {
    //        std::cout << "here\n";
    //        std::cout << objective.gap() << std::endl;
    //    }
  }

  displaySummary(std::cout, (objective.gap() > 0 ? "killed" : "optimal"));
}

template <typename T> boolean_state Scheduler<T>::satisfiable() {
  initializeSearch();
  return search();
}

template <typename T>
template <typename S>
void Scheduler<T>::updatedualBound(S &objective) {
  if (env.level() == init_level) {
    auto lb{objective.value()};
    if (lb > objective.dualBound()) {
      objective.setDual(lb);
      //      objective.display(std::cout);
      //      displayStats(std::cout, "");
    }
  }
}

template <typename T> boolean_state Scheduler<T>::search() {

  init_level = env.level();
  boolean_state satisfiability{Unknown};
  while (satisfiability == Unknown and
         not KillHandler::instance().signalReceived()) {

    ++num_choicepoints;
    try {
      propagate();
      //      updatedualBound(objective);

      // make a checkpoint
      saveState();

      // all resource constraints are accounted for => a solution has been found
      if (search_vars.empty()) {
        satisfiability = True;
      } else {
        ++num_choicepoints;

#ifdef DBG_TRACE
        if (DBG_BOUND) {
          std::cout << "--- search node (lvl=" << env.level()
                    << ") [i=" << num_choicepoints << "] ---\n";
          printTrace();
        }
#endif

        var x = heuristic->nextChoicePoint(*this);
        lit d = valueHeuristic->choosePolarity(x, *this);
        set(d);
      }
    } catch (const Failure &f) {
      try {
        backtrack(f.reason);
        if (num_fails > restart_limit) {
          restart();
        }
      } catch (const SearchExhausted &f) {
        satisfiability = False;
        //        break;
      }
    }
  }

  //    std::cout << "return " << satisfiability << " / " << Unknown <<
  //    std::endl;

  return satisfiability;
}

// template <typename T> void Scheduler<T>::search() {
//
//   assert(not hasSolution());
//
//     if(options.minimization >= 1000)
//         initialize_baseline();
//
//   heuristic.emplace(*this, options);
//
//   restart_policy->initialize(restart_limit);
//   start_time = cpu_time();
//   //      ground_facts = numLiterals();
//   //  bool SAT{false};
//
//   // initialisation
//   lb = lower(HORIZON);
//   ub = upper(HORIZON) + 1;
//
//   init_level = env.level();
//
//   //    size_t ground_arcs{domain.arcCount()};
//
//   while (lb < ub and not KillHandler::instance().signalReceived()) {
//
//     //        if(domain.arcCount() != (ground_arcs + numVariable() -
//     //        search_vars.size())) {
//     //            std::cout << domain.arcCount() << " / " << (ground_arcs +
//     //            numVariable() - search_vars.size()) << std::endl; exit(1);
//     //        }
//
//     ++num_choicepoints;
//
//#ifdef DBG_TRACE
//     if (DBG_BOUND) {
//       std::cout << "--- search node (lvl=" << env.level()
//                 << ") [i=" << num_choicepoints << "] ---\n";
//       printTrace();
//     }
//#endif
//
//     try {
//
//#ifdef DBG_SOL
//         var x{is_on_track()};
//         try {
//#endif
//
//             propagate();
//
//#ifdef DBG_SOL
//         } catch(Failure &f) {
//             if(x == NoVar) {
//                 std::cout << "failure while on track at (" <<
//                 num_choicepoints << ")\n"; exit(1);
//             }
//             throw f;
//         }
//         if(x == NoVar and (is_on_track() != NoVar)) {
//             std::cout << "get out of track after propag at (" <<
//             num_choicepoints << ")\n"; std::cout << "literal " <<
//             edges[LIT(x,value(x))] << " whereas it is "
//             << edges[LIT(x,ref_solution[x])] << " in sol.txt\n";
//             exit(1);
//         }
//#endif
//
//
//#ifdef DBG_TRACE
//       if (DBG_BOUND and (DBG_TRACE & PROPAGATION)) {
//         std::cout << "--- propagation ---\n";
//         printTrace();
//       }
//#endif
//
//       assert(propagation_queue.empty());
//
//       if (env.level() == init_level) {
//
//         if (lower(HORIZON) > lb) {
//           lb = lower(HORIZON);
//           if (options.verbosity >= Options::NORMAL)
//             displayStats(std::cout, " new lb");
//         }
//       }
//
//       // make a checkpoint
//       saveState();
//
//       //          auto cp = search_vars.any();
//
//       // all resource constraints are accounted for => a solution has been
//       found if (search_vars.empty()) {
//
//         //        SAT = true;
//         notifySolution();
//
//       } else {
//         ++num_choicepoints;
//
//         var cp = heuristic->nextChoicePoint(*this);
//         lit d;
////#ifdef DBG_SOL
////          if(not ref_solution.empty())
////              d = LIT(cp, ref_solution[cp]);
////          else if ((random() % 10) == 0) {
////              d = (random() % 2 ? POS(cp) : NEG(cp));
////            } else {
////              auto prec_a{getEdge(POS(cp))};
////              auto prec_b{getEdge(NEG(cp))};
////              auto gap_a = upper(prec_a.from) - lower(prec_a.to);
////              auto gap_b = upper(prec_b.from) - lower(prec_b.to);
////              d = (gap_a < gap_b ? NEG(cp) : POS(cp));
////            }
////
////#else
//        if ((random() % 10) == 0) {
//          d = (random() % 2 ? POS(cp) : NEG(cp));
//        } else {
//          auto prec_a{getEdge(POS(cp))};
//          auto prec_b{getEdge(NEG(cp))};
//          auto gap_a = upper(prec_a.from) - lower(prec_a.to);
//          auto gap_b = upper(prec_b.from) - lower(prec_b.to);
//          d = (gap_a < gap_b ? NEG(cp) : POS(cp));
//        }
////#endif
//
//#ifdef DBG_TRACE
//        if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
//          std::cout << *this << "\n-- new decision: " << edges[d] <<
//          std::endl;
//        }
//#endif
//
//        set(d);
//      }
//    } catch (const Failure &f) {
//      try {
//        backtrack(f.reason);
//      } catch (const SearchExhausted &f) {
//
//#ifdef DBG_TRACE
//        if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
//          if (hasSolution())
//            std::cout << " => optimal! (" << ub << ")\n";
//          else
//            std::cout << " => unfeasible!\n";
//        }
//#endif
//
//        lb = ub;
//        if (options.verbosity > Options::SILENT)
//          displayStats(std::cout,
//                       (hasSolution() ? "optimal" : "unhasSolution"));
//      }
//
//      if (num_fails > restart_limit) {
//        restart();
//      }
//    }
//  }
//
//#ifdef DBG_TRACE
//  std::cout << "--- end search ---\n";
//#endif
//
//  //  return hasSolution();
//}

template <typename T> bool Scheduler<T>::stoppingCondition() const {
  return conflict_set.size() == 1 and LTYPE(conflict_set[0]) == EDGE_LIT;
}

template <typename T> bool Scheduler<T>::needExplanation(const genlit l) const {
  return
      //        LTYPE(l) == BOUND_LIT or
      decisionLevel(l) == env.level();
}

template <typename T> bool Scheduler<T>::groundFact(const genlit l) const {
  return decisionLevel(l) <= init_level;
}

template <typename T> bool Scheduler<T>::inExplanation(const genlit l) const {
  auto lt{LTYPE(l)};
  if (lt == EDGE_LIT) {
    return visited_edge[VAR(FROM_GEN(l))];
  }
  return domain.bounds.visited[FROM_GEN(l)];
}

template <typename T> void Scheduler<T>::markVisited(const genlit l) {
  auto lt{LTYPE(l)};
  if (lt == EDGE_LIT) {
    assert(visited_edge[VAR(FROM_GEN(l))] == false);
    visited_edge[VAR(FROM_GEN(l))] = true;
    conflict_edges.push_back(VAR(FROM_GEN(l)));
  } else {
    domain.bounds.visited[FROM_GEN(l)] = true;
    conflict_bounds.push_back(FROM_GEN(l));
  }
}

template <typename T> Explanation Scheduler<T>::getExplanation(const genlit l) const {
  return LTYPE(l) == BOUND_LIT ? domain.bounds.getExplanation(FROM_GEN(l))
                               : reason[VAR(FROM_GEN(l))];
}

template <typename T>
bool Scheduler<T>::isResponsibleForEdge(const lit l,
                                        const Constraint *c) const {
  return reason[VAR(l)].expl == c;
}

template <typename T>
bool Scheduler<T>::isResponsibleForBound(const lit l,
                                         const Constraint *c) const {
  return domain.bounds.getExplanation(l).expl == c;
}

template <typename T> void Scheduler<T>::analyze(Explanation e) {

  if (env.level() == init_level)
    return;

#ifdef DBG_TRACE
  if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
    std::cout << "--- [conflict analysis] ---\n";
    printTrace();
  }
#endif

  conflict_edges.clear();
  conflict_bounds.clear();

  conflict.clear();
  conflict_set.clear();

#ifdef DBG_TRACE
  if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
    for (var x{0}; x < static_cast<var>(numVariable()); ++x) {
      assert(visited_edge[x] == false);
    }
      for(size_t i{0}; i<domain.bounds.visited.size(); ++i) {
          assert(domain.bounds.visited[i] == false);
      }

    std::cout << "analyze failure: " << e << "\n";
  }
#endif

  lit l{NoLit};
  Explanation el{e};

  do {

    resolve(l, el);

#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
      std::cout << "cflt:";
      for (auto il : conflict) {
        std::cout << " " << prettyLiteral(il);
        assert(not needExplanation(il));
      }
      std::cout << " / stack:";
      for (auto il : conflict_set) {
        std::cout << " " << prettyLiteral(il);
      }
      std::cout << std::endl;
    }
#endif

    if (stoppingCondition()) {
      break;
    } else {

      assert(not conflict_set.empty());

      l = conflict_set[0];
      heap::remove_min(conflict_set.begin(), conflict_set.end(), younger);
      conflict_set.pop_back();

      el = getExplanation(l);
    }

  } while (true);

  conflict.push_back(conflict_set[0]);

  //  for (auto x : conflict_edges)
  //    visited_edge[x] = false;
  //
  //  for (auto l : conflict_bounds)
  //    domain.bounds.visited[l] = false;
  //    clearVisited();
}

template <typename T> void Scheduler<T>::clearVisited() {
  for (auto x : conflict_edges)
    visited_edge[x] = false;

  for (auto l : conflict_bounds)
    domain.bounds.visited[l] = false;

#ifdef DBG_TRACE

  for (lit e{0}; e < static_cast<lit>(numBoundLiteral()); ++e) {
    assert(domain.bounds.visited[e] == false);
  }
  for (var x{0}; x < static_cast<var>(numVariable()); ++x) {
    assert(visited_edge[x] == false);
  }

#endif
}

template <typename T>
void Scheduler<T>::resolve(const lit l, Explanation e) {

#ifdef DBG_TRACE
  if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
    std::cout << "\nresolution step: " << prettyLiteral(l) << " <- " << e
              << std::endl;
  }
#endif

  // add l's reason to the conflict set (at the back)
  auto prev_size{conflict_set.size()};
  e.explain(l, conflict_set);
  auto beg_expl{conflict_set.begin() + prev_size};

  // the literals from previous levels are directly added to the learned
  // conflict
  auto last{conflict_set.end()};
  auto first{beg_expl};

#ifdef DBG_CLPLUS

  if (num_clauses++ > DBG_CL)
    exit(1);

  if (cl_file != NULL) {

    auto count{0};
    for (auto i{first}; i != last; ++i) {
      if (not groundFact(*i)) {
        ++count;
      }
    }

    *cl_file << "1 " << (count + (l != NoLit ? 2 : 1)) << " 0 1 "
             << upper(HORIZON);

    if (l != NoLit) {
      if (LTYPE(l) == EDGE_LIT) {
        assert(GNOT(l) == EDGE(NOT(FROM_GEN(l))));
        writeLiteral(GNOT(l));
      } else {
        writeConstraint(~(domain.bounds.getConstraint(FROM_GEN(l))));
      }
    }
  }
#endif

  for (auto i{first}; i != last; ++i) {

#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
      std::cout << "> " << prettyLiteral(*i);
    }
#endif

#ifdef DBG_CLPLUS
    if (cl_file != NULL) {
      if (not groundFact(*i)) {
        writeLiteral(*i);
      }
    }
#endif

    if (inExplanation(*i) or groundFact(*i)) {
#ifdef DBG_TRACE
      if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
        if (inExplanation(*i)) {
          std::cout << ": already in cut\n";
        } else {
          std::cout << ": ground fact\n";
        }
      }
#endif
      continue;
    }

    markVisited(*i);

    if (needExplanation(*i)) {
#ifdef DBG_TRACE
      if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
        std::cout << ": need further explanation (" << decisionLevel(*i) << "/"
                  << env.level() << ")\n";
      }
#endif

      *first++ = std::move(*i);
    } else {
#ifdef DBG_TRACE
      if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
        std::cout << ": in clause!\n";
      }
#endif

      conflict.push_back(*i);
    }
  }

#ifdef DBG_CLPLUS
  if (cl_file != NULL) {
    *cl_file << std::endl;
  }
#endif

  // resize the conflict set
  while (conflict_set.end() != first) {
    conflict_set.pop_back();
  }

  // reorder the literals
  while (beg_expl != conflict_set.end()) {
    ++beg_expl;
    heap::percolate_up(conflict_set.begin(), prev_size++, younger);
  }

  assert(prev_size == conflict_set.size());
}

template<typename T>
std::string Scheduler<T>::prettyEvent(const event e) const {
    if(e == ORIGIN)
        return "src";
    if(e == HORIZON)
        return "C_max";
    return std::string(1,etype(e)) + std::to_string(tasks[event2task_map[e]].id());
}

template<typename T>
std::string Scheduler<T>::prettyLiteral(const genlit el) const {
  if (el == NoLit) {
    return "failure";
  } else if (LTYPE(el) == EDGE_LIT) {
    std::stringstream ss;
      if (edges[FROM_GEN(el)] != Constant::NoEdge<T>) {
          ss << "(" << el << ") [" ;
          domain.displayConstraint(ss, edges[FROM_GEN(el)]);
          //<< edges[FROM_GEN(el)]
          ss << "]"; // {" << VAR(FROM_GEN(el)) << "}";
      } else
      ss << (SIGN(FROM_GEN(el)) ? "+" : "-") << "x_" << VAR(FROM_GEN(el));
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
      ss << " (" << getIndex(VAR(FROM_GEN(el))) << ")";
    }
#endif
    return ss.str();
  } else {
    assert(FROM_GEN(el) < static_cast<lit>(numBoundLiteral()));
    std::stringstream ss;
      ss << "(" << el << ") [" ;
      domain.bounds.displayLiteral(ss,FROM_GEN(el));
      ss << "]";
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
      ss << " (" << stamp(FROM_GEN(el)) << ")";
    }
#endif
    return ss.str();
  }
}

template <typename T>
std::ostream &Scheduler<T>::displayVariables(std::ostream &os) const {
  os << "choicepoints:\n";
  for (auto x : search_vars) {
    os << prettyLiteral(EDGE(POS(x))) << " xor " << prettyLiteral(EDGE(NEG(x)))
       << std::endl;
  }
  return os;
}

template<typename T>
std::ostream &Scheduler<T>::displayConstraints(std::ostream &os) const {
  if (evt_constraint_network.size() >= 2 * numEvent())
    for (lit l{0}; l < endBoundLit(); ++l)
      if (evt_constraint_network.has(l)) {
        if (evt_constraint_network.outdegree(l) > 0) {
          os << "triggers on " << prettyEventLit(l) << ":";
          for (auto c : evt_constraint_network[l]) {
            os << " " << *(constraints[c]);
          }
          os << std::endl;
        }
      }
  if (var_constraint_network.size() >= 2 * numVariable())
    for (lit l{0}; l < endEdgeLit(); ++l)
      if (var_constraint_network.has(l)) {
        if (var_constraint_network.outdegree(l) > 0) {
          os << "triggers on " << (SIGN(l) ? edges[VAR(l)] : ~edges[VAR(l)])
             << ":";
          for (auto c : var_constraint_network[l]) {
            os << " " << *(constraints[c]);
          }
          os << std::endl;
        }
      }
  return os;
}

template <typename T>
std::ostream &Scheduler<T>::displayBranch(std::ostream &os) const {
  os << "branch:";
  for (auto pi{search_vars.fbegin()}; pi != search_vars.fend(); ++pi) {
    //      if(edges[LIT(*pi, polarity[*pi])] != Constant::NoEdge<T>)
    //          os << " [" << edges[LIT(*pi, polarity[*pi])] << "]";
    //      else
    os << prettyLiteral(EDGE(LIT(*pi, polarity[*pi])));
  }
  os << std::endl;
  return os;
}

template <typename T>
std::ostream &Scheduler<T>::displayHeader(std::ostream &os,
                                          const int width) const {
  os << std::right << std::setw(width)
     << " objective | failures | branches | clauses |  size | cpu" << std::endl
     << std::left;
  os << std::setfill('=') << std::setw(width) << "=" << std::setfill(' ')
     << std::endl;
  return os;
}

template <typename T>
std::ostream &Scheduler<T>::displaySummary(std::ostream &os,
                                           std::string msg) const {

  os << std::setfill('=') << std::setw(59) << "\n" << std::setfill(' ');
  //    auto offset{msg.size()/2};
  //    os << std::setfill('=') << std::setw(36 + offset) << msg << std::setw(36
  //    - offset)<< "\n" << std::setfill(' ');
  os << std::setw(10) << msg;
  displayProgress(os);
  os << std::setfill('=') << std::setw(59) << "\n" << std::setfill(' ');
  return os;
}

template <typename T>
std::ostream &Scheduler<T>::displayProgress(std::ostream &os) const {

  os << "  " << std::setw(9) << num_fails << "  " << std::setw(9)
     << num_choicepoints
     << "  "
     //<< std::setw(12) << num_literals << "  "
     << std::setw(8) << clauses.size();
  if (clauses.size() == 0)
    os << "    n/a ";
  else
    os << "  " << std::setw(6) << std::setprecision(4)
       << static_cast<double>(clauses.volume()) /
              static_cast<double>(clauses.size());
  //      << clauses.volume() << " " << clauses.size();

  //    os << "  " << std::setw(6) << std::setprecision(4) << (gap_ratio /
  //    static_cast<double>(num_choicepoints)) ; os << "  " << std::setw(6) <<
  //    std::setprecision(4) << (static_cast<double>(num_tight) /
  //    static_cast<double>(num_choicepoints)) ;
  os << "   " << std::left << (cpu_time() - start_time) << std::right
     << std::endl;

  return os;
}

// template <typename T>
// std::ostream &Scheduler<T>::displayStats(std::ostream &os,
//                                          const char *msg) const {
//   os << msg << " fails=" << std::setw(7) << std::left << num_fails
//      << " literals=" << std::setw(12) << std::left << num_literals;
//   if (options.learning) {
//     os << " |cflct|=";
//     if (num_fails == 0)
//       os << "n/a  ";
//     else
//       os << std::setw(5) << std::left << std::setprecision(3)
//          << static_cast<double>(clauses.volume()) /
//                 static_cast<double>(clauses.size());
//   } else {
//     os << " #prop=" << std::setw(7) << std::left << num_cons_propagations;
//   }
//   os << " cpu=" << (cpu_time() - start_time) << "\n";
//   return os;
// }

template <typename T>
std::ostream &Scheduler<T>::display(std::ostream &os, const bool dom,
                                    const bool bra, const bool sva,
                                    const bool cla, const bool egr,
                                    const bool vgr, const bool con) const {

  if (bra)
    displayBranch(os);
    if (dom) {
      os << "domains:\n";
      for (auto t : tasks) {
        os << t << std::endl;
      }
    }
  if (sva)
    displayVariables(os);
  if (cla)
    os << clauses << std::endl;
  if (egr)
    os << evt_constraint_network << std::endl;
  if (vgr)
    os << var_constraint_network << std::endl;
  if (con)
    displayConstraints(os);
  return os;
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const Scheduler<T> &x) {
  return x.display(os);
}

// template<typename T>
// std::ostream &operator<<(std::ostream &os, const DistanceConstraint<T> &x) {
//   return x.display(os);
// }
}

#endif

