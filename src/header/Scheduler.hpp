#ifndef _TEMPO_SCHEDULER_HPP
#define _TEMPO_SCHEDULER_HPP


#include <sstream>
#include <fstream>

#include "ClauseBase.hpp"
#include "ConstraintQueue.hpp"
#include "DisjunctiveEdgeFinding.hpp"
#include "EdgeConstraint.hpp"
#include "Global.hpp"
#include "Heap.hpp"
#include "Options.hpp"
#include "Restart.hpp"
#include "TemporalNetwork.hpp"
#include "heuristics/HeuristicManager.hpp"
#include "util/KillHandler.hpp"
#include "util/SubscribableEvent.hpp"

namespace tempo {


// literals x - y <= k (with x,y pointing to vars and k a constant)
template <typename T> class DistanceConstraint {

public:
  DistanceConstraint(const event f, const event t, const T d)
      : from(f), to(t), distance(d) {}

  event from;
  event to;

  T distance;

  DistanceConstraint<T> operator~() const;

  static const DistanceConstraint<T> none;

  bool entails(const DistanceConstraint<T> &e) const;
  bool contradicts(const DistanceConstraint<T> &e) const;

  std::ostream &display(std::ostream &os) const;
};



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
  task newTask(const T min_dur, const T max_dur);
  var newVariable(const DistanceConstraint<T> &if_true,
                  const DistanceConstraint<T> &if_false);

  void newPrecedence(event x, event y, T delay);
  void newMaximumLag(event x, event y, T maxlag);

  template <typename ItTask, typename ItVar>
  void postEdgeFinding(const ItTask beg_task, const ItTask end_task,
                       const ItVar beg_var, const ItVar end_var);

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
  T minDuration(const task) const;
  T maxDuration(const task) const;
  T getMakespan() const;

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

  void setUpperBound(const T b);

  void trigger(const lit l);
  void propagate();

  void saveState();
  void restoreState(const int);
  void undo() override;
  void search();
  void restart(const bool on_solution = false);
  void notifySolution();
  void backtrack(Explanation e);
  void branchRight();
  void learnConflict(Explanation e);
  bool satisfiable() const;

  int decisionLevel(const genlit l) const;
  //    void incrementActivity(const var x);

  void getCriticalPath(std::vector<genlit> &clause);

  // explanation stuff
  void xplain(const lit l, const hint h, std::vector<lit> &Cl) override;
  std::ostream &print_reason(std::ostream &, const hint) const override;
  int getType() const override;

  std::ostream &display(std::ostream &os, const bool dom = true,
                        const bool bra = true, const bool cla = false,
                        const bool egr = false, const bool vgr = false,
                        const bool con = false) const;
  std::ostream &displayBranch(std::ostream &os) const;
  std::ostream &displayStats(std::ostream &os, const char *msg) const;
  std::ostream &displayConstraints(std::ostream &os) const;

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

private:
  Options options;

  BacktrackEnvironment env;

  /*
  // the current temporal graph (manages its own reversibility), used to compute
  the bounds, i.e., the shortest paths from an to the origin
  DirectedGraph<StampedLabeledEdge<T>> core;

  // the following two vectors are indexed with LOWERBOUND(e)/UPPERBOUND(e)
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

  TemporalNetwork<T> domain;

  ClauseBase<T> clauses;

  std::vector<DistanceConstraint<T>> edges;

  std::map<std::pair<event, event>, lit> edge_map;

  std::vector<T> min_duration;
  std::vector<T> max_duration;

  std::vector<Constraint *> constraints;

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

  void resize(const size_t n);

  //

  //    T lowerBound(const event x) const;
  //    T upperBound(const event x) const;

  std::optional<heuristics::HeuristicManager<T>> heuristic;
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

  // statistics
  T lb{0};
  T ub{INFTY};

  double start_time;

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
  void analyze(Explanation e);
  void resolve(const lit l, Explanation e);

#ifdef DBG_CL

  //  TemporalNetwork<T> cut;

  std::ofstream *cl_file{NULL};
  int num_clauses{0};

  void writeLiteral(const lit l) const;
  void writeConstraint(const BoundConstraint<T> &c) const;
#endif
};

/// DISTANCE

template <typename T>
bool operator==(const DistanceConstraint<T> &d1,
                const DistanceConstraint<T> &d2) {
  return d1.from == d2.from and d1.to == d2.to and d1.distance == d2.distance;
}

template <typename T>
const DistanceConstraint<T>
    DistanceConstraint<T>::none = DistanceConstraint<T>(-1, -1, -1);

template <typename T>
DistanceConstraint<T> DistanceConstraint<T>::operator~() const {
  return {to, from, -distance - Gap<T>::epsilon()};
}

template <typename T>
bool DistanceConstraint<T>::entails(const DistanceConstraint<T>& e) const {
  return e.from == from and e.to == to and distance <= e.distance;
}

template <typename T>
bool DistanceConstraint<T>::contradicts(const DistanceConstraint<T> &e) const {
  return e.from == to and e.to == from and e.distance + distance < 0;
}

template <typename T>
std::ostream &DistanceConstraint<T>::display(std::ostream &os) const {
  os << etype(to) << TASK(to) << " - " << etype(from) << TASK(from)
     << " <= " << distance;
  return os;
}



/// SCHEDULER
///
template <typename T>
Scheduler<T>::Scheduler(Options opt)
    : ReversibleObject(&env), options(std::move(opt)), domain(*this),
      clauses(*this), evt_constraint_network(&env),
      var_constraint_network(&env), propagation_queue(constraints),
      edge_propag_pointer(0, &env), bound_propag_pointer(0, &env) {
  resize(2);

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

template<typename T>
size_t Scheduler<T>::numTask() const {
  return min_duration.size();
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

template<typename T>
void Scheduler<T>::resize(const size_t n) {

//    assert(env.level() == 0);

domain.resize(n);
//    front_propag_pointer = 0;
edge_propag_pointer = 0;

//    visited_bound[LOWER].resize(n, INFTY);
//    visited_bound[UPPER].resize(n, INFTY);
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

template <typename T>
task Scheduler<T>::newTask(const T min_dur, const T max_dur) {
  assert(env.level() == 0);
  assert(numEvent() >= 2);

  task ti{static_cast<task>(numTask())};
  resize(END(ti) + 1);

  newMaximumLag(START(ti), END(ti), max_dur);
  newMaximumLag(END(ti), START(ti), -min_dur);

  min_duration.push_back(min_dur);
  max_duration.push_back(max_dur);

  return ti;
}

template <typename T>
var Scheduler<T>::newVariable(const DistanceConstraint<T> &if_true,
                              const DistanceConstraint<T> &if_false) {
  assert(env.level() == 0);

  var x{static_cast<var>(numVariable())};

  std::pair<event, event> ptrue{if_true.from, if_true.to};
  edge_map[ptrue] = static_cast<lit>(edges.size());
  edges.emplace_back(if_true);

  std::pair<event, event> pfalse{if_false.from, if_false.to};
  edge_map[pfalse] = static_cast<lit>(edges.size());
  edges.emplace_back(if_false);

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

  post(new EdgeConstraint<T>(*this, POS(x)));
  post(new EdgeConstraint<T>(*this, NEG(x)));

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

template <typename T> T Scheduler<T>::minDuration(const task t) const {
  return min_duration[t];
}

template <typename T> T Scheduler<T>::maxDuration(const task t) const {
  return max_duration[t];
}

template <typename T> T Scheduler<T>::getMakespan() const { return ub; }

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
template <typename ItTask, typename ItVar>
void Scheduler<T>::postEdgeFinding(const ItTask beg_task, const ItTask end_task,
                                   const ItVar beg_var, const ItVar end_var) {
  post(new DisjunctiveEdgeFinding(*this, beg_task, end_task, beg_var, end_var));
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

template<typename T>
void Scheduler<T>::setUpperBound(const T b) {
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
    std::cout << "new edge literal: " << edges[l];
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
      std::cout << " b/c " << e;
    }
    std::cout << std::endl;
  }
#endif

  var x{VAR(l)};
  auto val{SIGN(l)};
  if (isUndefined(x)) {
    search_vars.remove_front(x);
    polarity[x] = val;

    //        if (val) {
    //            branch.remove_front(x);
    //        } else {
    //            branch.remove_back(x);
    //        }
    //        propagation_stack.push_back(l);

    reason[x] = e;
    var_level[x] = env.level();
#ifdef RECOVER
    domain.newEdge(edges[l].from, edges[l].to, edges[l].distance, {this, l});
#else
    domain.newEdge(edges[l].from, edges[l].to, edges[l].distance);
#endif

//        std::cout << "ok?\n";

  } else if (value(x) != val) {

    //        assert(false);

#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & PROPAGATION)) {
      std::cout << "FAIL on setting edge!\n";
    }
#endif

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

#ifdef DBG_TRACE
  if (DBG_TRACE & QUEUE) {
    std::cout << "triggers for (" << l << ") " << edges[l] << std::endl;
  }
#endif

  clauses.unit_propagate(EDGE(NOT(l)));

  const std::vector<int> &cons = var_constraint_network[l];
  const std::vector<unsigned> &rank = var_constraint_network.rank(l);

  for (auto i{cons.size()}; i-- > 0;) {
    //    for (unsigned i{0}; i < cons.size(); ++i) {

#ifdef DBG_TRACE
    if (DBG_TRACE & QUEUE) {
      std::cout << " -" << *(constraints[cons[i]]) << std::endl;
    }
#endif

    propagation_queue.edge_triggers(l, rank[i], cons[i]);
  }
}

template<typename T>
void tempo::Scheduler<T>::propagate() {

  while (not propagation_queue.empty() or
         (search_vars.frontsize() > edge_propag_pointer) or
         (numBoundLiteral() > static_cast<size_t>(bound_propag_pointer))) {

    while (search_vars.frontsize() > edge_propag_pointer) {
      ++num_literals;
      lit l{LIT(search_vars[edge_propag_pointer],
                polarity[search_vars[edge_propag_pointer]])};
      //
      //            std::cout << " prop: (" <<
      //            LIT(search_vars[edge_propag_pointer],
      //            polarity[edge_propag_pointer]) << ") [" <<
      //            edges[LIT(search_vars[edge_propag_pointer],
      //            polarity[search_vars[edge_propag_pointer]])] << "] <" <<
      //            polarity[search_vars[edge_propag_pointer]] << ">\n";

      trigger(l);
      ++edge_propag_pointer;
    }

    while (numBoundLiteral() > bound_propag_pointer) {
      ++num_literals;
      lit l{getBoundLiteral(bound_propag_pointer)};

      const std::vector<int> &cons = evt_constraint_network[l];
      const std::vector<unsigned> &rank = evt_constraint_network.rank(l);

#ifdef DBG_TRACE
      if (DBG_TRACE & QUEUE) {
        std::cout << "triggers for " << prettyEventLit(l) << std::endl;
      }
#endif

      // important to visit in reverse order to be robust to relax
      for (auto i{cons.size()}; i-- > 0;) {

#ifdef DBG_TRACE
        if (DBG_TRACE & QUEUE) {
          std::cout << " -" << *(constraints[cons[i]]) << std::endl;
        }
#endif

        propagation_queue.bound_triggers(l, rank[i], cons[i]);
      }

      ++bound_propag_pointer;
    }

    if (not propagation_queue.empty()) {
      auto cons{propagation_queue.pop_front()};

#ifdef DBG_TRACE
      if (DBG_TRACE & QUEUE) {
        std::cout << "propagate " << cons << std::endl;
      }
#endif

      cons->propagate();
    }
  }
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
      //            q{BOUND(domain.bounds.getIndex(LOWERBOUND(edges[h].to)))};

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
      //            q{BOUND(domain.bounds.getIndex(UPPERBOUND(edges[h].from)))};

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

//    std::cout << "RESTORE FROM " << env.level() << " TO " << l << std::endl;

env.restore(l);
}

template<typename T>
void Scheduler<T>::undo() {
  search_vars.setStart(edge_propag_pointer);
}

template<typename T>
void Scheduler<T>::notifySolution() {

  assert(lower(HORIZON) < ub);

  ub = lower(HORIZON);

  if (options.verbosity >= Options::NORMAL)
    displayStats(std::cout, " new ub");

  restart(true);

  setUpperBound(ub - Gap<T>::epsilon());

  best_solution = polarity;
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

template<typename T>
void Scheduler<T>::learnConflict(Explanation e) {

  analyze(e);

  std::sort(conflict.begin(), conflict.end(), [&](const lit a, const lit b) {
    return decisionLevel(a) > decisionLevel(b);
  });

  ClauseAdded.trigger(conflict);

  assert(decisionLevel(conflict[0]) == env.level());

  int max_level{
      (conflict.size() > 1 ? decisionLevel(conflict[1]) : init_level)};

  assert(upper(HORIZON) == ub - Gap<T>::epsilon());

#ifdef DBG_CL
  if (num_clauses++ > DBG_CL)
    exit(1);
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
      clauses.learn(conflict.begin(), conflict.end());

#ifdef DBG_TRACE
  if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
    std::cout << "learn conflict";
    if (clauses.size() > 0)
      clauses.displayClause(std::cout, cl);
    std::cout << std::endl;
  }
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
  } else {
    restart_policy->reset(restart_limit);
    if (options.verbosity >= Options::YACKING)
      displayStats(std::cout, "restart");
  }

  SearchRestarted.trigger();
}

#ifdef DBG_TRACE
template<typename T>
void Scheduler<T>::printTrace() const {
  if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
    std::cout << "bounds = [" << lb << ".." << lower(HORIZON) << ".." << ub
              << "]\n";
    display(std::cout, (DBG_TRACE & DOMAINS), (DBG_TRACE & BRANCH),
            (DBG_TRACE & CLAUSES), false, false, false);
  }
}
#endif

template<typename T>
bool Scheduler<T>::satisfiable() const {
  return not best_solution.empty();
}

template <typename T> void Scheduler<T>::search() {

  assert(not satisfiable());

  heuristic.emplace(*this, options);

  restart_policy->initialize(restart_limit);
  start_time = cpu_time();
  //      ground_facts = numLiterals();
  //  bool SAT{false};

  // initialisation
  lb = lower(HORIZON);
  ub = upper(HORIZON) + 1;

  init_level = env.level();

  //    size_t ground_arcs{domain.arcCount()};

  while (lb < ub and not KillHandler::instance().signalReceived()) {

    //        if(domain.arcCount() != (ground_arcs + numVariable() -
    //        search_vars.size())) {
    //            std::cout << domain.arcCount() << " / " << (ground_arcs +
    //            numVariable() - search_vars.size()) << std::endl; exit(1);
    //        }

    ++num_choicepoints;

#ifdef DBG_TRACE
    if (DBG_BOUND) {
      std::cout << "--- search node (lvl=" << env.level()
                << ") [i=" << num_choicepoints << "] ---\n";
      printTrace();
    }
#endif

    try {

      propagate();

#ifdef DBG_TRACE
      if (DBG_BOUND and (DBG_TRACE & PROPAGATION)) {
        std::cout << "--- propagation ---\n";
        printTrace();
      }
#endif

      assert(propagation_queue.empty());

      if (env.level() == init_level) {

        if (lower(HORIZON) > lb) {
          lb = lower(HORIZON);
          if (options.verbosity >= Options::NORMAL)
            displayStats(std::cout, " new lb");
        }
      }

      // make a checkpoint
      saveState();

      //          auto cp = search_vars.any();

      // all resource constraints are accounted for => a solution has been found
      if (search_vars.empty()) {

        //        SAT = true;
        notifySolution();

      } else {
        ++num_choicepoints;

        var cp = heuristic->nextChoicePoint(*this);
        lit d;
        if ((random() % 10) == 0) {
          d = (random() % 2 ? POS(cp) : NEG(cp));
        } else {
          auto prec_a{getEdge(POS(cp))};
          auto prec_b{getEdge(NEG(cp))};
          auto gap_a = upper(prec_a.from) - lower(prec_a.to);
          auto gap_b = upper(prec_b.from) - lower(prec_b.to);
          d = (gap_a < gap_b ? NEG(cp) : POS(cp));
        }

#ifdef DBG_TRACE
        if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
          std::cout << "-- new decision: " << edges[d] << std::endl;
        }
#endif

        set(d);
      }
    } catch (const Failure &f) {
      try {
        backtrack(f.reason);
      } catch (const SearchExhausted &f) {

#ifdef DBG_TRACE
        if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
          if (satisfiable())
            std::cout << " => optimal! (" << ub << ")\n";
          else
            std::cout << " => unfeasible!\n";
        }
#endif

        lb = ub;
        if (options.verbosity > Options::SILENT)
          displayStats(std::cout,
                       (satisfiable() ? "optimal" : "unsatisfiable"));
      }

      if (num_fails > restart_limit) {
        restart();
      }
    }
  }

#ifdef DBG_TRACE
  std::cout << "--- end search ---\n";
#endif

  //  return satisfiable();
}

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
      l = conflict_set[0];
      heap::remove_min(conflict_set.begin(), conflict_set.end(), younger);
      conflict_set.pop_back();

      el = getExplanation(l);
    }

  } while (true);

  conflict.push_back(conflict_set[0]);

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
std::string Scheduler<T>::prettyLiteral(const genlit el) const {
  if (el == NoLit) {
    return "failure";
  } else if (LTYPE(el) == EDGE_LIT) {
    std::stringstream ss;
    ss << edges[FROM_GEN(el)];
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
      ss << " (" << getIndex(VAR(FROM_GEN(el))) << ")";
    }
#endif
    return ss.str();
  } else {
    assert(FROM_GEN(el) < static_cast<lit>(numBoundLiteral()));
    std::stringstream ss;
    ss << domain.bounds.getConstraint(FROM_GEN(el));
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
      ss << " (" << stamp(FROM_GEN(el)) << ")";
    }
#endif
    return ss.str();
  }
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
    os << " [" << edges[LIT(*pi, polarity[*pi])] << "]";
  }
  os << std::endl;
  return os;
}

template<typename T>
std::ostream &Scheduler<T>::displayStats(std::ostream &os, const char* msg) const {
  os << msg << ": [" << std::left << std::setw(5) << std::setfill('.') << lb
     << std::setfill(' ');
  if (ub < INFTY)
    os << std::right << std::setw(6) << std::setfill('.') << ub
       << std::setfill(' ');
  else
    os << ".infty";
#ifndef DEBUGCMP
  os << "] fails=" << std::setw(7) << std::left << num_fails
     << " literals=" << std::setw(12) << std::left << num_literals
     << " |cflct|=";
  if (num_fails == 0)
    os << "n/a  ";
  else
    os << std::setw(5) << std::left << std::setprecision(3)
       << static_cast<double>(clauses.volume()) /
              static_cast<double>(clauses.size());
  os << " cpu=" << (cpu_time() - start_time) << "\n";
#else
  os << std::endl;
#endif
  return os;
}


template<typename T>
std::ostream &Scheduler<T>::display(std::ostream &os, const bool dom, const bool bra, const bool cla, const bool egr, const bool vgr, const bool con) const {

  if (bra)
    displayBranch(os);
  if (dom)
    os << domain;
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

template<typename T>
std::ostream &operator<<(std::ostream &os, const DistanceConstraint<T> &x) {
  return x.display(os);
}

}

#endif

