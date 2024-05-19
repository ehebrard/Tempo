
#ifndef _TEMPO_SOLVER_HPP
#define _TEMPO_SOLVER_HPP

#include "ClauseBase.hpp"
#include "Constant.hpp"
#include "ConstraintQueue.hpp"
#include "DirectedGraph.hpp"
#include "DistanceConstraint.hpp"
#include "Failure.hpp"
#include "Global.hpp"
#include "Literal.hpp"
#include "Model.hpp"
//#include "Objective.hpp"
//#include "Restart.hpp"
//#include "TemporalNetwork.hpp"
//#include "constraints/DisjunctiveEdgeFinding.hpp"
#include "constraints/EdgeConstraint.hpp"
//#include "constraints/Transitivity.hpp"
//#include "heuristics/HeuristicManager.hpp"
//#include "heuristics/impl/DecayingEventActivityMap.hpp"
//#include "util/Heap.hpp"
//#include "util/KillHandler.hpp"
#include "util/Options.hpp"
#include "util/SubscribableEvent.hpp"


namespace tempo {

// using index_t = uint32_t;

template<typename T> class Solver;


template<typename T>
class BooleanStore {

public:
  BooleanStore(Solver<T> &s);
  ~BooleanStore() = default;

  bool value(const var_t x) const;
  bool isTrue(const var_t x) const;
  bool isFalse(const var_t x) const;
  bool isUndefined(const var_t x) const;
  bool falsified(const Literal<T> l) const;
  bool satisfied(const Literal<T> l) const;

  // declare a new Boolean variable
  BooleanVar<T> newVar(const info_t s = Constant::NoSemantic);
  // declare a new Boolean variable with a semantic (disjunction)
  DisjunctVar<T> newDisjunct(const DistanceConstraint<T> &d1,
                             const DistanceConstraint<T> &d2);

  size_t size() const;

  void set(Literal<T> l);
  void undo(Literal<T> l);

  Literal<T> getLiteral(const bool s, const var_t x) const;
  index_t litIndex(const Literal<T> l) const;

  const DistanceConstraint<T> &getEdge(const bool s, const var_t x) const;
  const DistanceConstraint<T> &getEdge(const Literal<T> l) const;

  bool hasSemantic(const var_t x) const;

protected:
  Solver<T> &solver;

  std::vector<bool> polarity;

  std::vector<info_t> edge_index;

  std::vector<DistanceConstraint<T>> edges;

  std::vector<index_t> propagation_level;
};

template<typename T>
class NumericStore {

public:
  NumericStore(Solver<T> &s);
  ~NumericStore() = default;

  bool falsified(const Literal<T> l) const;
  bool satisfied(const Literal<T> l) const;

  T upper(const var_t x) const;

  T lower(const var_t x) const;

  Literal<T> strongestLiteral(const bool s, const var_t x) const;

  index_t litIndex(const Literal<T> l) const;
  index_t lastLitIndex(const bool s, const var_t x) const;

  size_t size() const;

  // declare a new numeric variable
  NumericVar<T> newVar();
  // declare a new numeric variable with temporal semantic (can be involved in
  // disjunctions and precedences)
  //  TemporalVar<T> newTemporalVar();

  void set(Literal<T> l);
  void undo(Literal<T> l);

  const std::vector<T> &get(const int b) const;

private:
  Solver<T> &solver;

  //  // the stack of Literals reprensenting all the changes so far
  //  std::vector<Literal<T>> trail;
  //  // the reason for each propagation event
  //  std::vector<Explanation> reason;

  // [for each numeric signed_var] the current bounds (repeated for efficient
  // read access)
  std::vector<T> bound[2];
  // [for each numeric signed_var] the current index in the 'propagation_events'
  // stack
  std::vector<std::vector<index_t>> bound_index[2];
  // [for each Literal] pointer to the previous Literal of the same numeric
  // variable (useful for undoing and for searching implicants)
  //  std::vector<index_t> prev_bound;
  // a 'time stamp' for each literal (the index of the most recent Boolean
  // literal)
  //  std::vector<index_t> stamp;
  // used for clause minimization
  //  std::vector<bool> visited;
};

template <typename T = int>
class Solver : public ReversibleObject, public NewExplainer<T> {

public:
  /**
   * @name constructors
   */
  //@{
  Solver(Options opt);
  ~Solver() = default;
  //@}

  /**
   * @name count accessors
   */
  //@{
  /// Total Number of variables
  //  size_t numVariable() const;
  /// Number of variable of the numeric type
  //  size_t numNumericVariable() const;
  /// Number of variables of the Boolean type
  //  size_t numBooleanVariable() const;
  /// Number of  clauses
  //  size_t numClause() const;
  /// Number of constraints
  size_t numConstraint() const;
  /// Number of  changes
  size_t numLiteral() const;
  //@}

  /**
   * @name subscribable events
   */
  ///@{
  mutable SubscribableEvent<const std::vector<Literal<T>> &>
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

  /**
   * @name modelling methods
   */
  //@{
  //    var_t newDisjunctVar(const DistanceConstraint<T> &if_true,
  //                         const DistanceConstraint<T> &if_false);
  //@}

  /**
   * @name value accessors
   */
  //@{
  //  bool value(const var_t x) const;
  //  bool isTrue(const var_t x) const;
  //  bool isFalse(const var_t x) const;
  //  bool isUndefined(const var_t x) const;
  //  bool falsified(const Literal<T> l) const;
  //  bool satisfied(const Literal<T> l) const;

  //  T upper(const var_t) const;
  //  T lower(const var_t) const;
  //@}

  /**
   * @name Literal accessors
   */
  //@{
  // get the Literal corresponding to the i-th propagation event
  Literal<T> getLiteral(const index_t i) const;

  // get the most recent Literal that entails l
  Literal<T> getImplicant(const Literal<T> l) const;

  // get the index in the propagation queue of the last Literal involving
  // variable x
  index_t propagationLevel(const Literal<T> l) const;

  void set(Literal<T> l, const NewExplanation<T> &e = Constant::NewNoReason<T>);
  //  void setNumericNoUpdate(Literal<T> l,
  //                  const NewExplanation<T> &e);
  void setNumeric(Literal<T> l,
                  const NewExplanation<T> &e = Constant::NewNoReason<T>,
                  const bool do_update = true);
  void setBoolean(Literal<T> l,
                  const NewExplanation<T> &e = Constant::NewNoReason<T>);
  void set(const DistanceConstraint<T> &c, const index_t r = Constant::NoIndex);
  void boundClosure(const var_t x, const var_t y, const T d, const index_t r);
  template <typename G>
  void update(const bool bt, const int s, const G &neighbors);
  //@}

  void post(NewConstraint<T> *);
  void relax(NewConstraint<T> *);
  void wake_me_on(const Literal<T>, const int);
  template <typename X> void addToSearch(const X &x);

  /**
   * @name search
   */
  //@{
  void trigger(const Literal<T> l);
  void propagate();
  boolean_state search();
  void backtrack(NewExplanation<T> &e);
  void branchRight();
  void learnConflict(NewExplanation<T> &e);
  void analyze(NewExplanation<T> &e);
  //    int takeDecision(const Literal<T> l);

  int saveState();
  void restoreState(const int);
  void undo() override;
  //@}

  void xplain(const Literal<T>, const hint,
              std::vector<Literal<T>> &) override; // {}
  std::ostream &print_reason(std::ostream &os,
                             const hint) const override; // { return os; }
  //    int getType() const override;// { return CYCLEEXPL; }

  BacktrackEnvironment &getEnv() { return env; }
  const Options &getOptions() const { return options; }

  /**
   * @name printing and trace
   */
  //@{
  std::ostream &display(std::ostream &os, const bool dom = true,
                        const bool bra = true, const bool sva = true,
                        const bool pre = false, const bool cla = false,
                        const bool bgr = false, const bool ngr = false,
                        const bool con = false, const bool trl = false) const;
  std::ostream &displayDomains(std::ostream &os) const;
  std::ostream &displayBranches(std::ostream &os) const;
  std::ostream &displayVariables(std::ostream &os) const;
  std::ostream &displayConstraints(std::ostream &os) const;
  //@}

  BooleanStore<T> boolean;
  BooleanVar<T> newBoolean();
  DisjunctVar<T> newDisjunct(const DistanceConstraint<T> &,
                             const DistanceConstraint<T> &);

  NumericStore<T> numeric;
  NumericVar<T> newNumeric();
  TemporalVar<T> newTemporal(const T offset = 0);
  Job<T> newJob(const T mindur = 0, const T maxdur = Constant::Infinity<T>);

  // all the clauses (learnt or from the base problem)
  NewClauseBase<T> clauses;

  // graph with all the known edges
  DirectedGraph<StampedLabeledEdge<T, index_t>> core;

  //      DifferenceLogicStore<T> precedences;

private:
  Options options;

  BacktrackEnvironment env;

  std::vector<Literal<T>> decisions;

  // the stack of Literals reprensenting all the changes so far
  std::vector<Literal<T>> trail;
  // the reason for each propagation event
  std::vector<NewExplanation<T>> reason;

  /**
   * @name domains
   */
  //@{

  // a reversible pointer to the most recent preopagation event that is not
  // yet propagated
  Reversible<index_t> propag_pointer;
  //  Reversible<index_t> var_pointer;
  //  Reversible<index_t> num_propag_pointer;
  //@}

  /**
   * @name constraints
   */
  //@{
  // data structure used to implement the overall propagation (parameter is
  // the number of priority classes)
  NewConstraintQueue<T, 3> propagation_queue;
  // all of the posted constraints
  std::vector<NewConstraint<T> *> constraints;
  // dependency graph variables/constraints
  DirectedGraph<int> boolean_constraints;
  DirectedGraph<int> numeric_constraints;
  // @}

  /**
   * @name search @TODO: search only on Boolean variables for now
   */
  //@{
  // the set of variables remaining to fix

public:
  SparseSet<var_t, Reversible<size_t>> boolean_search_vars;
  SparseSet<var_t, Reversible<size_t>> numeric_search_vars;

private:
  // current polarity of the Boolean variables
  //  std::vector<bool> polarity;
  //  // copy of the best solution so far
  //  std::vector<bool> best_solution;
  // level at which a variable has been decided
  //  std::vector<index_t> var_level;
  std::vector<bool> explored;
  //@}

  // buffers
  SparseSet<> changed;
  //    std::vector<Literal<T>> lit_buffer;
  std::vector<Literal<T>> conflict;

#ifdef DBG_TRACE
  void printTrace() const;
#endif

  //    heuristics::impl::EventActivityMap<T> *activityMap{NULL};

public:
  //    void setActivityMap(heuristics::impl::EventActivityMap<T> *map) {
  //      activityMap = map;
  //    }
  //    heuristics::impl::EventActivityMap<T> *getActivityMap() {
  //      return activityMap;
  //    }

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

  int init_level{0};

  static const Literal<T> Contradiction;

  double looseness(const Literal<T> &l) const;
};

template <typename T>
const Literal<T> Solver<T>::Contradiction = Literal<T>(false, Constant::NoVarx,
                                                       info_t(0));

#ifdef DBG_TRACE
template <typename T> void Solver<T>::printTrace() const {
  if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
    display(std::cout, (DBG_TRACE & DOMAINS), (DBG_TRACE & BRANCH), false,
            false, (DBG_TRACE & CLAUSES), false, false, false, false);
  }
}
#endif

template <typename T>
Literal<T> BooleanStore<T>::getLiteral(const bool s, const var_t x) const {
  return Literal<T>(s, x, edge_index[x]);
}

template <typename T>
index_t BooleanStore<T>::litIndex(const Literal<T> l) const {
  return propagation_level[l.variable()];
}

template <typename T>
const DistanceConstraint<T> &BooleanStore<T>::getEdge(const bool s,
                                                      const var_t x) const {
  return edges[edge_index[x] + s];
}

template <typename T>
const DistanceConstraint<T> &
BooleanStore<T>::getEdge(const Literal<T> l) const {
  return edges[l.constraint()];
}

template <typename T> bool BooleanStore<T>::hasSemantic(const var_t x) const {
  return edge_index[x] != Constant::NoSemantic;
}

template <typename T> BooleanStore<T>::BooleanStore(Solver<T> &s) : solver(s) {
  edges.push_back(Constant::NoEdge<T>);
  edges.push_back(Constant::NoEdge<T>);
}

template <typename T> size_t BooleanStore<T>::size() const {
  return polarity.size() / 2;
}

template <typename T> BooleanVar<T> BooleanStore<T>::newVar(const info_t s) {
  BooleanVar<T> x{static_cast<var_t>(size())};

  propagation_level.push_back(Constant::NoIndex);

  polarity.push_back(false);
  polarity.push_back(false);

  edge_index.push_back(s);

  return x;
}

template <typename T>
DisjunctVar<T> BooleanStore<T>::newDisjunct(const DistanceConstraint<T> &d1,
                                            const DistanceConstraint<T> &d2) {
  info_t d{static_cast<info_t>(edges.size())};
  DisjunctVar<T> x{newVar(d).id(), d};
  edges.push_back(d1);
  edges.push_back(d2);
  return x;
}

template <typename T> void BooleanStore<T>::set(Literal<T> l) {
  propagation_level[l.variable()] = solver.numLiteral();
  polarity[l] = true;
  if (l.hasSemantic()) {
    assert(l.constraint() == (edge_index[l.variable()] + l.sign()));
    solver.set(edges[l.constraint()],
               static_cast<index_t>(solver.numLiteral()));
  }
  assert(l.hasSemantic() or edge_index[l.variable()] == Constant::NoSemantic);
}

template <typename T> void BooleanStore<T>::undo(Literal<T> l) {
  polarity[l] = false;
}

template <typename T> bool BooleanStore<T>::isTrue(const var_t x) const {
  return polarity[Literal<T>::index(true, x)];
}

template <typename T> bool BooleanStore<T>::isFalse(const var_t x) const {
  return polarity[Literal<T>::index(false, x)];
}

template <typename T> bool BooleanStore<T>::isUndefined(const var_t x) const {
  return not(polarity[Literal<T>::index(true, x)] or
             polarity[Literal<T>::index(false, x)]);
}

template <typename T>
bool BooleanStore<T>::satisfied(const Literal<T> l) const {
  return polarity[l];
}

template <typename T>
bool BooleanStore<T>::falsified(const Literal<T> l) const {
  return polarity[~l];
}

// template<typename T>
// class DifferenceLogicStore {
//
// public:
//     DifferenceLogicStore(Solver<T> &s);
//   ~DifferenceLogicStore() = default;
//
// private:
//
//     // graph with all the known edges
//     DirectedGraph<StampedLabeledEdge<T, Literal<T>>> core;
//
//     std::vector<DistanceConstraint<T> edges;
//
// };

template <typename T> NumericStore<T>::NumericStore(Solver<T> &s) : solver(s) {
  //  resize(n);
}

template <typename T> size_t NumericStore<T>::size() const {
  return bound[bound::lower].size();
}

template <typename T>
const std::vector<T> &NumericStore<T>::get(const int b) const {
  return bound[b];
}

template <typename T> NumericVar<T> NumericStore<T>::newVar() {
  NumericVar<T> x{static_cast<var_t>(size())};

  bound[bound::lower].push_back(Constant::Infinity<T>);
  bound[bound::upper].push_back(Constant::Infinity<T>);

  //    bound_index[bound::lower].emplace_back(std::initializer_list<
  //    index_t>{Constant::IndexOfMin});
  //    bound_index[bound::upper].emplace_back(std::initializer_list<
  //    index_t>{Constant::IndexOfMax});

  bound_index[bound::lower].resize(size());
  bound_index[bound::upper].resize(size());
  bound_index[bound::lower].back().push_back(Constant::InfIndex);
  bound_index[bound::upper].back().push_back(Constant::InfIndex);

  return x;
}

// template <typename T> TemporalVar<T> NumericStore<T>::newTemporalVar() {
//   TemporalVar<T> x{newVar().id(), 0};
//   return x;
// }

template <typename T> void NumericStore<T>::set(Literal<T> l) {
  auto s{l.sign()};
  auto v{l.variable()};
  //  if (bound[s][v] > l.value()) {

  assert(bound[s][v] > l.value());

  bound[s][v] = l.value();
  bound_index[s][v].push_back(static_cast<index_t>(solver.numLiteral()));
  //    return true;
  //  }
  //  return false;
}

template <typename T> void NumericStore<T>::undo(Literal<T> l) {
  auto s{l.sign()};
  auto v{l.variable()};

  //    std::cout << "undo [" << l << "] " << s << " " << v << std::endl;

  //    std::cout << bound_index[s][v].size() << std::endl;

  //    for(auto i : bound_index[s][v]) {
  //        std::cout << solver.getLiteral(i) << std::endl;
  //    }

  bound_index[s][v].pop_back();
  bound[s][v] = solver.getLiteral(bound_index[s][v].back()).value();
}

// template <typename T> void NumericStore<T>::resize(const size_t n) {
//   if (n == 0)
//     return;
//
//   bound[bound::lower].resize(n, -INFTY);
//   bound[bound::upper].resize(n, INFTY);
//
//     bound_index[bound::lower].resize(n);
//     bound_index[bound::upper].resize(n);
//
// }

template <typename T> T NumericStore<T>::upper(const var_t x) const {
  return bound[bound::upper][x];
}

template <typename T> T NumericStore<T>::lower(const var_t x) const {
  return -bound[bound::lower][x];
}

template <typename T>
Literal<T> NumericStore<T>::strongestLiteral(const bool s,
                                             const var_t x) const {
  return solver.getLiteral(bound_index[s][x].back());
}

template <typename T>
index_t NumericStore<T>::lastLitIndex(const bool s, const var_t x) const {
  return bound_index[s][x].back();
}

template <typename T>
index_t NumericStore<T>::litIndex(const Literal<T> l) const {
  auto i{bound_index[l.sign()][l.variable()].rbegin()};
  while (solver.getLiteral(*i).value() < l.value())
    ++i;
  return *i;
}

template <typename T>
bool NumericStore<T>::falsified(const Literal<T> l) const {
  return -(l.value()) > bound[not l.sign()][l.variable()];
}

template <typename T>
bool NumericStore<T>::satisfied(const Literal<T> l) const {
  return l.value() >= bound[l.sign()][l.variable()];
}

template <typename T>
Solver<T>::Solver(Options opt)
    : ReversibleObject(&env), boolean(*this), numeric(*this), clauses(*this),
      core(&env), options(std::move(opt)), propag_pointer(1, &env),
      propagation_queue(constraints), boolean_constraints(&env),
      numeric_constraints(&env), boolean_search_vars(0, &env),
      numeric_search_vars(0, &env) {
  trail.emplace_back(Constant::NoVarx, Constant::Infinity<T>);
  reason.push_back(Constant::NewNoReason<T>);
  seed(options.seed);
}

template <typename T> BooleanVar<T> Solver<T>::newBoolean() {
  auto x{boolean.newVar()};
  clauses.newBooleanVar(x.id());
  //    propagation_level.resize(boolean.size());
  boolean_constraints.resize(std::max(numConstraint(), 2 * boolean.size()));
  return x;
}

template <typename T>
DisjunctVar<T> Solver<T>::newDisjunct(const DistanceConstraint<T> &d1,
                                      const DistanceConstraint<T> &d2) {
  auto x{boolean.newDisjunct(d1, d2)};
  boolean_constraints.resize(std::max(numConstraint(), 2 * boolean.size()));

  post(new NewEdgeConstraint<T>(*this, boolean.getLiteral(true, x)));
  post(new NewEdgeConstraint<T>(*this, boolean.getLiteral(false, x)));

  return x;
}

template <typename T> NumericVar<T> Solver<T>::newNumeric() {
  auto x{numeric.newVar()};
  changed.reserve(numeric.size());
  clauses.newNumericVar(x.id());
  numeric_constraints.resize(std::max(numConstraint(), 2 * numeric.size()));
  return x;
}

template <typename T> TemporalVar<T> Solver<T>::newTemporal(const T offset) {
  TemporalVar<T> x{newNumeric().id(), offset};
  core.newVertex(x);
  return x;
}

template <typename T> Job<T> Solver<T>::newJob(const T mindur, const T maxdur) {
  return Job<T>(*this, mindur, maxdur);
}

template <typename T>
size_t Solver<T>::numConstraint() const {
  return constraints.size();
}

template <typename T>
size_t Solver<T>::numLiteral() const {
  return trail.size();
}

template <typename T> Literal<T> Solver<T>::getLiteral(const index_t i) const {
  return trail[i];
}

template <typename T>
void Solver<T>::set(const DistanceConstraint<T> &c, const index_t r) {

#ifdef DBG_TRACE
  if (DBG_BOUND and (DBG_TRACE & PROPAGATION)) {
    std::cout << "set constraint " << c << std::endl;
  }
#endif

  core.emplace_edge(c.from, c.to, c.distance, r);
  boundClosure(c.from, c.to, c.distance, r);
}

template <typename T>
void Solver<T>::set(Literal<T> l, const NewExplanation<T> &e) {

#ifdef DBG_TRACE
  if (DBG_BOUND and (DBG_TRACE & PROPAGATION)) {
    std::cout << "set literal " << l << " b/c " << e << std::endl;
  }
#endif

  if (l.isNumeric()) {
    setNumeric(l, e);
  } else {
    setBoolean(l, e);
  }
}

template <typename T>
void Solver<T>::setNumeric(Literal<T> l, const NewExplanation<T> &e,
                           const bool do_update) {

  if (not numeric.satisfied(l)) {
    if (numeric.falsified(l))
      throw NewFailure<T>(e);

    numeric.set(l);
    reason.emplace_back(e);
    trail.push_back(l);

    if (do_update) {
      if (l.sign() == bound::upper) {
        update(bound::upper, l.variable(), core);
      } else {
        update(bound::lower, l.variable(), core.backward());
      }
    }
  }
}

template <typename T>
void Solver<T>::setBoolean(Literal<T> l, const NewExplanation<T> &e) {
  assert(not boolean.satisfied(l));

  if (boolean.falsified(l))
    throw NewFailure<T>(e);

  boolean.set(l);
  reason.emplace_back(e);
  trail.push_back(l);

  if (boolean_search_vars.has(l.variable()))
    boolean_search_vars.remove_back(l.variable());
}

template <typename T>
void Solver<T>::boundClosure(const var_t x, const var_t y, const T d,
                             const index_t r) {
  // closure w.r.t. 0 (0 -> x -(d)-> y -> 0)

  NewExplanation<T> e{this, static_cast<hint>(r)};
  if (r == Constant::NoIndex)
    e = Constant::NewNoReason<T>;

  if (numeric.lower(y) != -Constant::Infinity<T>) {
    setNumeric(geq<T>(x, numeric.lower(y) - d), e);
  }

  if (numeric.upper(x) != Constant::Infinity<T>) {
    setNumeric(leq<T>(y, numeric.upper(x) + d), e);
  }
}

// template<typename T>
// void Solver<T>::restart(const bool on_solution) {
//   env.restore(init_level);
//   undo();
//
//   if (on_solution) {
//     restart_policy->initialize(restart_limit);
//     restart_limit += num_fails;
//     //      std::cout << num_fails << " / " << restart_limit << std::endl;
//   } else {
//     restart_policy->reset(restart_limit);
//     //      std::cout << num_fails << " / " << restart_limit << std::endl;
//     //    displayStats(std::cout, "             ");
//   }
//
////  SearchRestarted.trigger();
//}

template <typename T> void Solver<T>::backtrack(NewExplanation<T> &e) {

  ++num_fails;

#ifdef DBG_TRACE
  if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
    //    conflict_set.clear();
    std::cout << "failure @level " << env.level() << "/" << init_level
              << " b/c " << e << ":\n";
    //    e.explain(NoLit, conflict_set);
    //    for (auto gl : conflict_set) {
    //      std::cout << " " << prettyLiteral(gl);
    //    }
    //    std::cout << std::endl;
  }
#endif

  //  ConflictEncountered.trigger(e);
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
  } catch (NewFailure<T> &f) {
    backtrack(f.reason);
  }
}

template <typename T>
index_t Solver<T>::propagationLevel(const Literal<T> l) const {
  if (l.isNumeric()) {
    return numeric.litIndex(l);
  } else {
    return boolean.litIndex(l);
  }
}

// template<typename T>
// void Solver<T>::markExplored(const Literal<T> &l) {
//     explored[propagationLevel(l)] = true;
// }

template <typename T> void Solver<T>::analyze(NewExplanation<T> &e) {
  explored.resize(trail.size(), false);
  conflict.clear();
  auto decision_lvl{propagationLevel(decisions.back())};

  int num_lit{0};
  index_t li{static_cast<index_t>(trail.size() - 1)};
  Literal<T> l{Contradiction};

#ifdef DBG_TRACE
  if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
    std::cout << "analyze conflict: ";
    displayBranches(std::cout);
  }
#endif

  NewExplanation<T> &exp = e;
  do {
    int csize{static_cast<int>(conflict.size())};
    exp.explain(l, conflict);

#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
      std::cout << "resolve ";
      if (l == Contradiction) {
        std::cout << "contradiction";
      } else {
        std::cout << l;
      }
      std::cout << " by " << exp << std::endl;
    }
#endif

    //      while(not lit_buffer.empty()) {
    //          auto p{lit_buffer.back()};
    //          auto p_lvl{propagationLevel(p)};
    //
    //          if(p_lvl < decision_lvl) {
    //              conflict.push_back(p);
    //          } else {
    //              explored[p_lvl] = true;
    //              ++num_lit;
    //          }
    //      }

    for (int i{static_cast<int>(conflict.size()) - 1}; i >= csize;) {

      auto p{conflict[i]};
      auto p_lvl{propagationLevel(p)};

#ifdef DBG_TRACE
      if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
        std::cout << " ** " << p << " (" << p_lvl << "/" << decision_lvl << ")";
      }
#endif

      if (not explored[p_lvl]) {
        if (p_lvl < decision_lvl) {

#ifdef DBG_TRACE
          if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
            std::cout << " => keep in conflict!\n";
          }
#endif

          std::swap(conflict[csize], conflict[i]);
          ++csize;
          //              conflict.push_back(p);
        } else {

#ifdef DBG_TRACE
          if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
            std::cout << " => to explore\n";
          }
#endif

          ++num_lit;
          --i;
        }
      } else
        --i;

      explored[p_lvl] = true;
    }
    conflict.resize(csize);

    while (not explored[li--])
      ;
    l = trail[li + 1];
    exp = reason[li + 1];

    //#ifdef DBG_LEARNING
    //      std::cout << " #=" << (num_lit) << std::endl ;
    //#endif

    // explored.reset(VAR(l));
    explored[li + 1] = false;
    --num_lit;

  } while (num_lit > 0);

  index_t lvl{0};

  for (auto l : conflict) {
    explored[propagationLevel(l)] = false;
  }

  //  if (not conflict.empty())
  //    lvl = var_level[VAR(
  //        *std::max_element(conflict_clause.begin(), conflict_clause.end(),
  //                          [&](const lit p, const lit q) {
  //                            return var_level[VAR(p)] < var_level[VAR(q)];
  //                          }))];

  // cout << "UIP is " << TODIMACS(l) << endl;

  //  compute_level_count();

  //  // minimization happens here
  //  if (options.minimization > 0)
  //    minimize_conflict();

  //  clear_explored();

  conflict.push_back(~l);
  std::sort(conflict.begin(), conflict.end(),
            [&](const Literal<T> a, const Literal<T> b) {
              return propagationLevel(a) > propagationLevel(b);
            });
  //  sort(conflict_clause.begin(), conflict_clause.end());

  //    for(auto p : conflict) {
  //        std::cout << p << ": " << propagationLevel(p) << std::endl;
  //    }
  //    exit(1);

  //  update_activity();

  // ++level_count[level()];

  //  return lvl;

//  for (auto l : conflict) {
//    std::cout << " " << l;
//  }
//  std::cout << std::endl;
//
//  for (size_t i{0}; i < explored.size(); ++i) {
//    assert(not explored[i]);
//  }
//
//  std::vector<bool> duplicate(2 * boolean.size(), false);
//  for (auto l : conflict) {
//    if (not l.isNumeric()) {
//      assert(not duplicate[l]);
//      duplicate[l] = true;
//    }
//  }
}

template <typename T> void Solver<T>::learnConflict(NewExplanation<T> &e) {

  analyze(e);

  //  if (options.minimization >= 0) {
  //    minimization(static_cast<size_t>(options.minimization));
  //  }

  //  clearVisited();

  //  std::sort(conflict.begin(), conflict.end(), [&](const lit a, const lit b)
  //  {
  //    return decisionLevel(a) > decisionLevel(b);
  //  });

  ClauseAdded.trigger(conflict);

  assert(propagationLevel(conflict[0]) >= propagationLevel(decisions.back()));

  auto uip_lvl{propagationLevel(conflict[1])};
  int jump{1};
  while (propagationLevel(decisions[decisions.size() - 1 - jump]) > uip_lvl) {
    //        std::cout << decisions[decisions.size() - 1 - jump] << " ("
    //        << propagationLevel(decisions[decisions.size() - 1 - jump]) <<
    //        ") > "
    //        << uip_lvl << std::endl;
    ++jump;
    decisions.pop_back();
  }
  decisions.pop_back();

#ifdef DBG_TRACE
  if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
    std::cout << "backjump " << jump << " levels\n";
  }
#endif

  //  int max_level{
  //      (conflict.size() > 1 ? decisionLevel(conflict[1]) : init_level)};

  //  assert(upper(HORIZON) == ub - Gap<T>::epsilon());

  //#ifdef DBG_CL
  ////  if (++num_clauses > DBG_CL)
  ////    exit(1);
  //  if (cl_file != NULL) {
  //    *cl_file << "0 " << (conflict.size() + 1) << " 0 1 " << upper(HORIZON);
  //  }
  //#endif

  //  for (size_t i{0}; i < conflict.size(); ++i) {
  //
  //#ifdef DBG_CL
  //    if (cl_file != NULL) {
  //      writeLiteral(conflict[i]);
  //    }
  //#endif
  //
  //    conflict[i] = clauses.newNegLiteral(conflict[i]);
  //  }

  //#ifdef DBG_CL
  //  if (cl_file != NULL)
  //    *cl_file << std::endl;
  //#endif

  restoreState(env.level() - jump);

#ifdef DBG_TRACE
  auto cl =
#endif
      clauses.add(conflict.begin(), conflict.end(), true);

#ifdef DBG_TRACE
  if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
    if (clauses.size() > 0 and cl != NULL) {
      std::cout << "learn conflict" << *cl << std::endl;
    }
  }
#endif

#ifdef DBG_CL
  if (++num_clauses > DBG_CL)
    exit(1);
#endif
}

template <typename T> void Solver<T>::branchRight() {

  restoreState(env.level() - 1);
  //  lit deduction_var{search_vars[edge_propag_pointer]};
  auto deduction{~decisions.back()};
  decisions.pop_back();
  //  bool sign{not polarity[deduction_var]};

#ifdef DBG_TRACE
  if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
    std::cout << "-- backtrack to lvl " << env.level() << " & deduce "
              << deduction << " [i=" << num_choicepoints << "] --\n";
  }
#endif

  set(deduction);
}

template <typename T> boolean_state Solver<T>::search() {

  init_level = env.level();
  boolean_state satisfiability{Unknown};
  while (satisfiability == Unknown) {
    //        and not KillHandler::instance().signalReceived()) {

    ++num_choicepoints;
    try {
      propagate();

      // make a checkpoint
      saveState();

      // all resource constraints are accounted for => a solution has been found
      if (boolean_search_vars.empty()) {
        satisfiability = True;
      } else {
        ++num_choicepoints;

        var_t x = *(
            boolean_search_vars.begin()); // heuristic->nextChoicePoint(*this);
        //        lit d = valueHeuristic->choosePolarity(x, *this);
        Literal<T> d{boolean.getLiteral((random() % 2), x)};
        decisions.push_back(d);

#ifdef DBG_TRACE
        if (DBG_BOUND) {
          std::cout << "--- search node (lvl=" << env.level()
                    << ") [i=" << num_choicepoints << "] ---\n";
          std::cout << " ** take decision " << d << " **\n";
          printTrace();
        }
#endif

        set(d);
      }
    } catch (NewFailure<T> &f) {
      try {
        backtrack(f.reason);
        //        if (num_fails > restart_limit) {
        //          restart();
        //        }
      } catch (const SearchExhausted &f) {
        satisfiability = False;
      }
    }
  }

  return satisfiability;
}

template <typename T> void tempo::Solver<T>::propagate() {

  index_t p_index{static_cast<index_t>(propag_pointer)};

  while (not propagation_queue.empty() or trail.size() > p_index) {

    while (trail.size() > p_index) {
      ++num_literals;
      Literal<T> l{trail[p_index]};
      auto culprit{reason[p_index].expl};

#ifdef DBG_TRACE
      if (DBG_TRACE & QUEUE) {
        std::cout << "triggers for (" << l << ") b/c " << culprit->id() << "/"
                  << reason[p_index] << std::endl;
      }
#endif

      clauses.unit_propagate(l);

      //              std::cout << "propagate " << info_t(l) << std::endl;
      //
      //        std::cout << boolean_constraints << std::endl;

      const std::vector<int> &cons =
          (l.isNumeric() ? numeric_constraints[l] : boolean_constraints[l]);
      const std::vector<unsigned> &rank =
          (l.isNumeric() ? numeric_constraints.rank(l)
                         : boolean_constraints.rank(l));

      //              std::cout << cons.size() << std::endl;
      //              exit(1);

      for (auto i{cons.size()}; i-- > 0;) {
        auto c{constraints[cons[i]]};
        if (c->idempotent and culprit == c)
          continue;

#ifdef DBG_TRACE
        if (DBG_TRACE & QUEUE) {
          std::cout << " -" << *c << " (" << c->id() << ")" << std::endl;
        }
#endif

        propagation_queue.triggers(l, rank[i], cons[i]);
      }

      ++p_index;
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

  propag_pointer = p_index;
}

template<typename T>
int Solver<T>::saveState() {
//    if(propag_pointer < trail.size()) {
//        propag_pointer = trail.size();
//    }
assert(propag_pointer == static_cast<index_t>(trail.size()));

int lvl{env.level()};
env.save();
ReversibleObject::save();
return lvl;
//    return env.level();
}

template<typename T>
void Solver<T>::restoreState(const int l) {
  env.restore(l);
}

template <typename T> void Solver<T>::undo() {
  //  boolen_search_vars.setStart(boolean_var_pointer);
  //    boolen_search_vars.setStart(numeric_var_pointer);
  size_t n{propag_pointer};
  while (trail.size() > n) {
    auto l{trail.back()};
    if (l.isNumeric()) {
      numeric.undo(l);
    } else {
      boolean.undo(l);
    }
    trail.pop_back();
    reason.pop_back();
  }
}

template <typename T>
template <typename G>
void Solver<T>::update(const bool bounds, const int s, const G &neighbors) {

  const std::vector<T> &shortest_path{numeric.get(bounds)};

  changed.clear();
  changed.add(s);

#ifdef DBG_BELLMAN
  int max_iter = 1000;
  if (DBG_BELLMAN) {
    std::cout << core << "\nstart explore from " << s
              << (bounds == bound::lower ? " (backward)" : " (forward)")
              << std::endl;
  }
#endif

  while (not changed.empty()) {

    auto u{changed.front()};
    changed.pop_front();

#ifdef DBG_BELLMAN
    if (max_iter-- < 0)
      exit(1);
    if (DBG_BELLMAN) {
      std::cout << "pop " << u << " q=(";
      for (auto evt : changed)
        std::cout << " " << evt;
      std::cout << " )" << std::endl;
    }
#endif

    for (auto edge : neighbors[u]) {
      int v{edge};
      auto w{edge.label()};

      if (shortest_path[u] + w < shortest_path[v]) {

#ifdef DBG_BELLMAN
        if (DBG_BELLMAN) {
          std::cout << " * shorter path "
                    << (bounds == bound::lower ? "from " : "to ") << v
                    << std::endl;
        }
#endif

        if (v == s) {

#ifdef DBG_BELLMAN
          if (DBG_BELLMAN) {
            std::cout << " negative cyle\n";
          }
#endif

          throw NewFailure<T>(
              {this, static_cast<hint>(Literal<T>::index(bounds, s))});
        }
        setNumeric(Literal<T>(bounds, v, shortest_path[u] + w),
                   {this, static_cast<hint>(edge.stamp())}, false);

        //                               bt, v, shortest_path[u] + w,
        //                               {this, (edge.stamp() >= 0
        //                                           ? EDGE(edge.stamp())
        //                                           :
        //                                           BOUND(bounds.getIndex(LIT(u,
        //                                           bt))))});

        if (not changed.has(v))
          changed.add(v);
      }
#ifdef DBG_BELLMAN
      else if (DBG_BELLMAN) {
        std::cout << " ignore " << v << std::endl;
      }
#endif
    }
  }
}

template <typename T> void Solver<T>::post(NewConstraint<T> *con) {

  constraints.push_back(con);
  propagation_queue.resize(constraints.size());

  boolean_constraints.resize(std::max(2 * boolean.size(), numConstraint()));
  numeric_constraints.resize(std::max(2 * numeric.size(), numConstraint()));

  con->post(numConstraint() - 1);
}

template <typename T> void Solver<T>::relax(NewConstraint<T> *con) {

  if (boolean_constraints.indegree(con->id()) > 0)
    boolean_constraints.remove(con->id(), IN);
  if (numeric_constraints.indegree(con->id()) > 0)
    numeric_constraints.remove(con->id(), IN);
}

template <typename T>
template <typename X>
void Solver<T>::addToSearch(const X &x) {
  var_t var_id{x.id()};
  if (x.isNumeric()) {
    numeric_search_vars.reserve(var_id + 1);
    if (not numeric_search_vars.has(var_id))
      numeric_search_vars.add(var_id);
  } else {
    boolean_search_vars.reserve(var_id + 1);
    if (not boolean_search_vars.has(var_id))
      boolean_search_vars.add(var_id);
  }
}

template <typename T>
void Solver<T>::wake_me_on(const Literal<T> l, const int c) {
  if (l.isNumeric()) {
    numeric_constraints.add(l, c);
  } else {
    boolean_constraints.add(l, c);
  }
}

template <typename T>
void Solver<T>::xplain(const Literal<T>, const hint,
                       std::vector<Literal<T>> &) {}

template <typename T>
std::ostream &Solver<T>::print_reason(std::ostream &os, const hint) const {
  os << "shortest-path";
  return os;
}

template <typename T> double Solver<T>::looseness(const Literal<T> &l) const {
  if (l.isNumeric()) {
    auto lb{numeric.lower(l.variable())};
    auto ub{numeric.upper(l.variable())};

    if (l.sign() == bound::lower) {
      auto b{-l.value()};
      assert(b >= lb);
      assert(b <= ub);
      return static_cast<double>(b - lb) / static_cast<double>(ub - lb);
    } else {
      auto b{l.value()};
      assert(b >= lb);
      assert(b <= ub);
      return static_cast<double>(ub - b) / static_cast<double>(ub - lb);
    }
  } else if (l.hasSemantic()) {
    auto c{boolean.getEdge(l)};
    return static_cast<double>(numeric.upper(c.from) - numeric.lower(c.to) +
                               c.distance) /
           static_cast<double>(numeric.upper(c.to) - numeric.lower(c.from) -
                               c.distance);
  }
  return .5;
}

template <typename T>
std::ostream &Solver<T>::displayDomains(std::ostream &os) const {
  for (var_t x{0}; x < numeric.size(); ++x) {
    os << "x" << x << ": [" << numeric.lower(x) << ".." << numeric.upper(x)
       << "]\n";
  }
  return os;
}

template <typename T>
std::ostream &Solver<T>::displayBranches(std::ostream &os) const {
//    os << "branch:";
//      for(var_t x{0}; x<boolean.size(); ++x) {
//          if(boolean.isTrue(x)) {
//              os << " " << Literal<T>(true,x);
//          } else if(boolean.isFalse(x)) {
//              os << " " << Literal<T>(false,x);
//          }
//      }
size_t i{0};
size_t j{1};
while (j < numLiteral()) {
  auto l{getLiteral(j)};
  if (i < decisions.size() and decisions[i] == l) {
    os << " [" << l << "]";
    ++i;
  } else {
    os << " " << l;
  }
  ++j;
}
os << std::endl;

//  for (auto l : decisions) {
//    os << " " << l;
//    if (l.hasSemantic()) {
//      os << " (" << boolean.getEdge(l) << ")";
//    }
//  }
//  os << std::endl;
return os;
}

template <typename T>
std::ostream &Solver<T>::displayVariables(std::ostream &os) const {
  for (auto x : boolean_search_vars) {
    os << " b" << x;
  }
  os << std::endl;
  return os;
}

template <typename T>
std::ostream &Solver<T>::displayConstraints(std::ostream &os) const {
  return os;
}

template <typename T>
std::ostream &Solver<T>::display(std::ostream &os, const bool dom,
                                 const bool bra, const bool sva, const bool pre,
                                 const bool cla, const bool bgr, const bool ngr,
                                 const bool con, const bool trl) const {
  if (bra) {
    //    os << "decisions:";
    os << "branch:";
    displayBranches(os);
    //      os << "branch:";
    ////      for(var_t x{0}; x<boolean.size(); ++x) {
    ////          if(boolean.isTrue(x)) {
    ////              os << " " << Literal<T>(true,x);
    ////          } else if(boolean.isFalse(x)) {
    ////              os << " " << Literal<T>(false,x);
    ////          }
    ////      }
    //      size_t i{0};
    //      size_t j{1};
    //      while(j < numLiteral()) {
    //          auto l{getLiteral(j)};
    //          if(i<decisions.size() and decisions[i] == l) {
    //              os << " [" << l << "]";
    //              ++i;
    //          } else {
    //              os << " " << l;
    //          }
    //          ++j;
    //      }
    //      os << std::endl;
  }
  if (dom) {
    os << "domains:\n";
    displayDomains(os);
  }
  if (sva) {
    os << "search vars:\n";
    displayVariables(os);
  }
  if (pre) {
    os << "precedence graph:\n" << core << std::endl;
  }
  if (cla) {
    //    os << "clauses:\n" << clauses << std::endl;
  }
  if (bgr) {
    os << "boolean triggers:\n";
    for (var_t x{0}; x < static_cast<var_t>(boolean.size()); ++x) {

      auto l{Literal<T>::index(true, x)};
      if (boolean_constraints.has(l) and boolean_constraints.outdegree(l) > 0) {
        os << "b" << x << ":";
        for (auto c : boolean_constraints[l]) {
          os << " " << c;
        }
        os << std::endl;
      }

      l = Literal<T>::index(false, x);
      if (boolean_constraints.has(l) and boolean_constraints.outdegree(l) > 0) {
        os << "¬b" << x << ":";
        for (auto c : boolean_constraints[l]) {
          os << " " << c;
        }
        os << std::endl;
      }
    }
  }
  if (ngr) {
    os << "numeric triggers:\n";
    for (var_t x{0}; x < static_cast<var_t>(numeric.size()); ++x) {

      auto l{Literal<T>::index(bound::lower, x)};
      if (numeric_constraints.has(l) and numeric_constraints.outdegree(l) > 0) {
        os << "lb(x" << x << "):";
        for (auto c : numeric_constraints[l]) {
          os << " " << c;
        }
        os << std::endl;
      }

      l = Literal<T>::index(bound::upper, x);
      if (numeric_constraints.has(l) and numeric_constraints.outdegree(l) > 0) {
        os << "ub(x" << x << "):";
        for (auto c : numeric_constraints[l]) {
          os << " " << c;
        }
        os << std::endl;
      }
    }
  }
  if (trl) {
    os << "trail (" << numLiteral() << " literals):\n";
    index_t i{0};
    for (auto l : trail) {
      os << l << " b/c " << reason[i++] << std::endl;
    }
  }
  if (con) {
    os << "constraints:\n";
    for (auto c : constraints)
      os << *c << std::endl;
  }

  //
  //  os << boolean.size() << " boolean vars:\n";
  //  for (var_t x{0}; x < boolean.size(); ++x) {
  //    os << "b" << x << ": ";
  //    if (boolean.hasSemantic(x)) {
  //      if (boolean.isTrue(x)) {
  //        os << boolean.getEdge(true, x) << std::endl;
  //      } else if (boolean.isFalse(x)) {
  //        os << boolean.getEdge(false, x) << std::endl;
  //      } else {
  //        os << boolean.getEdge(true, x) << " or " << boolean.getEdge(false,
  //        x)
  //           << std::endl;
  //      }
  //    } else {
  //      if (boolean.isTrue(x)) {
  //        os << "true\n";
  //      } else if (boolean.isFalse(x)) {
  //        os << "false\n";
  //      } else {
  //        os << "undef\n";
  //      }
  //    }
  //  }
  //  os << numeric.size() << " numeric vars:\n";
  //  for (var_t x{0}; x < numeric.size(); ++x) {
  //    os << "x" << x << ": [" << numeric.lower(x) << ".." << numeric.upper(x)
  //       << "]\n";
  //  }
  //  os << numLiteral() << " literals:\n";
  //  index_t i{0};
  //  for (auto l : trail) {
  //    os << l << " b/c " << reason[i++] << std::endl;
  //  }
  //  os << " precedence graph:\n" << core << std::endl;
  //  os << "constraints:\n";
  //  for (auto c : constraints) {
  //    os << c->id() << ": " << *c << std::endl;
  //  }
  //  os << "numeric triggers:\n";
  //  for (var_t x{0}; x < static_cast<var_t>(numeric.size()); ++x) {
  //
  //    auto l{Literal<T>::index(bound::lower, x)};
  //    if (numeric_constraints.has(l) and numeric_constraints.outdegree(l) > 0)
  //    {
  //      os << "lb(x" << x << "):";
  //      for (auto c : numeric_constraints[l]) {
  //        os << " " << c;
  //      }
  //      os << std::endl;
  //    }
  //
  //    l = Literal<T>::index(bound::upper, x);
  //    if (numeric_constraints.has(l) and numeric_constraints.outdegree(l) > 0)
  //    {
  //      os << "ub(x" << x << "):";
  //      for (auto c : numeric_constraints[l]) {
  //        os << " " << c;
  //      }
  //      os << std::endl;
  //    }
  //  }
  //
  //  os << "boolean search variables:";
  //  for (auto x : boolean_search_vars) {
  //    os << " b" << x;
  //  }
  //  os << std::endl;
  //  os << "numeric search variables:";
  //  for (auto x : numeric_search_vars) {
  //    os << " x" << x;
  //  }
  //  os << std::endl;
  return os;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const Solver<T> &x) {
  return x.display(os);
}
}

#endif

