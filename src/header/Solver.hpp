
#ifndef _TEMPO_SOLVER_HPP
#define _TEMPO_SOLVER_HPP

#include "DirectedGraph.hpp"
#include "Literal.hpp"
#include "Model.hpp"
#include "Failure.hpp"
//#include "ClauseBase.hpp"
#include "Constant.hpp"
#include "ConstraintQueue.hpp"
#include "DistanceConstraint.hpp"
#include "Global.hpp"
//#include "Objective.hpp"
//#include "Restart.hpp"
//#include "TemporalNetwork.hpp"
//#include "constraints/DisjunctiveEdgeFinding.hpp"
#include "constraints/EdgeConstraint.hpp"
//#include "constraints/Transitivity.hpp"
//#include "heuristics/HeuristicManager.hpp"
//#include "util/Heap.hpp"
//#include "util/KillHandler.hpp"
#include "util/Options.hpp"
//#include "util/SubscribableEvent.hpp"


namespace tempo {

using stamp_t = uint32_t;

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
  const DistanceConstraint<T> &getEdge(const bool s, const var_t x) const;
  const DistanceConstraint<T> &getEdge(const Literal<T> l) const;

  bool hasSemantic(const var_t x) const;

protected:
  Solver<T> &solver;

  std::vector<bool> polarity;

  std::vector<info_t> edge_index;

  std::vector<DistanceConstraint<T>> edges;
};

template <typename T>
Literal<T> BooleanStore<T>::getLiteral(const bool s, const var_t x) const {
  return Literal<T>(s, x, edge_index[x] + s);
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

template <typename T>
size_t BooleanStore<T>::size() const {
  return polarity.size() / 2;
}

template <typename T> BooleanVar<T> BooleanStore<T>::newVar(const info_t s) {
  BooleanVar<T> x{static_cast<var_t>(size())};

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

template <typename T>
void BooleanStore<T>::set(Literal<T> l) {
  polarity[l] = true;
  if (l.hasSemantic()) {
    assert(l.constraint() == (edge_index[l.variable()] + l.sign()));
    solver.set(edges[l.constraint()],
               static_cast<index_t>(solver.numLiteral() - 1));
  }
  assert(l.hasSemantic() or edge_index[l.variable()] == Constant::NoSemantic);
}

template <typename T>
void BooleanStore<T>::undo(Literal<T> l) {
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

template <typename T> bool BooleanStore<T>::satisfied(const Literal<T> l) const {
  return polarity[l];
}

template <typename T> bool BooleanStore<T>::falsified(const Literal<T> l) const {
  return polarity[~l];
}


//template<typename T>
//class DifferenceLogicStore {
//
//public:
//    DifferenceLogicStore(Solver<T> &s);
//  ~DifferenceLogicStore() = default;
//
//private:
//
//    // graph with all the known edges
//    DirectedGraph<StampedLabeledEdge<T, Literal<T>>> core;
//
//    std::vector<DistanceConstraint<T> edges;
//
//};



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

  index_t litIndex(const bool s, const var_t x) const;

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
  std::vector<std::vector<stamp_t>> bound_index[2];
  // [for each Literal] pointer to the previous Literal of the same numeric
  // variable (useful for undoing and for searching implicants)
  //  std::vector<stamp_t> prev_bound;
  // a 'time stamp' for each literal (the index of the most recent Boolean
  // literal)
  //  std::vector<stamp_t> stamp;
  // used for clause minimization
  //  std::vector<bool> visited;
};

template <typename T>
NumericStore<T>::NumericStore(Solver<T> &s) : solver(s) {
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
    bound_index[s][v].push_back(static_cast<stamp_t>(solver.numLiteral() - 1));
//    return true;
//  }
//  return false;
}

template <typename T>
void NumericStore<T>::undo(Literal<T> l) {
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

//template <typename T> void NumericStore<T>::resize(const size_t n) {
//  if (n == 0)
//    return;
//
//  bound[bound::lower].resize(n, -INFTY);
//  bound[bound::upper].resize(n, INFTY);
//
//    bound_index[bound::lower].resize(n);
//    bound_index[bound::upper].resize(n);
//
//}

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
index_t NumericStore<T>::litIndex(const bool s, const var_t x) const {
  return bound_index[s][x].back();
}

template <typename T>
bool NumericStore<T>::falsified(const Literal<T> l) const {
  return -(l.value()) > bound[~(l.sign())][l.variable()];
}

template <typename T> bool NumericStore<T>::satisfied(const Literal<T> l) const {
  return l.value() >= bound[l.sign()][l.variable()];
}

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
  Literal<T> getLiteral(const stamp_t i) const;

  // get the most recent Literal that entails l
  Literal<T> getImplicant(const Literal<T> l) const;

  // get the index in the propagation queue of the last Literal involving
  // variable x
  stamp_t getPropagationLevel(const var_t x) const;

  void propagate();
  void set(Literal<T> l, const NewExplanation<T> &e = Constant::NewNoReason<T>);
  bool setNumeric(Literal<T> l,
                  const NewExplanation<T> &e = Constant::NewNoReason<T>);
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

  /**
   * @name search
   */
  //@{
  int saveState();
  void restoreState(const int);
  void undo() override;
  //@}

  void xplain(const Literal<T>, const hint,
              std::vector<Literal<T>> &) override; // {}
  std::ostream &print_reason(std::ostream &os,
                             const hint) const override; // { return os; }
  //    int getType() const override;// { return CYCLEEXPL; }

  /**
   * @name printing and trace
   */
  //@{
  std::ostream &display(std::ostream &os) const;
  //@}

  BooleanStore<T> boolean;
  BooleanVar<T> newBoolean();
  DisjunctVar<T> newDisjunct(const DistanceConstraint<T> &,
                             const DistanceConstraint<T> &);

  NumericStore<T> numeric;
  NumericVar<T> newNumeric();
  TemporalVar<T> newTemporal(const T offset = 0);
  Job<T> newJob(const T mindur = 0, const T maxdur = Constant::Infinity<T>);

  // graph with all the known edges
  DirectedGraph<StampedLabeledEdge<T, index_t>> core;

  //      DifferenceLogicStore<T> precedences;

private:
  Options options;

  BacktrackEnvironment env;

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
  Reversible<size_t> propag_pointer;
  Reversible<size_t> var_pointer;
  //  Reversible<stamp_t> num_propag_pointer;
  //@}

  /**
   * @name constraints
   */
  //@{
  // all the clauses (learnt or from the base problem)
  //    ClauseBase<T> clauses;
  // data structure used to implement the overall propagation (parameter is
  // the number of priority classes)
  NewConstraintQueue<T, 3> propagation_queue;
  // all of the posted constraints
  std::vector<NewConstraint<T> *> constraints;
  // dependency graph variables/constraints
  DirectedGraph<int> boolean_constraint_network;
  DirectedGraph<int> numeric_constraint_network;
  // @}

  /**
   * @name search @TODO: search only on Boolean variables for now
   */
  //@{
  // the set of variables remaining to fix
  SparseSet<var_t, Reversible<size_t>> search_vars;
  // current polarity of the Boolean variables
  //  std::vector<bool> polarity;
  //  // copy of the best solution so far
  //  std::vector<bool> best_solution;
  // level at which a variable has been decided
  std::vector<stamp_t> var_level;
  //@}

  // buffers
  SparseSet<> changed;

public:
  long unsigned int num_fails{0};
  long unsigned int num_choicepoints{0};
  long unsigned int num_backtracks{0};
  long unsigned int num_Literals{0};
  long unsigned int num_updates{0};
  long unsigned int num_clause_propagations{0};
  long unsigned int num_cons_propagations{0};
  long unsigned int num_edge_prunings{0};
  long unsigned int num_cons_prunings{0};
  long unsigned int learnt_size{0};
};

template <typename T>
Solver<T>::Solver(Options opt)
    : ReversibleObject(&env), boolean(*this), numeric(*this), core(&env),
      options(std::move(opt))
      //, clauses(*this)
      ,
      propag_pointer(0, &env), propagation_queue(constraints),
      boolean_constraint_network(&env), numeric_constraint_network(&env),
      search_vars(0, &env) {
  trail.emplace_back(Constant::NoVarx, Constant::Infinity<T>);
  reason.push_back(Constant::NewNoReason<T>);
  seed(options.seed);
}

template <typename T> BooleanVar<T> Solver<T>::newBoolean() {
  auto x{boolean.newVar()};
  boolean_constraint_network.resize(std::max(numConstraint(), boolean.size()));
  return x;
}

template <typename T>
DisjunctVar<T> Solver<T>::newDisjunct(const DistanceConstraint<T> &d1,
                                      const DistanceConstraint<T> &d2) {
  auto x{boolean.newDisjunct(d1, d2)};
  boolean_constraint_network.resize(std::max(numConstraint(), boolean.size()));

  post(new NewEdgeConstraint<T>(*this, boolean.getLiteral(true, x)));
  post(new NewEdgeConstraint<T>(*this, boolean.getLiteral(false, x)));

  return x;
}

template <typename T> NumericVar<T> Solver<T>::newNumeric() {
  auto x{numeric.newVar()};
  changed.reserve(numeric.size());
  numeric_constraint_network.resize(std::max(numConstraint(), numeric.size()));
  return x;
}

template <typename T> TemporalVar<T> Solver<T>::newTemporal(const T offset) {
  TemporalVar<T> x{newNumeric().id(), offset};
  //  auto x{numeric.newVar()};
  core.newVertex(x);
//  numeric_constraint_network.resize(std::max(numConstraint(), numeric.size()));
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

template <typename T>
Literal<T> Solver<T>::getLiteral(const stamp_t i) const {
  return trail[i];
}

template <typename T> void tempo::Solver<T>::propagate() {
  propag_pointer = trail.size();
}

template <typename T>
void Solver<T>::set(const DistanceConstraint<T> &c, const index_t r) {
  core.emplace_edge(c.from, c.to, c.distance, r);
  boundClosure(c.from, c.to, c.distance, r);
}

template <typename T>
void Solver<T>::set(Literal<T> l, const NewExplanation<T> &e) {
  if (l.isNumeric()) {
    setNumeric(l,e);
  } else {
    setBoolean(l,e);
  }
}

template <typename T>
bool Solver<T>::setNumeric(Literal<T> l, const NewExplanation<T> &e) {
  //
  //  //    std::cout << *this ;
  //  std::cout << "set " << l << ": " << numeric.satisfied(l) << std::endl;

  if (not numeric.satisfied(l)) {
    reason.emplace_back(e);
    trail.push_back(l);
    numeric.set(l);
//
//    std::cout << "sign: " << l.sign() << "/" << bound::lower << "/"
//              << bound::upper << std::endl;

    if (l.sign() == bound::upper) {
//
//      std::cout << "update w.r.t. ub(x" << l.variable() << ")\n";
//      std::cout << core << std::endl;

      update(bound::upper, l.variable(), core);
    } else {
//
//      std::cout << "update w.r.t. lb(x" << l.variable() << ")\n";
//      std::cout << core << std::endl;
//      std::cout << "here\n";

      update(bound::lower, l.variable(), core.backward());

//      std::cout << "there\n";
    }

    return true;
    }
    return false;
}

template <typename T>
void Solver<T>::setBoolean(Literal<T> l, const NewExplanation<T> &e) {
  assert(not boolean.satisfied(l));
  reason.emplace_back(e);
  trail.push_back(l);
    boolean.set(l);
}

template <typename T>
void Solver<T>::boundClosure(const var_t x, const var_t y, const T d,
                             const index_t r) {
  // closure w.r.t. 0 (0 -> x -(d)-> y -> 0)

  NewExplanation<T> e{this, static_cast<hint>(r)};
  if (r == Constant::NoIndex)
    e = Constant::NewNoReason<T>;

  //    std::cout << "lower: " << numeric.lower(y) << std::endl;
  //    std::cout << "BC: " << geq<T>(x, numeric.lower(y) - d) << std::endl;
  //

  // reduce the lower bound of x and precessor
  if (numeric.lower(y) != -Constant::Infinity<T>) {
//    std::cout << geq<T>(x, numeric.lower(y) - d).sign() << "/" << bound::lower
//              << std::endl;
    setNumeric(geq<T>(x, numeric.lower(y) - d), e);
  }
  //    {
  //    update(bound::lower, static_cast<int>(x), core.backward());
  //  }

  //    std::cout << "upper: "<< numeric.upper(x) << std::endl;
  //    std::cout << "BC: " << leq<T>(y, numeric.upper(x) + d) << std::endl;

  // increase the upper bound of y and successors
  if (numeric.upper(x) != Constant::Infinity<T>) {
//    std::cout << leq<T>(y, numeric.upper(x) + d) << "/" << bound::upper
//              << std::endl;
    setNumeric(leq<T>(y, numeric.upper(x) + d), e);
  }
  //  {
  //    update(bound::upper, static_cast<int>(y), core);
  //  }
}

template<typename T>
int Solver<T>::saveState() {
//    if(propag_pointer < trail.size()) {
//        propag_pointer = trail.size();
//    }
assert(propag_pointer == trail.size());

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
  search_vars.setStart(var_pointer);
  size_t n{propag_pointer};
  while (trail.size() > n) {
    auto l{trail.back()};
    if (l.isNumeric()) {
      numeric.undo(l);
    } else {
      boolean.undo(l);
    }
    trail.pop_back();
  }
}

template <typename T>
template <typename G>
void Solver<T>::update(const bool bounds, const int s, const G &neighbors) {

  const std::vector<T> &shortest_path{numeric.get(bounds)};

  //                                    int max_iter{1000};

  changed.clear();
    
//    std::cout << s << " / " << changed.capacity() << std::endl;
    
  changed.add(s);

#ifdef DBG_BELLMAN
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
                    {this, static_cast<hint>(edge.stamp())});

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

  boolean_constraint_network.resize(
      std::max(2 * boolean.size(), numConstraint()));
  numeric_constraint_network.resize(
      std::max(2 * numeric.size(), numConstraint()));

  con->post(numConstraint() - 1);
}

template <typename T> void Solver<T>::relax(NewConstraint<T> *con) {

  if (boolean_constraint_network.indegree(con->id()) > 0)
    boolean_constraint_network.remove(con->id(), IN);
  if (numeric_constraint_network.indegree(con->id()) > 0)
    numeric_constraint_network.remove(con->id(), IN);
}

template <typename T>
void Solver<T>::wake_me_on(const Literal<T> l, const int c) {
  if (l.isNumeric()) {
    numeric_constraint_network.add(l, c);
  } else {
    boolean_constraint_network.add(l, c);
  }
}

template <typename T>
void Solver<T>::xplain(const Literal<T>, const hint,
                       std::vector<Literal<T>> &) {}

template <typename T>
std::ostream &Solver<T>::print_reason(std::ostream &os, const hint) const {
  os << " shortest path";
  return os;
}

template <typename T> std::ostream &Solver<T>::display(std::ostream &os) const {
  os << boolean.size() << " boolean vars:\n";
  for (var_t x{0}; x < boolean.size(); ++x) {
    os << "b" << x << ": ";
    if (boolean.hasSemantic(x)) {
      if (boolean.isTrue(x)) {
        os << boolean.getEdge(true, x) << std::endl;
      } else if (boolean.isFalse(x)) {
        os << boolean.getEdge(false, x) << std::endl;
      } else {
        os << boolean.getEdge(true, x) << " or " << boolean.getEdge(false, x)
           << std::endl;
      }
    } else {
      if (boolean.isTrue(x)) {
        os << "true\n";
      } else if (boolean.isFalse(x)) {
        os << "false\n";
      } else {
        os << "undef\n";
      }
    }
  }
  os << numeric.size() << " numeric vars:\n";
  for (var_t x{0}; x < numeric.size(); ++x) {
    os << "x" << x << ": [" << numeric.lower(x) << ".." << numeric.upper(x)
       << "]\n";
  }
  os << numLiteral() << " literals:\n";
  index_t i{0};
  for (auto l : trail) {
    os << l << " b/c " << reason[i++] << std::endl;
  }
  os << " precedence graph:\n" << core << std::endl;
  os << "constraints:\n";
  for (auto c : constraints) {
    os << c->id() << ": " << *c << std::endl;
  }
  os << "numeric triggers:\n";
  for (var_t x{0}; x < static_cast<var_t>(numeric.size()); ++x) {

    auto l{Literal<T>::index(bound::lower, x)};
    if (numeric_constraint_network.has(l) and
        numeric_constraint_network.outdegree(l) > 0) {
      os << "lb(x" << x << "):";
      for (auto c : numeric_constraint_network[l]) {
        os << " " << c;
      }
      os << std::endl;
    }

    l = Literal<T>::index(bound::upper, x);
    if (numeric_constraint_network.has(l) and
        numeric_constraint_network.outdegree(l) > 0) {
      os << "ub(x" << x << "):";
      for (auto c : numeric_constraint_network[l]) {
        os << " " << c;
      }
      os << std::endl;
    }
  }
  os << "boolean triggers:\n";
  for (var_t x{0}; x < static_cast<var_t>(boolean.size()); ++x) {

    auto l{Literal<T>::index(true, x)};
    if (boolean_constraint_network.has(l) and
        boolean_constraint_network.outdegree(l) > 0) {
      os << "x" << x << ":";
      for (auto c : boolean_constraint_network[l]) {
        os << " " << c;
      }
      os << std::endl;
    }

    l = Literal<T>::index(false, x);
    if (boolean_constraint_network.has(l) and
        boolean_constraint_network.outdegree(l) > 0) {
      os << "¬x" << x << ":";
      for (auto c : boolean_constraint_network[l]) {
        os << " " << c;
      }
      os << std::endl;
    }
  }
  return os;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const Solver<T> &x) {
  return x.display(os);
}
}

#endif

