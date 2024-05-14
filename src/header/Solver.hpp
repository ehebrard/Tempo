
#ifndef _TEMPO_SOLVER_HPP
#define _TEMPO_SOLVER_HPP


#include "Literal.hpp"
#include "DirectedGraph.hpp"
//#include "ClauseBase.hpp"
#include "Constant.hpp"
#include "ConstraintQueue.hpp"
#include "DistanceConstraint.hpp"
//#include "Global.hpp"
//#include "Objective.hpp"
//#include "Restart.hpp"
//#include "TemporalNetwork.hpp"
//#include "constraints/DisjunctiveEdgeFinding.hpp"
//#include "constraints/EdgeConstraint.hpp"
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
    var_t newVar(const info_t s=Constant::NoSemantic);
    // declare a new Boolean variable with a semantic (disjunction)
    var_t newDisjunct(DistanceConstraint<T> d1, DistanceConstraint<T> d2);
    
    size_t size() const;
    
    void set(Literal<T> l);
    void undo(Literal<T> l);
        
    Literal<T> getLiteral(const bool s, const var_t x) const;
    
    
protected:
  Solver<T> &solver;

    std::vector<bool> polarity;
    
    std::vector<info_t> edge_index;
    
        std::vector<DistanceConstraint<T>> edges;
    
};

template <typename T>
Literal<T> BooleanStore<T>::getLiteral(const bool s, const var_t x) const {
    return Literal<T>(s, x, edge_index[x]-s);
}


template <typename T>
BooleanStore<T>::BooleanStore(Solver<T> &s) : solver(s) {}

template <typename T>
size_t BooleanStore<T>::size() const {
    return polarity.size() / 2;
}

template <typename T>
var_t BooleanStore<T>::newVar(const info_t s) {
    var_t x{static_cast<var_t>(size())};
    
    polarity.push_back(false);
    polarity.push_back(false);
    
    edge_index.push_back(s);
    
    return x;
}

template <typename T>
var_t BooleanStore<T>::newDisjunct(DistanceConstraint<T> d1, DistanceConstraint<T> d2) {
    var_t x{newVar(static_cast<info_t>(edges.size()+1))};
    edges.push_back(d1);
    edges.push_back(d2);
}

template <typename T>
void BooleanStore<T>::set(Literal<T> l) {
    polarity[l] = true;
}

template <typename T>
void BooleanStore<T>::undo(Literal<T> l) {
    polarity[l] = false;
}

template <typename T> bool BooleanStore<T>::isTrue(const var_t x) const {
    return polarity[Literal<T>::index(true,x)];
}

template <typename T> bool BooleanStore<T>::isFalse(const var_t x) const {
    return polarity[Literal<T>::index(false,x)];
}

template <typename T> bool BooleanStore<T>::isUndefined(const var_t x) const {
    return not (polarity[Literal<T>::index(true,x)] or polarity[Literal<T>::index(false,x)]);
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
    
    size_t size() const;

    // declare a new numeric variable
    var_t newVar();
    // declare a new numeric variable with temporal semantic (can be involved in disjunctions and precedences)
    var_t newTemporalVar();
    
    void set(Literal<T> l);
    void undo(Literal<T> l);
    
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
var_t NumericStore<T>::newVar() {
    var_t x{static_cast<var_t>(size())};
    
    bound[bound::lower].push_back(Constant::Infinity<T>);
    bound[bound::upper].push_back(Constant::Infinity<T>);
    
//    bound_index[bound::lower].emplace_back(std::initializer_list< index_t>{Constant::IndexOfMin});
//    bound_index[bound::upper].emplace_back(std::initializer_list< index_t>{Constant::IndexOfMax});
    
    bound_index[bound::lower].resize(size());
    bound_index[bound::upper].resize(size());
    bound_index[bound::lower].back().push_back(Constant::InfIndex);
    bound_index[bound::upper].back().push_back(Constant::InfIndex);
    
    return x;
}

template <typename T>
var_t NumericStore<T>::newTemporalVar() {
    var_t x{newVar()};
    solver.core.newVertex(x);
    return x;
}

template <typename T>
void NumericStore<T>::set(Literal<T> l) {
    auto s{l.sign()};
    auto v{l.variable()};
    bound[s][v] = l.value();
    bound_index[s][v].push_back(static_cast<stamp_t>(solver.numLiteral()-1));
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
bool NumericStore<T>::falsified(const Literal<T> l) const {
  l.value() <= -bound[~l.sign()][l.variable()];
}

template <typename T> bool NumericStore<T>::satisfied(const Literal<T> l) const {
  l.value() > -bound[~l.sign()][l.variable()];
}


template <typename T> class Solver : public ReversibleObject {

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
  void set(Literal<T> l, Explanation e = Constant::NoReason);
//@}

  /**
   * @name search
   */
  //@{
  int saveState();
  void restoreState(const int);
  void undo() override;
  //@}

  /**
   * @name printing and trace
   */
  //@{
  std::ostream &display(std::ostream &os) const;
  //@}

   
    BooleanStore<T> boolean;
    var_t newBoolean();
    var_t newDisjunct(DistanceConstraint<T>, DistanceConstraint<T>);
    
    NumericStore<T> numeric;
    var_t newNumeric();
    var_t newTemporal();
    
            // graph with all the known edges
            DirectedGraph<StampedLabeledEdge<T, Literal<T>>> core;
    
      
//      DifferenceLogicStore<T> precedences;
    
private:
  Options options;

  BacktrackEnvironment env;

    
    // the stack of Literals reprensenting all the changes so far
    std::vector<Literal<T>> trail;
    // the reason for each propagation event
    std::vector<Explanation> reason;
    
    
  /**
   * @name domains
   */
  //@{

  // a reversible pointer to the most recent preopagation event that is not yet
  // propagated
  Reversible<size_t> propag_pointer;
    Reversible<size_t>  var_pointer;
//  Reversible<stamp_t> num_propag_pointer;
  //@}

  /**
   * @name constraints
   */
  //@{
  // all the clauses (learnt or from the base problem)
  //    ClauseBase<T> clauses;
  // data structure used to implement the overall propagation (parameter is the
  // number of priority classes)
  ConstraintQueue<3> propagation_queue;
  // all of the posted constraints
  std::vector<Constraint *> constraints;
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
    : ReversibleObject(&env),
boolean(*this),
numeric(*this),
options(std::move(opt))
      //, clauses(*this)
      ,
      propag_pointer(0, &env), propagation_queue(constraints),
      boolean_constraint_network(&env), numeric_constraint_network(&env), search_vars(0, &env) {
          trail.emplace_back(Constant::NoVarx, Constant::Infinity<T>);
          reason.push_back(Constant::NoReason);
          seed(options.seed);
      }

template <typename T>
var_t Solver<T>::newBoolean() {
    auto x{boolean.newVar()};
    boolean_constraint_network.resize(std::max(numConstraint(), boolean.size()));
    return x;
}

template <typename T>
var_t Solver<T>::newDisjunct(DistanceConstraint<T> d1, DistanceConstraint<T> d2) {
    auto x{boolean.newDisjunct(d1,d2)};
    boolean_constraint_network.resize(std::max(numConstraint(), boolean.size()));
    return x;
}

template <typename T>
var_t Solver<T>::newNumeric() {
    auto x{numeric.newVar()};
    numeric_constraint_network.resize(std::max(numConstraint(), numeric.size()));
    return x;
}

template <typename T>
var_t Solver<T>::newTemporal() {
    auto x{numeric.newTemporalVar()};
    numeric_constraint_network.resize(std::max(numConstraint(), numeric.size()));
    return x;
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
void Solver<T>::set(Literal<T> l, Explanation e) {
    reason.push_back(e);
    trail.push_back(l);
    if(l.isNumeric()) {
        numeric.set(l);
    } else {
        boolean.set(l);
    }
    
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
    while(trail.size() > n) {
        auto l{trail.back()};
        
        std::cout << "undo " << l << std::endl;
        
        if(l.isNumeric()) {
            numeric.undo(l);
        } else {
            boolean.undo(l);
        }
        trail.pop_back();
    }
}

template <typename T> std::ostream &Solver<T>::display(std::ostream &os) const {
    std::cout << boolean.size() << " boolean vars:\n";
    for(var_t x{0}; x<boolean.size(); ++x) {
        std::cout << "b" << x << ": ";
        if(boolean.isTrue(x)) {
            std::cout << "true\n";
        } else if(boolean.isFalse(x)) {
            std::cout << "false\n";
        } else {
            std::cout << "undef\n";
        }
    }
    std::cout << numeric.size() << " numeric vars:\n";
    for(var_t x{0}; x<numeric.size(); ++x) {
        std::cout << "x" << x << ": [" << numeric.lower(x) << ".." << numeric.upper(x) << "]\n";
    }
    std::cout << numLiteral() << " literals:\n";
    index_t i{0};
    for(auto l : trail) {
        std::cout << l << " b/c " 
        << reason[i++]
        << std::endl;
    }
    std::cout << " precedence graph:\n" << core << std::endl;
  return os;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const Solver<T> &x) {
  return x.display(os);
}
}

#endif

