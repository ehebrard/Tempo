
#ifndef _TEMPO_SOLVER_HPP
#define _TEMPO_SOLVER_HPP


#include "Literal.hpp"
#include "ClauseBase.hpp"
//#include "Constant.hpp"
#include "ConstraintQueue.hpp"
//#include "DistanceConstraint.hpp"
//#include "Global.hpp"
#include "Objective.hpp"
#include "Restart.hpp"
//#include "TemporalNetwork.hpp"
#include "constraints/DisjunctiveEdgeFinding.hpp"
#include "constraints/EdgeConstraint.hpp"
#include "constraints/Transitivity.hpp"
#include "heuristics/HeuristicManager.hpp"
#include "util/Heap.hpp"
#include "util/KillHandler.hpp"
#include "util/Options.hpp"
#include "util/SubscribableEvent.hpp"


namespace tempo {

using index_t = uint32_t;
using Lit = Literal<T>;


//
//class NumericStore<T> {
//    
//private:
//    // the current bounds (repeated for efficient read access)
//    std::vector<T> bound;
//    
//    // pointers
//    std::vector<index_t> prev_bound;
//};
//
//
//class TemporalNetwork<T> {
//    
//};

template<typename T>
class Solver
{

public:
  /**
   * @name constructors
   */
  //@{
  Solver(Options opt);
  ~Solver();
  //@}

  /**
   * @name count accessors
   */
  //@{
  /// Total Number of variables
  size_t numVariable() const;
  /// Number of variable of the numeric type
    size_t numNumericVariable() const;
  /// Number of variables of the Boolean type
  size_t numBooleanVariable() const;
  /// Number of  clauses
  size_t numClause() const;
  /// Number of constraints
  size_t numConstraint() const;
  /// Number of  changes
  size_t numLiteral() const;
  //@}

  /**
   * @name modelling methods
   */
  //@{
  var_t newNumericVar();
  var_t newBooleanVar();
    var_t newDisjunctVar(const DistanceConstraint<T> &if_true,
                         const DistanceConstraint<T> &if_false);
    //@}

    /**
     * @name value accessors
     */
    //@{
  bool value(const var_t x) const;
  bool isTrue(const var_t x) const;
  bool isFalse(const var_t x) const;
  bool isUndefined(const var_t x) const;
  bool falsified(const Lit l) const;
  bool satisfied(const Lit l) const;

  T upper(const var_t) const;
  T lower(const var_t) const;
    //@}

    /**
     * @name literal accessors
     */
    //@{
    // get the literal corresponding to the i-th propagation event
    Lit getLiteral(const index_t i) const;
    
    // get the most recent literal that entails l
    Lit getImplicant(const Lit l) const;
    
    // get the index in the propagation queue of the last Literal implying variable x
    index_t getPropagationLevel(const var_t x) const;

  void set(const var_t x, const bool sign, const T val = 0,
           Explanation e = Constant::NoReason);
//@}


  private:
    Options options;

    BacktrackEnvironment env;
    
    /**
     * @name domains
     */
    //@{
    // [for each numeric signed_var] the current bounds (repeated for efficient read access)
    std::vector<T> bound;
    // [for each numeric signed_var] the current index in the 'propagation_events' stack
    std::vector<index_t> bound_index;
        // [for each literal] pointer to the previous literal of the same numeric variable (useful for undoing and for searching implicants)
        std::vector<index_t> prev_bound;
    
    // the stack of literals reprensenting all the changes so far
    std::vector<Lit> propagation_events;
    // the reason for each propagation event
    std::vector<Explanation> reason;
    // a reversible pointer to the most recent preopagation event that is not yet propagated
    Reversible<index_t> propag_pointer;
    // graph with all the known edges
    DirectedGraph<StampedLabeledEdge<T,lit>> core;
    //@}

    /**
     * @name constraints
     */
    //@{
    // all the clauses (learnt or from the base problem)
    ClauseBase<T> clauses;
    // data structure used to implement the overall propagation (parameter is the number of priority classes)
    ConstraintQueue<3> propagation_queue;
// all of the posted constraints
    std::vector<Constraint *> constraints;
// dependency graph variables/constraints
    DirectedGraph<int> constraint_network;
// @}
    
    /**
     * @name search
     */
    //@{
    // the set of variables remaining to fix
    SparseSet<var_t> search_vars;
    // current polarity @TODO: check if using front/back in search_vars to implement polarity is more efficient
    std::vector<bool> polarity;
    // copy of the best solution so far
    std::vector<bool> best_solution;
    // level at which a variable has been decided @TODO: not sure what is means for numeric variables right now
    std::vector<index_t> var_level;
    //@}

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

};

}

#endif

