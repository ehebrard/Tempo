/************************************************
 * Tempo Solver.hpp
 *
 * Copyright 2024 Emmanuel Hebrard
 *
 * Tempo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 * Tempo is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Tempo.  If not, see <http://www.gnu.org/licenses/>.
 *
 ***********************************************/

#ifndef _TEMPO_SOLVER_HPP
#define _TEMPO_SOLVER_HPP

#include <fstream>
#include <iostream>

#include "ClauseBase.hpp"
#include "Constant.hpp"
#include "ConstraintQueue.hpp"
#include "DirectedGraph.hpp"
#include "DistanceConstraint.hpp"
#include "Failure.hpp"
#include "Global.hpp"
#include "Literal.hpp"
#include "Model.hpp"
#include "Objective.hpp"
#include "Restart.hpp"
#include "constraints/Cardinality.hpp"
#include "constraints/Incrementality.hpp"
#include "constraints/CumulativeCheck.hpp"
#include "constraints/CumulativeEdgeFinding.hpp"
#include "constraints/CumulativeOverlapFinding.hpp"
#include "constraints/CumulativeTimetabling.hpp"
#include "constraints/DisjunctiveEdgeFinding.hpp"
#include "constraints/EdgeConstraint.hpp"
#include "constraints/FullTransitivity.hpp"
#include "constraints/PseudoBoolean.hpp"
#include "constraints/Transitivity.hpp"
#include "heuristics/LNS/relaxation_interface.hpp"
#include "heuristics/heuristic_factories.hpp"
#include "heuristics/impl/DecayingEventActivityMap.hpp"
#include "heuristics/impl/ActivityMap.hpp"
#include "util/KillHandler.hpp"
#include "util/Options.hpp"
#include "util/Profiler.hpp"
#include "util/SubscribableEvent.hpp"
#include "util/traits.hpp"

// #define LEARNING_RATE_STUFF true
//  #define DBG_SHRINK
// #define NEW_ANALYZE

namespace tempo {

//! T is the numeric variable domain type
template<typename T> class Solver;

template <typename T> struct ConflictSet;

template <typename T> struct TrailItem {
    Literal<T> lit;
    int lvl;
    Explanation<T> reason;
    
    operator Literal<T>() const { return lit; }
    const Explanation<T> &getExplanation() const { return reason; }
    int level() const { return lvl; }
};

//! Boolean variables and literals manager
/*!
 Responsible for:
 Storage (memory)
 Read/write access
 */
template<typename T>
class StaticBooleanStore {
    
public:
    /**
     * @name constructors
     */
    //@{
    StaticBooleanStore();
    ~StaticBooleanStore() = default;
    //@}
    
    /**
     * @name Boolean variable constructors
     */
    //@{
    // declare a new Boolean variable
    BooleanVar<T> newVar(const info_t s = Constant::NoSemantic);
    // declare a new Boolean variable with a semantic (disjunction)
    BooleanVar<T> newDisjunct(const DistanceConstraint<T> &d1,
                              const DistanceConstraint<T> &d2);
    //@}
    
    /**
     * @name utils
     */
    //@{
    // number of Boolean variables
    size_t size() const;
    void resize(const size_t n);
    
    // returns a literal from a sign and a variable
    Literal<T> getLiteral(const bool s, const var_t x) const;
    
    // returns the difference logic constraint corresponding to a sign and a
    // variable [or Constant::NoEdge<T> if the variable has no semantic]
    const DistanceConstraint<T> &getEdge(const bool s, const var_t x) const;
    
    // returns the difference logic constraint corresponding to a literal [or
    // Constant::NoEdge<T> if the literal has no semantic]
    const DistanceConstraint<T> &getEdge(const Literal<T> l) const;
    
    // returns the k-th difference logic constraint
    const DistanceConstraint<T> &getEdge(const index_t k) const;
    
    // whether variable x has a difference logic constraint attached
    bool hasSemantic(const var_t x) const;
    //@}
    
    const std::vector<info_t>& getEdgeInfo() {return edge_index;}
    
    // the list of difference logic constraints
    const std::vector<DistanceConstraint<T>>& getEdges() { return edges; }
    
protected:
    
    // the rank of the difference logic constraint in "edges" for each Boolean
    // variable
    std::vector<info_t> edge_index;
    
    // the list of difference logic constraints
    std::vector<DistanceConstraint<T>> edges;
    
};




//! Boolean variables and literals manager
/*!
 Responsible for:
 Storage (memory)
 Read/write access
 */
template<typename T>
class BooleanStore : public StaticBooleanStore<T>
{
    
public:
    /**
     * @name constructors
     */
    //@{
    BooleanStore(Solver<T> &s);
    ~BooleanStore() = default;
    void copy(StaticBooleanStore<T>& s);
    //@}
    
    /**
     * @name value accessors
     */
    //@{
    // value in the best solution!! (use 'equal' within search)
    bool value(const BooleanVar<T> x) const;
    
    // value in the current branch
    bool value(const var_t x) const;
    bool equal(const var_t x, const bool v) const;
    bool isTrue(const var_t x) const;
    bool isFalse(const var_t x) const;
    bool isUndefined(const var_t x) const;
    bool falsified(const Literal<T> l) const;
    bool satisfied(const Literal<T> l) const;
    //@}
    
    /**
     * @name Boolean variable constructors
     */
    //@{
    // declare a new Boolean variable
    BooleanVar<T> newVar(const info_t s = Constant::NoSemantic);
    // declare a new Boolean variable with a semantic (disjunction)
    BooleanVar<T> newDisjunct(const DistanceConstraint<T> &d1,
                              const DistanceConstraint<T> &d2);
    void reserveVarMemory();
    void resize(const size_t n);
    //@}
    
    /**
     * @name utils
     */
    //@{
    // set literal l to true
    void set(Literal<T> l);
    
    // unset literal l
    void undo(Literal<T> l);
    
    // returns the rank of literal l in the trail
    index_t litIndex(const Literal<T> l) const;
    
    // saves the current solution
    void saveSolution() { best_solution = polarity; }
    bool hasSolution() const { return not best_solution.empty(); }
    auto bestSolution() const noexcept -> const std::vector<bool> & { return best_solution; }
    //@}
    
#ifdef LEARNING_RATE_STUFF
    //@{
    // learning rate stuff
    double getLearningRate(const var_t x) const;
    void updateLearningRate(const var_t x);
    void updateActivity(const var_t x);
    //@}
#endif
    
protected:
    Solver<T> &solver;
    
    // [for each literal] the current polarity (x is undefined if x and ~x are
    // both false)
    std::vector<bool> polarity;
    
    // [for each literal] the polarity in the best solution
    std::vector<bool> best_solution;
    
    // the rank of each literal in the trail (Constant::NoIndex if the literal
    // is not on the trail)
    std::vector<index_t> propagation_stamp;
    
#ifdef LEARNING_RATE_STUFF
    // learning rate stuff
    double alpha{.4};
    
    // [for each variable] the number of times it participated to a conflict
    std::vector<long unsigned int> participated;
    
    // [for each variable] the number of conflicts when it was assigned
    std::vector<long unsigned int> assigned_at;
    
    // [for each variable] its current learning rate
    std::vector<double> learning_rate;
#endif
};


//! Numeric variables and literals manager
/*!
 Responsible for:
 Storage (memory)
 Read/write access
 */
template<typename T>
class StaticNumericStore {
    
public:
    /**
     * @name constructors
     */
    //@{
    StaticNumericStore();
    ~StaticNumericStore() = default;
    //@}
    
    /**
     * @name value accessors
     */
    //@{
    // for use in search
    bool falsified(const Literal<T> l) const;
    bool satisfied(const Literal<T> l) const;
    T upper(const var_t x) const;
    T lower(const var_t x) const;
    //@}
    
    /**
     * @name Boolean variable constructors
     */
    //@{
    // declare a new numeric variable
    NumericVar<T> newVar(const T b = Constant::Infinity<T>);
    //@}
    
    /**
     * @name utils
     */
    //@{
    size_t size() const;
    void resize(const size_t n);
    void set(Literal<T> l);
    const std::vector<T> &get(const int b) const;
    //@}
    
protected:
    // [for each numeric signed_var] the current bounds
    std::vector<T> bound[2];
    
};



//! Numeric variables and literals manager
/*!
 Responsible for:
 Storage (memory)
 Read/write access
 */
template<typename T>
class NumericStore : public StaticNumericStore<T> {
    
public:
    /**
     * @name constructors
     */
    //@{
    NumericStore(Solver<T> &s);
    ~NumericStore() = default;
    //@}
    
    /**
     * @name value accessors
     */
    //@{
    // solution (do not use in search)
    T solutionUpper(const NumericVar<T> x) const;
    T solutionLower(const NumericVar<T> x) const;
    
    //    // for use in search
    //    bool falsified(const Literal<T> l) const;
    //    bool satisfied(const Literal<T> l) const;
    //    T upper(const var_t x) const;
    //    T lower(const var_t x) const;
    //    //@}
    
    /**
     * @name Boolean variable constructors
     */
    //@{
    // declare a new numeric variable
    NumericVar<T> newVar(const T b = Constant::Infinity<T>);
    void reserveVarMemory();
    void resize(const size_t n);
    //    NumericVar<T> newVar(const T lb, const T ub);
    //@}
    
    /**
     * @name utils
     */
    //@{
    Literal<T> getLiteral(const bool s, const var_t x) const;
    Literal<T> previousBound(const Literal<T> l) const;
    
    index_t litIndex(const Literal<T> l) const;
    index_t lastLitIndex(const bool s, const var_t x) const;
    
    //    size_t size() const;
    
    void set(Literal<T> l);
    void undo(Literal<T> l);
    
    //    const std::vector<T> &get(const int b) const;
    //@}
    
    /**
     * @name debug
     */
    //@{
    std::ostream &displayLiteralTrail(std::ostream &os, const bool s,
                                      const var_t x) const {
        for (auto i : bound_index[s][x]) {
            os << " @" << i << ": " << solver.getLiteral(i);
        }
        return os;
    }
    //@}
    
    // saves the current solution
    void saveSolution() {
        best_solution[bound::lower] = StaticNumericStore<T>::bound[bound::lower];
        best_solution[bound::upper] = StaticNumericStore<T>::bound[bound::upper];
    }
    bool hasSolution() const { return not best_solution[bound::lower].empty(); }
    auto bestSolution(const int b) const noexcept -> const std::vector<T> & {
        return best_solution[b];
    }
    //@}
    
#ifdef LEARNING_RATE_STUFF
    //@{
    // learning rate stuff
    double getLearningRate(const var_t x) const;
    void updateLearningRate(const var_t x);
    void updateActivity(const var_t x);
    //@}
#endif
    
private:
    Solver<T> &solver;
    
    // [for each numeric signed_var] the current bounds
    //    std::vector<T> bound[2];
    
    // [for each numeric signed_var] the bound in the best solution
    std::vector<T> best_solution[2];
    
    // [for each numeric signed_var] the current index in the 'propagation_events'
    // stack
    std::vector<std::vector<index_t>> bound_index[2];
    
#ifdef LEARNING_RATE_STUFF
    // learning rate stuff
    double alpha{.4};
    
    // [for each variable] the number of times it participated to a conflict
    std::vector<std::vector<long unsigned int>> participated;
    
    // [for each variable] the number of conflicts when it was assigned
    std::vector<std::vector<long unsigned int>> assigned_at;
    
    // [for each variable] its current learning rate
    std::vector<double> learning_rate;
#endif
};


//! Explainer for literals from difference logic
template <typename T = int> class GraphExplainer : public Explainer<T> {
    
public:
    GraphExplainer(Solver<T> &s);
    void xplain(const Literal<T>, const hint, std::vector<Literal<T>> &) override;
    std::ostream &print_reason(std::ostream &os, const hint) const override;
    
private:
    Solver<T> &solver;
};

//! Explainer for contradictions from bound collapses
template <typename T = int> class BoundExplainer : public Explainer<T> {
    
public:
    BoundExplainer(Solver<T> &s);
    void xplain(const Literal<T>, const hint, std::vector<Literal<T>> &) override;
    std::ostream &print_reason(std::ostream &os, const hint) const override;
    //    int getType() const override;
    
private:
    Solver<T> &solver;
};

//! Solver
/*!
 All solving algorithms and methods
 */
template <typename T = int> class Solver : public ReversibleObject {
    
public:
    /**
     * @name constructors
     */
    //@{
    Solver();
    Solver(Options opt);
    ~Solver(); // = default;
    //@}
    
    /**
     * @name count accessors
     */
    //@{
    /// Number of constraints
    size_t numConstraint() const;
    /// Number of  changes
    size_t numLiteral() const;
    /// Number of  decisions
    size_t numDecision() const;
    /// Decision level
    int level() const {return env.level();}
    //@}
    
    /**
     * @name subscribable events
     */
    ///@{
    mutable SubscribableEvent<Literal<T>> ChoicePoint; ///< triggered on choicepoints
    //    mutable SubscribableEvent<const std::vector<Literal<T>> &>
    //    ClauseAdded; ///< triggered when a new clause is learned
    mutable SubscribableEvent<const Solver<T> &>
        ClauseAdded; ///< triggered when a new clause is learned
    mutable SubscribableEvent<Literal<T>>
    DeductionMade; ///< triggered when branching right
    mutable SubscribableEvent<Explanation<T> &>
    ConflictEncountered; ///< triggered when a conflict is encountered
    mutable SubscribableEvent<> BackTrackCompleted; ///< triggered after a successful backtrack
    mutable SubscribableEvent<bool> SearchRestarted; ///< triggered on restart
    //    mutable SubscribableEvent<> FailureDetected; ///< triggered on failure
    mutable SubscribableEvent<Solver<T> &> SolutionFound; ///< triggered when a solution is found
    mutable SubscribableEvent<const Solver<T> &> PropagationCompleted; ///< triggered after a successful propagation
    mutable SubscribableEvent<const Solver<T> &> PropagationInitiated; ///< triggered before propagation
    ///@}
    
    /**
     * @name modelling methods
     */
    //@{
    // create an internal Boolean variable and return a model object pointing to
    // it
    BooleanVar<T> newBoolean();
    // create an internal Boolean variable with a difference logic semantic,
    // post the channelling constraints, and return a model object pointing to
    // it
    BooleanVar<T> newDisjunct(const DistanceConstraint<T> &,
                              const DistanceConstraint<T> &);
    // returns the constant 0
    static NumericVar<T> zero() { return NumericVar<T>(Constant::K, 0); }
    // returns the constant true
    static BooleanVar<T> truism() { return BooleanVar<T>(0); }
    // create an internal numeric variable and return a model object pointing to it
    NumericVar<T> newNumeric(const T lb = -Constant::Infinity<T>,
                             const T ub = Constant::Infinity<T>);
    
protected:
    // force to not use views for constants
    NumericVar<T> _newNumeric_(const T lb = -Constant::Infinity<T>,
                               const T ub = Constant::Infinity<T>);
    
public:
    NumericVar<T> newConstant(const T k);
    NumericVar<T> newOffset(NumericVar<T> &x, const T k);
    //    // create an internal temporal variable and return a model object pointing
    //    // to it
    //    //    TemporalVar<T> newTemporal(const T offset = 0);
    //    //    NumericVar<T> newTemporal(const T offset = 0);
    //    // create the internal variables (depending on the type of Interval) and
    //    // return a model object pointing to them
    //    Interval<T> newInterval(const T mindur = 0,
    //                            const T maxdur = Constant::Infinity<T>,
    //                            const T earliest_start = -Constant::Infinity<T>,
    //                            const T latest_start = Constant::Infinity<T>,
    //                            const T earliest_end = -Constant::Infinity<T>,
    //                            const T latest_end = Constant::Infinity<T>,
    //                            const BooleanVar<T> opt = Constant::True
    //                            );
    
    Interval<T> maybe_between(const NumericVar<T> s, const NumericVar<T> e);
    Interval<T> maybe_continuefor(const NumericVar<T> s, const NumericVar<T> d);
    
    Interval<T> between(const NumericVar<T> s, const NumericVar<T> e);
    Interval<T> continuefor(const NumericVar<T> s, const NumericVar<T> d);
    
    Interval<T> between(const NumericVar<T> s, const NumericVar<T> e, const BooleanVar<T> opt);
    Interval<T> continuefor(const NumericVar<T> s, const NumericVar<T> d, const BooleanVar<T> opt);
    //@}
    
    /**
     * @name Literal and variable accessors
     */
    //@{
    // get the Literal corresponding to the i-th propagation event
    Literal<T> getLiteral(const index_t i) const;
    Explanation<T> getReason(const index_t i) const;
    int getLevel(const index_t i) const;
    
    // get the explanation for the i-th literal
    Explanation<T> getReason(const Literal<T> l) const;
    
    // get the index in the propagation queue of the last Literal involving
    // variable x (to be used parcimoniously for numeric lits, not so efficient)
    index_t propagationStamp(const Literal<T> l) const;
    
    // get the decision level of the last Literal involving variable x (calls
    // propagationStamp, so conditions apply)
    index_t decisionLevel(const Literal<T> l) const;
    
    // set literal l true with explanation e
    void set(Literal<T> l, const Explanation<T> &e = Constant::NoReason<T>);
    void post(Literal<T> l) {
        set(l);
        propagate();
    }
    
    // add a difference logic constraint (r is for explanation purpose, useless
    // externally)
    void set(const DistanceConstraint<T> &c, const index_t r = Constant::NoIndex);
    void post(const DistanceConstraint<T> &c);
    
    // add x to the list of variable that must be given a value
    //    template <typename X> void addToSearch(const X &x);
    void addToSearch(const BooleanVar<T> &x);
    void addToSearch(const NumericVar<T> &x);
    
    // some measure of tightness/looseness for bound literals or literals with a
    // difference logic semantic
    double looseness(const Literal<T> &l) const;
    //@}
    
    /**
     * @name Constraint posting
     */
    //@{
    // add a new constraint to the model
    void post(Constraint<T> *);
    
    // add a new expression-tree constraint to the model
    void post(BooleanVar<T> con);
    
    // remove a constraint from the model
    void relax(Constraint<T> *);
    
    // notify the solver that propagator of given id should be called when a
    // literal becomes true. returns the rank of the trigger for the constraint's viewpoint. beware, triggers are not repeated, hence more than one can have the same rank
    size_t wake_me_on(const Literal<T>, const int id);
    
    // create and post a new cardinality propagator
    template <typename ItVar>
    void postCardinality(const ItVar beg_var, const ItVar end_var,
                         const bool sign, const T bound);
    template <typename ItLit>
    void postCardinality(const ItLit beg_lit, const ItLit end_lit,
                         const T bound);
    
    template <typename ItLit, typename ItW>
    void postPseudoBoolean(const ItLit beg_lit, const ItLit end_lit, ItW w,
                           const T bound);
    
    // create and post a new edge-finding propagator
    template <concepts::typed_range<Interval<T>> Tasks>
    void postEdgeFinding(Interval<T> &schedule, Tasks &&taskRange, Matrix<Literal<T>> lits);
    
    // create and post a new precedence reasoning propagator
    template <concepts::typed_range<Interval<T>> Tasks>
    void postTransitivity(Interval<T> &schedule, Tasks &&taskRange, Matrix<Literal<T>> lits);
    
    // create and post a new full transitivity propagator
    template <typename ItRes>
    FullTransitivity<T> *postFullTransitivity(const ItRes beg_res,
                                              const ItRes end_res);
    
    // create and post a propagator that defines "solved" regions, for incrementality purpose
    template <concepts::typed_range<Interval<T>> Tasks>
    void postCumulativeIncrementality(Tasks &&tasks, Matrix<Literal<T>> lits);
    
    // create and post a checker based on intersection graph for the cumulative constraint
    template <concepts::typed_range<Interval<T>> Tasks, concepts::typed_range<NumericVar<T>> Demands>
    void postCumulative(const NumericVar<T> c, Tasks &&tasks, Demands &&demands, Matrix<Literal<T>> lits);
    
    // create and post the strong edge-finding propagator for the cumulative constraint
    template <typename ItTask, typename ItNVar>
    void postStrongEdgeFinding(const Interval<T> s, const NumericVar<T> c,
                               const ItTask beg_task, const ItTask end_task,
                               const ItNVar beg_dem, const bool tt,
                               Incrementality<T> *b, const int approx);
    
    // create and post the strong edge-finding propagator for the cumulative constraint
    template <typename ItTask, typename ItNVar>
    void postOverlapFinding(const Interval<T> s, const NumericVar<T> c,
                            const ItTask beg_task, const ItTask end_task,
                            const ItNVar beg_dem, Matrix<Literal<T>> lits);
    
    // create and post the time-tabling propagator for the cumulative constraint
    template <typename ItTask, typename ItNVar>
    void postTimetabling(const NumericVar<T> c, const ItTask beg_task,
                         const ItTask end_task, const ItNVar beg_dem);
    //@}
    
    /**
     * @name search
     */
    //@{
    void trigger(const Literal<T> l);
    void propagate();
    
    // compute the bound-consistent closure of numeric variables w.r.t. the
    // precedence graph, after adding y - x <= d
    void boundClosure(const var_t x, const var_t y, const T d, const index_t r);
    
    // update in direction bt in {bound::lower,bound::upper} in the graph G
    template <typename G>
    void update(const bool bt, const int s, const G &neighbors);
    
    // must be called before the first call to 'search()'
    void initializeSearch();
    
    template<heuristics::heuristic<T> H>
    void setBranchingHeuristic(H &&h);
    // record the backtrack-environment level when calling 'initializeSearch()'
    // level 0 is for ground truth
    // if there are assumptions, level 1 is for assumptions and level 2 and
    // above are for search else level 1 and above are for search
    int init_level{0};
    
    // index of the first literal that is not a ground truth
    index_t ground_stamp{1};
    // index of the first literal that is not an assumption nor a ground truth
    index_t assumption_stamp{1};
    
    void updateGroundTruth() {
        if(env.level() == 0) {
            ground_stamp = numLiteral();
        }
    }
    
    void saveSolution();
    
    /**
     * set literals without any checks. May throw
     * @tparam L literal range type
     * @param literals range containing literals
     */
    template <concepts::typed_range<Literal<T>> L>
    void makeAssumptions(const L &literals);
    void makeAssumption(const Literal<T> l);
    //@}
    
    /**
     * @name search
     */
    //@{
    boolean_state search();
    
    template <typename S> void optimize(S &objective);
    
    template <typename S, lns::relaxation_policy P>
    void largeNeighborhoodSearch(S &objective, P &&relaxationPolicy);
    
    boolean_state satisfiable();
    void minimize(const NumericVar<T> x);
    void maximize(const NumericVar<T> x);
    
    void restart(const bool on_solution = false);
    void backtrack(Explanation<T> &e);
    void branchRight();
    void learnConflict(Explanation<T> &e);
    void analyze(Explanation<T> &e, const bool only_boolean = false);
    void analyzeDecisions(Explanation<T> &e);
    void not_analyze(Explanation<T> &e);
    void decisionCut(Explanation<T> &e);
    void minimizeClause();
    template <typename Iter> void minimizeSlice(Iter beg, Iter stop);
    void shrinkClause(); // std::vector<Literal<T>> &c);
    template <typename Iter> bool shrinkSlice(Iter beg, Iter stop);
    
    // return the stamp of the literal that can replace p in the cut
    // * when the return value is ref_stamp, there is no replacement
    // * when the return value is 0, p is irrelevant (redundant)
    // * otherwise p can be replaced the literal at the returned value position
    // in trail
    index_t getRelevance(const std::pair<index_t, Literal<T>> &lit);
    index_t getRelevantBoundRec(const index_t ref_stamp, const Literal<T> p,
                                const index_t p_stamp, const int depth,
                                const index_t stamp);
    
    const SparseSet<var_t, Reversible<size_t>> &getBranch() const {
        return boolean_search_vars;
    }
    //@}
    
    /**
     * @name reversibility
     */
    //@{
    int saveState();
    void restoreState(const int);
    void undo() override;
    
    BacktrackEnvironment &getEnv() { return env; }
    const Options &getOptions() const { return options; }
    //@}
    
    /**
     * @name printing and trace
     */
    //@{
    std::ostream &display(std::ostream &os, const bool dom = true,
                          const bool bra = true, const bool sva = false,
                          const bool pre = true, const bool cla = true,
                          const bool bgr = false, const bool ngr = false,
                          const bool con = true, const bool trl = false) const;
    std::ostream &displayTrail(std::ostream &os) const;
    std::ostream &displayDomains(std::ostream &os) const;
    std::ostream &displayBranches(std::ostream &os) const;
    std::ostream &displayVariables(std::ostream &os) const;
    std::ostream &displayConstraints(std::ostream &os) const;
    std::ostream &displayPrecedences(std::ostream &os) const;
    
    std::ostream &displayProgress(std::ostream &os) const;
    std::ostream &displayHeader(std::ostream &os, const int width = 62) const;
    std::ostream &displaySummary(std::ostream &os, std::string msg) const;
    
    void check_clauses(const char* msg);
    void check_clause(const index_t i);
    
    std::string pretty(const Literal<T> l) const;
    //@}
    
#ifdef LEARNING_RATE_STUFF
    //@{
    // learning rate stuff
    void updateActivity(const Literal<T> l);
    //@}
#endif
    
private:
    // reversible strutures
    BacktrackEnvironment env;
    
    Reversible<size_t> boolean_size;
    Reversible<size_t> numeric_size;
    Reversible<size_t> constraint_size;
    //    std::vector<size_t> boolean_size_trail;
    //    std::vector<size_t> numeric_size_trail;
    //    std::vector<size_t> constraint_size_trail;
    
public:
    /**
     * @name domains
     */
    //@{
    // everything about Boolean variable
    BooleanStore<T> boolean;
    
    // everything about numeric variable
    NumericStore<T> numeric;
    
    // all the clauses (learnt or from the base problem)
    ClauseBase<T> clauses;
    
    // graph with all the known edges
    DirectedGraph<StampedLabeledEdge<T, index_t>> core;
    //@}
    
    /**
     * @name search variables
     */
    //@{
    SparseSet<var_t, Reversible<size_t>> boolean_search_vars;
    SparseSet<var_t, Reversible<size_t>> numeric_search_vars;
    //@}
    
    const std::vector<Literal<T>>& getDecisions() const {return decisions;}
    
private:
    // solver options
    Options options;
    
    // decision stack
//    std::vector<Literal<T>> decisions;
    ReversibleVector<Literal<T>> decisions;
    
    // the stack of Literals reprensenting all the changes so far
    std::vector<TrailItem<T>> trail;
    
    // a reversible pointer to the most recent preopagation event that is not
    // yet propagated
    Reversible<index_t> propag_pointer;
    
    /**
     * @name constraints
     */
    //@{
    // data structure used to implement the overall propagation (parameter is
    // the number of priority classes)
    ConstraintQueue<T, 3> propagation_queue;
    // all of the posted constraints
    std::vector<Constraint<T> *> constraints;
    // dependency graph variables/constraints
    DirectedGraph<int> boolean_constraints;
    DirectedGraph<int> numeric_constraints;
    // @}
    
    /**
     * @name search strategies
     */
    //@{
    heuristics::PolymorphicHeuristic<T> heuristic;
    RestartManager<Solver<T>> restartPolicy;
    // @}
    
    //  private:
    // explanation for bounds from difference logic
    GraphExplainer<T> graph_exp;
    
    // explanation for cound collapses
    BoundExplainer<T> bound_exp;
    
    // helper for Bellman-Ford
    SparseSet<> changed;
    
    
    // helper for search methods
    var_t objective_var{Constant::NoVar};
    
    // helper for 'initializeSearch()' (since it need to be called only once)
    bool initialized{false};
    
    bool searchCancelled = false;
    
    /**
     * @name helpers for conflict-analysis
     */
    //@{

    std::vector<Literal<T>> lit_buffer;
    std::vector<Literal<T>> learnt_clause;
    
    util::StopWatch stopWatch;
    
public:
  ConflictSet<T> cut;

  const std::vector<int> &getNumericScope(const int cons_idx) const {
    return numeric_constraints.backward()[cons_idx];
    }
    const std::vector<int>& getBooleanScope(const int cons_idx) const {
        return boolean_constraints.backward()[cons_idx];
    }

    const std::vector<Literal<T>> &lastLearnt() const { return learnt_clause; }
    std::vector<Literal<T>>::iterator begin_learnt() {
        return learnt_clause.begin();
    }
    std::vector<Literal<T>>::iterator end_learnt() {
        return learnt_clause.end();
    }
    
    void cancelSearch() noexcept {
        searchCancelled = true;
    }
    // @}
    
private:
    // specialisation to numeric literals
    void setNumeric(Literal<T> l, const Explanation<T> &e = Constant::NoReason<T>,
                    const bool do_update = true);
    
    // specialisation to Boolean literals
    void setBoolean(Literal<T> l,
                    const Explanation<T> &e = Constant::NoReason<T>);

#ifdef OLDVSIDS
    heuristics::impl::EventActivityMap *activityMap{nullptr};
#endif
    heuristics::impl::ActivityMap numericActivityMap{options.vsids_decay};
    heuristics::impl::ActivityMap booleanActivityMap{options.vsids_decay};
//#endif

  public:
#ifdef OLDVSIDS
    void setActivityMap(heuristics::impl::EventActivityMap *map) {
        activityMap = map;
    }
    //    heuristics::impl::EventActivityMap *getActivityMap() {
    //        return activityMap;
    //    }
    double getLiteralActivity(const Literal<T> l) const {
      return activityMap->get(l, *this);
    }
#endif
    heuristics::impl::ActivityMap &getNumericActivity() {
      return numericActivityMap;
    }
    heuristics::impl::ActivityMap &getBooleanActivity() {
      return booleanActivityMap;
    }
    double getActivity(const Literal<T> l) const {
        auto x{l.variable()};
        if(l.isNumeric()) {
            if(numericActivityMap.size() > static_cast<size_t>(x))
                return numericActivityMap[x];
        } else {
            if(booleanActivityMap.size() > static_cast<size_t>(x))
                return booleanActivityMap[x];
        }
        return heuristics::impl::ActivityMap::baseIncrement;
    }


    /**
     * @name statistics
     */
    //@{
    // cpu time recorded on the first call to 'initializeSearch()'
    double start_time;
    // total number of failures
    long unsigned int num_fails{0};
    // total number of restarts
    long unsigned int num_restarts{0};
    // total number of solutions
    long unsigned int num_solutions{0};
    // total number of branches in the search tree
    long unsigned int num_choicepoints{0};
    // total number of literal generated (i.e., pruning)
    long unsigned int num_literals{0};
    // number of calls to a queued propagator
    long unsigned int num_cons_propagations{0};
    // number of calls to a queued propagator
    long unsigned int num_unit_propagations{0};
    // average depth of the search tree
    double avg_fail_level{0};
    //@}
    
//public:
    /**
     * @name debug
     */
    //@{
    bool isAssertive(std::vector<Literal<T>> &conf) const;
    //    bool isAssertive();
    
#ifdef DBG_TRACE
    int debug_flag{0};
    void printTrace() const;
#endif
    
#ifdef DBG_CL
    std::ofstream *cl_file{NULL};
    int num_clauses{0};
    
    void writeLiteral(const Literal<T> l) const;
    void writeConflict() const;
    void writeClause() const;
    void writeExplanation(const Literal<T> l) const;
    void checkClauseLimit();
#endif
    
};

template <typename T>
struct ConflictSet : public std::vector<std::pair<index_t, Literal<T>>> {
    
    // whether the conflict is sorted (to avoid doing it needlessly)
    bool sorted{false};
    
    // for numeric vars : we keep the index of every var in the list of literal
    // (Constant::NoIndex otherwise)
    std::vector<index_t> conflict_index[2];
    
    // for boolean vars : whether the var is in the conflict
    std::vector<bool> in_conflict;
    
    //
    std::vector<bool> cache_;
    std::vector<index_t> cached_;
    
    // add a new literal to the conflict set (check if it is already directly
    // entailed by, or directly entails a previous literal)
    bool add(const Literal<T> l, const index_t s = Constant::NoIndex);
    bool add(const std::pair<index_t, Literal<T>> &p);
    
    void change(const Literal<T> l, const index_t i);
    
    template <typename Iter> // std::vector<std::pair<index_t,
    // Literal<T>>>::iterator
    void remove(Iter i);
    
    // clear and resize the data structures
    void clear(Solver<T> *solver);
    void sort();
    
    void get(std::vector<Literal<T>> &clause);
    
    void getCached(std::vector<Literal<T>> &lits, Solver<T> &solver) {
        for (auto p : *this) {
            if (not cached(p.first)) {
                lits.push_back(p.second);
            }
        }
        for (auto i : cached_) {
            lits.push_back(solver.getLiteral(i));
        }
    }
    
    index_t getConflictIndex(const Literal<T> l) const {
        return conflict_index[l.sign()][l.variable()];
    }
    
    void setConflictIndex(const Literal<T> l, index_t i) {
        conflict_index[l.sign()][l.variable()] = i;
    }
    
    bool cached(const index_t stamp) const { return cache_[stamp]; }
    
    void cache(const index_t stamp) {
        if (not cached(stamp)) {
            cache_[stamp] = true;
            cached_.push_back(stamp);
        }
    }
    
    void uncache(const index_t stamp) { cache_[stamp] = false; }
    
    size_t saveCache() { return cached_.size(); }
    void restoreCache(const size_t n) {
        while (cached_.size() > n) {
            cache_[cached_.back()] = false;
            cached_.pop_back();
        }
    }
    
    bool has(const Literal<T> p) const {
        if (p.isNumeric()) {
            return getConflictIndex(p) != Constant::NoIndex;
        } else {
            return in_conflict[p.variable()];
        }
    }
    
    bool entailed(const Literal<T> p) const {
        if (p.isNumeric()) {
            auto p_idx{getConflictIndex(p)};
            if (p_idx != Constant::NoIndex) {
                return this->operator[](p_idx).second.value() <= p.value_unsafe();
            }
            return false;
        } else {
            return in_conflict[p.variable()];
        }
    }
    
    int backtrackLevel(Solver<T> &solver);
    
    int glueScore(Solver<T> &solver);
    
    std::ostream &display(std::ostream &os) const;
    
    bool verifyIntegrity();
};

template <typename T>
std::ostream &ConflictSet<T>::display(std::ostream &os) const {
    for (auto l : *this) {
        os << " [" << l.second << "]@" << l.first;
    }
    return os;
}

template <typename T>
template <typename Iter>
void ConflictSet<T>::remove(Iter i) {
    i->first = 0;
    if (i->second.isNumeric()) {
        setConflictIndex(i->second, Constant::NoIndex);
    } else {
        in_conflict[i->second.variable()] = false;
    }
}

template <typename T>
void ConflictSet<T>::change(const Literal<T> l, const index_t s) {
    
    auto i{getConflictIndex(l)};
    assert(this->operator[](i).first > s);
    
    this->operator[](i) = {s, l};
    
    sorted &= (((i == 0) or (this->operator[](i - 1).first >= s)) and
               ((i == this->size() - 1) or (this->operator[](i + 1).first <= s)));
}

template <typename T>
bool ConflictSet<T>::add(const std::pair<index_t, Literal<T>> &p) {
    return add(p.second, p.first);
}

template <typename T>
bool ConflictSet<T>::add(const Literal<T> l, const index_t s) {
    if (l.isNumeric()) {
        auto i{getConflictIndex(l)};
        if (i != Constant::NoIndex) {
            auto value{this->operator[](i).second.value()};
            if (value > l.value()) {
                this->operator[](i) = {s, l};
                sorted &=
                (((i == 0) or (this->operator[](i - 1).first >= s)) and
                 ((i == this->size() - 1) or (this->operator[](i + 1).first <= s)));
            }
        } else {
            setConflictIndex(l, this->size());
            auto i{this->size()};
            this->emplace_back(s, l);
            sorted &= (((i == 0) or (this->operator[](i - 1).first >= s)));
            return true;
        }
    } else if (not in_conflict[l.variable()]) {
        auto i{this->size()};
        this->emplace_back(s, l);
        in_conflict[l.variable()] = true;
        
        sorted &= (((i == 0) or (this->operator[](i - 1).first >= s)));
        return true;
    }
    
    return false;
}

template <typename T> void ConflictSet<T>::clear(Solver<T> *solver) {
    
    while (not cached_.empty()) {
        cache_[cached_.back()] = false;
        cached_.pop_back();
    }
    cache_.resize(solver->numLiteral(), false);
    
    conflict_index[bound::lower].resize(solver->numeric.size(),
                                        Constant::NoIndex);
    conflict_index[bound::upper].resize(solver->numeric.size(),
                                        Constant::NoIndex);
    in_conflict.resize(solver->boolean.size());
    sorted = true;
    
    while (not this->empty()) {
        auto p{this->back()};
        if (p.second.isNumeric()) {
            setConflictIndex(p.second, Constant::NoIndex);
        } else {
            in_conflict[p.second.variable()] = false;
        }
        this->pop_back();
    }
}
template <typename T> void ConflictSet<T>::sort() {
    if (not sorted) {
        std::sort(this->begin(), this->end(), std::greater<>());
        index_t i{0};
        for (auto p : *this) {
            if (p.second.isNumeric() and p.first != 0) {
                setConflictIndex(p.second, i);
            }
            ++i;
        }
        sorted = true;
    }
}

template <typename T>
void ConflictSet<T>::get(std::vector<Literal<T>> &clause) {
    sort();
    for (auto p : *this) {
        if (p.first != 0 and has(p.second)) {
            clause.push_back(~(p.second));
        }
    }
}

template <typename T> int ConflictSet<T>::glueScore(Solver<T> &solver) {

  //    std::cout << "compute glue score of";
  //    for(auto p : *this) {
  //        std::cout << " <" << p.first << "|" << p.second << ">";
  //    }
  //    std::cout << std::endl;

  if (this->size() < 2) {
    return 1;
    } else {
        sort();
    }
    int g{0};
    int pl{-1};
    for(auto p : *this) {
        if(p.first != 0) {

          //            std::cout << " <" << p.first << "|" << p.second <<
          //            ">\n";

          auto l{solver.getLevel(p.first)};
          if (l != pl) {
            ++g;
            pl = l;
            }
        }
    }
    return g;
}

template <typename T> int ConflictSet<T>::backtrackLevel(Solver<T> &solver) {
    
//    assert(this->size() >= 2);
    
//    std::cout << "comp. bt level (sz=" << this->size() << ")\n";
    
    
    if (this->size() < 2) {
        return 0;
    } else {
        sort();
    }
//    auto bl{solver.decisionLevel(this->operator[](1).second)};
    
    int bl{0};
    for(auto p : *this | std::views::drop(1)) {
        if(p.first != 0) {
            bl = solver.getLevel(p.first);
            break;
        }
    }
    
    //    std::cout << "BT LEVEL = " << bl << " (LIT = " <<
    //    this->operator[](1).second << " / " << this->operator[](1).first <<
    //    ")\n";
    
    return bl;
}

template <typename T> bool ConflictSet<T>::verifyIntegrity() {
    size_t count_bool;
    for (size_t i{0}; i < this->size(); ++i) {
        auto l{this->operator[](i)};
        if (l.isNumeric()) {
            auto lidx{getConflictIndex()};
            if (lidx != Constant::NoIndex) {
                if (lidx != i) {
                    std::cout << "bug mapping numeric lit\n";
                    return false;
                }
            }
        } else {
            if (in_conflict[l.variable()])
                ++count_bool;
        }
    }
    for (size_t i{0}; 0 < in_conflict.size(); ++i) {
        if (not in_conflict[i])
            ++count_bool;
    }
    if (count_bool != in_conflict.size()) {
        std::cout << "bug counting boolean lit\n";
        return false;
    }
    for (var_t i{0}; 0 < static_cast<var_t>(conflict_index[0].size()); ++i) {
        for (auto s{0}; s < 2; ++s) {
            if (conflict_index[s][i] != Constant::NoIndex) {
                auto l{this->operator[](conflict_index[s][i])};
                if (l.variable() != i or l.sign() != s) {
                    std::cout << "bug mapping numeric lit (<-)\n";
                    return false;
                }
            }
        }
    }
    return true;
}

/*!
 GraphExplainer implementation
 */
template <typename T>
void GraphExplainer<T>::xplain(const Literal<T> l, const hint h,
                               std::vector<Literal<T>> &Cl) {
    
    if (l == Contradiction<T>) {
        
        auto s{Literal<T>::sgn(h)};
        auto x{Literal<T>::var(h)};
        auto end_cycle{x};
        
        do {
            auto le{solver.getLiteral(
                                      solver.getReason(solver.numeric.lastLitIndex(s, x)).the_hint)};
            
            if (le.isNumeric()) {
                x = le.variable();
                assert(s == le.sign());
            } else {
                
                Cl.push_back(le);
                auto c{solver.boolean.getEdge(le)};
                x = (s == bound::lower ? c.to : c.from);
            }
        } while (x != end_cycle);
        
    } else {
        
        auto r_idx{static_cast<index_t>(h)};
        auto le{solver.getLiteral(r_idx)};
        
        Cl.push_back(le);
        
        if (not le.isNumeric()) {
            auto c{solver.boolean.getEdge(le)};
            Cl.emplace_back(l.sign(), (l.sign() == bound::lower ? c.to : c.from),
                            l.value() - c.distance, detail::Numeric{});
        }
    }
}

template <typename T>
std::ostream &GraphExplainer<T>::print_reason(std::ostream &os,
                                              const hint) const {
    os << "precedence graph";
    return os;
}

template <typename T>
GraphExplainer<T>::GraphExplainer(Solver<T> &s) : solver(s) {}

template <typename T>
void BoundExplainer<T>::xplain(const Literal<T> l, const hint h,
                               std::vector<Literal<T>> &Cl) {
    
    if (l == Contradiction<T>) {
        var_t x{static_cast<var_t>(h)};
        
#ifdef DBG_TRACE
        if (DBG_CBOUND and (DBG_TRACE & LEARNING)) {
            std::cout << "explain contradiction: wipe out on numeric var x" << x
            << " in [" << solver.numeric.lower(x) << ".."
            << solver.numeric.upper(x) << "]" << std::endl;
        }
#endif
        
        auto lidx{solver.numeric.lastLitIndex(bound::lower, x)};
        auto uidx{solver.numeric.lastLitIndex(bound::upper, x)};
        
#ifdef DBG_TRACE
        if (DBG_CBOUND and (DBG_TRACE & LEARNING)) {
            
            std::cout << lidx  << " AND " << uidx << ": ";
            
            if(lidx > 0)
                std::cout << solver.getLiteral(lidx);
            else
                std::cout << "ground truth\n";
            
            std::cout << " AND ";
            
            if(uidx > 0)
                std::cout << solver.getLiteral(uidx);
            else
                std::cout << "ground truth\n";
            
            std::cout << std::endl;
        }
#endif
        
        Literal<T> le{Truism<T>};
        Literal<T> nle;
        Explanation<T> exp{Constant::NoReason<T>};
        
        if (lidx < uidx) {
            if(lidx > 0)
                le = solver.getLiteral(lidx);
            else
                nle = solver.getLiteral(uidx);
            exp = solver.getReason(uidx);
        } else {
            if(uidx > 0)
                le = solver.getLiteral(uidx);
            else
                nle = solver.getLiteral(lidx);
            exp = solver.getReason(lidx);
        }
        
#ifdef DBG_TRACE
        if (DBG_CBOUND and (DBG_TRACE & LEARNING)) {
            if(le != Truism<T>)
                std::cout << " => " << le << " AND " << exp << std::endl;
            else
                std::cout << " => " << nle << " AND " << exp << std::endl;
        }
#endif
        
        
        if(le != Truism<T>) {
            Cl.push_back(le);
            nle = ~le;
        }
        exp.explain(nle, Cl);
        
    } else {
        std::cout << "explain lit " << l << " due to constraint "
        << solver.boolean.getEdge(static_cast<index_t>(h)) << std::endl;
        
        exit(1);
    }
}

/*!
 BoundExplainer implementation
 */
template <typename T>
std::ostream &BoundExplainer<T>::print_reason(std::ostream &os,
                                              const hint) const {
    os << "collapsed bounds";
    return os;
}

template <typename T>
BoundExplainer<T>::BoundExplainer(Solver<T> &s) : solver(s) {}



/*!
 StaticBooleanStore implementation
 */

template <typename T>
Literal<T> StaticBooleanStore<T>::getLiteral(const bool s, const var_t x) const {
    return makeBooleanLiteral<T>(s, x, edge_index[x]);
}

template <typename T>
const DistanceConstraint<T> &StaticBooleanStore<T>::getEdge(const bool s,
                                                            const var_t x) const {
    return edges[edge_index[x] + s];
}

template <typename T>
const DistanceConstraint<T> &
StaticBooleanStore<T>::getEdge(const Literal<T> l) const {
    return edges[l.constraint()];
}

template <typename T>
const DistanceConstraint<T> &StaticBooleanStore<T>::getEdge(const index_t i) const {
    return edges[i];
}

template <typename T> bool StaticBooleanStore<T>::hasSemantic(const var_t x) const {
    return edge_index[x] >= Constant::SomeSemantic;
}

template <typename T> StaticBooleanStore<T>::StaticBooleanStore() {
    // I don't remember what's that for
    edges.push_back(Constant::NoEdge<T>);
    edges.push_back(Constant::NoEdge<T>);
    
    // this is to have the constants true/false
    auto x{newVar()};
}

template <typename T> size_t StaticBooleanStore<T>::size() const {
    //    return polarity.size() / 2;
    return edge_index.size();
}

template <typename T>
void StaticBooleanStore<T>::resize(const size_t n) {
    edge_index.resize(n, 0);
    edges.resize(2*n, Constant::NoEdge<T>);
}

template <typename T> BooleanVar<T> StaticBooleanStore<T>::newVar(const info_t s) {
    BooleanVar<T> x{static_cast<var_t>(size())};
    edge_index.push_back(s);
    return x;
}

template <typename T>
BooleanVar<T> StaticBooleanStore<T>::newDisjunct(const DistanceConstraint<T> &d1,
                                                 const DistanceConstraint<T> &d2) {
    info_t d{static_cast<info_t>(edges.size())};
    BooleanVar<T> x{newVar(d).id(), d};
    edges.push_back(d1);
    edges.push_back(d2);
    return x;
}



/*!
 BooleanStore implementation
 */

template <typename T>
index_t BooleanStore<T>::litIndex(const Literal<T> l) const {
    return propagation_stamp[l.variable()];
}

// deep copy of static boolean store
template <typename T> void BooleanStore<T>::copy(StaticBooleanStore<T>& s) {
    StaticBooleanStore<T>::edge_index = s.getEdgeIndex();
    StaticBooleanStore<T>::edges = s.getEdges();
    resize(s.size());
    set(makeBooleanLiteral<T>(true, 0));
}

template <typename T> BooleanStore<T>::BooleanStore(Solver<T> &s) : StaticBooleanStore<T>(), solver(s) {
    reserveVarMemory();
    set(makeBooleanLiteral<T>(true, 0));
}

template <typename T> void BooleanStore<T>::reserveVarMemory() {
    propagation_stamp.push_back(Constant::NoIndex);
    
    polarity.push_back(false);
    polarity.push_back(false);
    
#ifdef LEARNING_RATE_STUFF
    // learning rate stuff
    participated.push_back(0);
    assigned_at.push_back(0);
    learning_rate.push_back(0.0);
#endif
}

template <typename T> void BooleanStore<T>::resize(const size_t n) {
    
    //    StaticBooleanStore<T>::resize(n);
    
    propagation_stamp.resize(n,Constant::NoIndex);
    polarity.resize(2*n,false);
    
#ifdef LEARNING_RATE_STUFF
    // learning rate stuff
    participated.resize(n,0);
    assigned_at.resize(n,0);
    learning_rate.resize(n,0.0);
#endif
}

template <typename T> BooleanVar<T> BooleanStore<T>::newVar(const info_t s) {
    reserveVarMemory();
    return StaticBooleanStore<T>::newVar(s);
}

template <typename T>
BooleanVar<T> BooleanStore<T>::newDisjunct(const DistanceConstraint<T> &d1,
                                           const DistanceConstraint<T> &d2) {
    reserveVarMemory();
    return StaticBooleanStore<T>::newDisjunct(d1,d2);
}

template <typename T> void BooleanStore<T>::set(Literal<T> l) {
    propagation_stamp[l.variable()] = (solver.numLiteral() - 1);
    polarity[l] = true;
    
#ifdef LEARNING_RATE_STUFF
    // learning rate stuff
    auto x{l.variable()};
    assigned_at[x] = solver.num_fails;
    participated[x] = 0;
    //
#endif
    
    
    if (l.hasSemantic()) {
        assert(l.constraint() == (StaticBooleanStore<T>::edge_index[l.variable()] + l.sign()));
        auto e{StaticBooleanStore<T>::edges[l.constraint()]};
        if (e != Constant::NoEdge<T>)
            solver.set(e, static_cast<index_t>(solver.numLiteral() - 1));
    }
    assert(l.hasSemantic() or StaticBooleanStore<T>::edge_index[l.variable()] == Constant::NoSemantic);
}

#ifdef LEARNING_RATE_STUFF

template <typename T> double BooleanStore<T>::getLearningRate(const var_t x) const {
    return learning_rate[x];
}

// learning rate stuff
template <typename T> void BooleanStore<T>::updateLearningRate(const var_t x) {
    if(solver.num_fails == 0)
        return;
    learning_rate[x] *= (1.0 - alpha);
    learning_rate[x] +=
    (static_cast<double>(participated[x]) /
     static_cast<double>(solver.num_fails - assigned_at[x] + 1)) *
    alpha;
}

template <typename T> void BooleanStore<T>::updateActivity(const var_t x) {
    ++participated[x];
}
#endif

template <typename T> void BooleanStore<T>::undo(Literal<T> l) {
    polarity[l] = false;
    
#ifdef LEARNING_RATE_STUFF
    // learning rate stuff
    updateLearningRate(l.variable());
#endif
}

template <typename T> bool BooleanStore<T>::value(const BooleanVar<T> x) const {
    return best_solution[Literal<T>::index(true, x.id())];
}

template <typename T> bool BooleanStore<T>::value(const var_t x) const {
    return polarity[Literal<T>::index(true, x)];
}

template <typename T>
bool BooleanStore<T>::equal(const var_t x, const bool v) const {
    return polarity[Literal<T>::index(v, x)];
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




/*!
 StaticNumericStore implementation
 */
template <typename T> StaticNumericStore<T>::StaticNumericStore() {}

template <typename T> size_t StaticNumericStore<T>::size() const {
    return bound[bound::lower].size();
}

template <typename T>
void StaticNumericStore<T>::resize(const size_t n) {
    bound[bound::lower].resize(n, Constant::Infinity<T>);
    bound[bound::upper].resize(n, Constant::Infinity<T>);
}

template <typename T> NumericVar<T> StaticNumericStore<T>::newVar(const T b) {
    NumericVar<T> x{static_cast<var_t>(size())};
    
    bound[bound::lower].push_back(b);
    bound[bound::upper].push_back(b);
    
    return x;
}


template <typename T> void StaticNumericStore<T>::set(Literal<T> l) {
    auto s{l.sign()};
    auto v{l.variable()};
    
    //    std::cout << "set lit " << l << std::endl;
    //
    //    std::cout << "was " << bound[s][v] << " set to " << l.value() << std::endl;
    //
    //    assert(bound[s][v] > l.value());
    
    
    if(bound[s][v] > l.value())
        bound[s][v] = l.value();
}

template <typename T>
const std::vector<T> &StaticNumericStore<T>::get(const int b) const {
    return bound[b];
}

//template <typename T> void StaticNumericStore<T>::setLowerBound(const NumericVar<T> var, const T val) {
//    bound[bound::lower][var.id()] = -val;
//}
//
//template <typename T> void StaticNumericStore<T>::setUpperBound(const NumericVar<T> var, const T val) {
//    bound[bound::upper][var.id()] = val;
//}

template <typename T> T StaticNumericStore<T>::upper(const var_t x) const {
    return bound[bound::upper][x];
}

template <typename T> T StaticNumericStore<T>::lower(const var_t x) const {
    return -bound[bound::lower][x];
}

template <typename T>
bool StaticNumericStore<T>::falsified(const Literal<T> l) const {
    return -(l.value()) > bound[not l.sign()][l.variable()];
}

template <typename T>
bool StaticNumericStore<T>::satisfied(const Literal<T> l) const {
    return l.value() >= bound[l.sign()][l.variable()];
}


///*!
// NumericStore implementation
// */
//template <typename T> NumericStore<T>::NumericStore(Solver<T> &s) : solver(s) {}
//
//template <typename T> size_t NumericStore<T>::size() const {
//    return bound[bound::lower].size();
//}
//
//template <typename T>
//const std::vector<T> &NumericStore<T>::get(const int b) const {
//    return bound[b];
//}
//
//template <typename T> NumericVar<T> NumericStore<T>::newVar(const T b) {
//    NumericVar<T> x{static_cast<var_t>(size())};
//
//    bound[bound::lower].push_back(b);
//    bound[bound::upper].push_back(b);
//
//    bound_index[bound::lower].resize(size());
//    bound_index[bound::upper].resize(size());
//
//    bound_index[bound::lower].back().push_back(Constant::InfIndex);
//    bound_index[bound::upper].back().push_back(Constant::InfIndex);
//
//#ifdef LEARNING_RATE_STUFF
//    // learning rate stuff
//    participated.resize(size());
//    assigned_at.resize(size());
//    learning_rate.push_back(0.0);
//#endif
//
//    return x;
//}
//
//template <typename T> void NumericStore<T>::set(Literal<T> l) {
//    auto s{l.sign()};
//    auto v{l.variable()};
//
//    assert(bound[s][v] > l.value());
//
//#ifdef LEARNING_RATE_STUFF
//    // learning rate stuff
//    assigned_at[v].push_back(solver.num_fails);
//    participated[v].push_back(0);
//    //
//#endif
//
//    bound[s][v] = l.value();
//    bound_index[s][v].push_back(static_cast<index_t>(solver.numLiteral() - 1));
//}
//
//template <typename T> void NumericStore<T>::undo(Literal<T> l) {
//
//    auto s{l.sign()};
//    auto v{l.variable()};
//
//    bound_index[s][v].pop_back();
//    bound[s][v] = solver.getLiteral(bound_index[s][v].back()).value();
//
//#ifdef LEARNING_RATE_STUFF
//    // learning rate stuff
//    updateLearningRate(l.variable());
//#endif
//}
//
//#ifdef LEARNING_RATE_STUFF
//
//template <typename T>
//double NumericStore<T>::getLearningRate(const var_t x) const {
//  return learning_rate[x];
//}
//
//// learning rate stuff
//template <typename T> void NumericStore<T>::updateLearningRate(const var_t x) {
//  if (solver.num_fails == 0)
//    return;
//
//  learning_rate[x] *= (1.0 - alpha);
//
//  learning_rate[x] +=
//      (static_cast<double>(participated[x].back()) /
//       static_cast<double>(solver.num_fails - assigned_at[x].back() + 1)) *
//      alpha;
//
//  participated[x].pop_back();
//  assigned_at[x].pop_back();
//}
//
//template <typename T> void NumericStore<T>::updateActivity(const var_t x) {
//  ++participated[x].back();
//}
//#endif
//
//template <typename T> T NumericStore<T>::upper(const NumericVar<T> x) const {
//
//    assert(best_solution[bound::upper].size() == bound[bound::upper].size());
//
//    assert(best_solution[bound::upper].size() > x.id());
//
//    return best_solution[bound::upper][x.id()] + x.offset();
//}
//template <typename T> T NumericStore<T>::lower(const NumericVar<T> x) const {
//
//    assert(best_solution[bound::lower].size() == bound[bound::lower].size());
//
//    assert(best_solution[bound::lower].size() > x.id());
//
//    return -best_solution[bound::lower][x.id()] + x.offset();
//}
//
//template <typename T> T NumericStore<T>::upper(const var_t x) const {
//    return bound[bound::upper][x];
//}
//
//template <typename T> T NumericStore<T>::lower(const var_t x) const {
//    return -bound[bound::lower][x];
//}
//
//template <typename T>
//Literal<T> NumericStore<T>::getLiteral(const bool s, const var_t x) const {
//    return solver.getLiteral(bound_index[s][x].back());
//}
//
//template <typename T>
//Literal<T> NumericStore<T>::previousBound(const Literal<T> l) const {
//  auto i{bound_index[l.sign()][l.variable()].rbegin()};
//  assert(solver.getLiteral(*i) == l);
//  --i;
//  return solver.getLiteral(*i);
//}
//
//template <typename T>
//index_t NumericStore<T>::lastLitIndex(const bool s, const var_t x) const {
//    return bound_index[s][x].back();
//}
//
//// get the highest index in the literal stack of literal p directly entailing literal l
//template <typename T>
//index_t NumericStore<T>::litIndex(const Literal<T> l) const {
//  auto i{bound_index[l.sign()][l.variable()].rbegin()};
//  while (solver.getLiteral(*(i + 1)).value() <= l.value()) {
//    ++i;
//  }
//
//  return *i;
//}
//
//template <typename T>
//bool NumericStore<T>::falsified(const Literal<T> l) const {
//    return -(l.value()) > bound[not l.sign()][l.variable()];
//}
//
//template <typename T>
//bool NumericStore<T>::satisfied(const Literal<T> l) const {
//    return l.value() >= bound[l.sign()][l.variable()];
//}


/*!
 NumericStore implementation
 */
template <typename T> NumericStore<T>::NumericStore(Solver<T> &s) : solver(s) {}

template <typename T> void NumericStore<T>::reserveVarMemory() {
    resize(bound_index[bound::lower].size()+1);
}

template <typename T> void NumericStore<T>::resize(const size_t n) {
    
    //    StaticNumericStore<T>::resize(n);
    
    auto i{bound_index[bound::lower].size()};
    
    bound_index[bound::lower].resize(n);
    bound_index[bound::upper].resize(n);
    
    while(i < n) {
        bound_index[bound::lower][i].push_back(Constant::InfIndex);
        bound_index[bound::upper][i].push_back(Constant::InfIndex);
        ++i;
    }
    
#ifdef LEARNING_RATE_STUFF
    // learning rate stuff
    participated.resize(n);
    assigned_at.resize(n);
    learning_rate.resize(n,0.0);
#endif
}

template <typename T> NumericVar<T> NumericStore<T>::newVar(const T b) {
    reserveVarMemory();
    return StaticNumericStore<T>::newVar(b);
}

template <typename T> void NumericStore<T>::set(Literal<T> l) {
    auto s{l.sign()};
    auto v{l.variable()};
    
    assert(StaticNumericStore<T>::bound[s][v] > l.value());
    
#ifdef LEARNING_RATE_STUFF
    // learning rate stuff
    assigned_at[v].push_back(solver.num_fails);
    participated[v].push_back(0);
    //
#endif
    
    StaticNumericStore<T>::bound[s][v] = l.value();
    bound_index[s][v].push_back(static_cast<index_t>(solver.numLiteral() - 1));
}

template <typename T> void NumericStore<T>::undo(Literal<T> l) {
    
    auto s{l.sign()};
    auto v{l.variable()};
    
    bound_index[s][v].pop_back();
    StaticNumericStore<T>::bound[s][v] = solver.getLiteral(bound_index[s][v].back()).value();
    
#ifdef LEARNING_RATE_STUFF
    // learning rate stuff
    updateLearningRate(l.variable());
#endif
}

#ifdef LEARNING_RATE_STUFF

template <typename T>
double NumericStore<T>::getLearningRate(const var_t x) const {
    return learning_rate[x];
}

// learning rate stuff
template <typename T> void NumericStore<T>::updateLearningRate(const var_t x) {
    if (solver.num_fails == 0)
        return;
    
    learning_rate[x] *= (1.0 - alpha);
    
    learning_rate[x] +=
    (static_cast<double>(participated[x].back()) /
     static_cast<double>(solver.num_fails - assigned_at[x].back() + 1)) *
    alpha;
    
    participated[x].pop_back();
    assigned_at[x].pop_back();
}

template <typename T> void NumericStore<T>::updateActivity(const var_t x) {
    ++participated[x].back();
}
#endif

template <typename T> T NumericStore<T>::solutionUpper(const NumericVar<T> x) const {
    
    assert(best_solution[bound::upper].size() == StaticNumericStore<T>::bound[bound::upper].size());
    
    assert(best_solution[bound::upper].size() > x.id());
    
    return best_solution[bound::upper][x.id()] + x.offset();
}
template <typename T> T NumericStore<T>::solutionLower(const NumericVar<T> x) const {
    
    assert(best_solution[bound::lower].size() == StaticNumericStore<T>::bound[bound::lower].size());
    
    assert(best_solution[bound::lower].size() > x.id());
    
    return -best_solution[bound::lower][x.id()] + x.offset();
}

//template <typename T> T NumericStore<T>::upper(const var_t x) const {
//    return bound[bound::upper][x];
//}
//
//template <typename T> T NumericStore<T>::lower(const var_t x) const {
//    return -bound[bound::lower][x];
//}

template <typename T>
Literal<T> NumericStore<T>::getLiteral(const bool s, const var_t x) const {
    return solver.getLiteral(bound_index[s][x].back());
}

template <typename T>
Literal<T> NumericStore<T>::previousBound(const Literal<T> l) const {
    auto i{bound_index[l.sign()][l.variable()].rbegin()};
    assert(solver.getLiteral(*i) == l);
    --i;
    return solver.getLiteral(*i);
}

template <typename T>
index_t NumericStore<T>::lastLitIndex(const bool s, const var_t x) const {
    return bound_index[s][x].back();
}

// get the highest index in the literal stack of literal p directly entailing literal l
template <typename T>
index_t NumericStore<T>::litIndex(const Literal<T> l) const {
    auto i{bound_index[l.sign()][l.variable()].rbegin()};
    while (solver.getLiteral(*(i + 1)).value() <= l.value()) {
        ++i;
    }
    
    return *i;
}

//template <typename T>
//bool NumericStore<T>::falsified(const Literal<T> l) const {
//    return -(l.value()) > bound[not l.sign()][l.variable()];
//}
//
//template <typename T>
//bool NumericStore<T>::satisfied(const Literal<T> l) const {
//    return l.value() >= bound[l.sign()][l.variable()];
//}

/*!
 Solver implementation
 */
template <typename T>
Solver<T>::Solver()
: ReversibleObject(&env),
boolean_size(1, &env),
numeric_size(0, &env),
constraint_size(0, &env),
boolean(*this), numeric(*this), clauses(*this),
core(&env), boolean_search_vars(0, &env), numeric_search_vars(0, &env),
propag_pointer(1, &env), propagation_queue(constraints),
boolean_constraints(&env), numeric_constraints(&env),
restartPolicy(*this), graph_exp(*this), bound_exp(*this) {
    
    // sentinel literal for initial bounds
    Literal<T> l{Constant::NoVar, Constant::Infinity<T>, detail::Numeric{}};
    trail.emplace_back(l, level(), Constant::NoReason<T>);
    
    // pointed-to by all constants
    _newNumeric_(0, 0);
    seed(options.seed);
    
    post(&clauses);
    
    //          boolean_size_trail.push_back(0);
    //          numeric_size_trail.push_back(0);
    //          constraint_size_trail.push_back(0);
}

template <typename T> Solver<T>::~Solver() {
    for (auto c : constraints|std::views::drop(1)) {
        assert(c != &clauses);
        delete c;
    }
}

/*!
 Solver implementation
 */
template <typename T>
Solver<T>::Solver(Options opt)
: ReversibleObject(&env),
boolean_size(1, &env),
numeric_size(0, &env),
constraint_size(0, &env),
boolean(*this), numeric(*this), clauses(*this),
core(&env), boolean_search_vars(0, &env), numeric_search_vars(0, &env),
options(std::move(opt)), decisions(&env), propag_pointer(1, &env),
propagation_queue(constraints), boolean_constraints(&env),
numeric_constraints(&env), restartPolicy(*this), graph_exp(*this),
bound_exp(*this) {
    
    // sentinel literal for initial bounds
    Literal<T> l{Constant::NoVar, Constant::Infinity<T>, detail::Numeric{}};
    trail.emplace_back(l, level(), Constant::NoReason<T>);
    
    // pointed-to by all constants
    _newNumeric_(0, 0);
    seed(options.seed);
    
    post(&clauses);
    
    //          boolean_size_trail.push_back(0);
    //          numeric_size_trail.push_back(0);
    //          constraint_size_trail.push_back(0);
    
#ifdef DBG_CL
    if (options.dbg_file != "")
        cl_file = new std::ofstream(options.dbg_file, std::ofstream::out);
#endif
}

template <typename T> void Solver<T>::saveSolution() {
    
    //    std::cout << "save solution\n";
    
    boolean.saveSolution();
    numeric.saveSolution();
    
    //
    //    for(var_t x{0}; x<boolean.size(); ++x) {
    //        if(boolean.hasSemantic(x)) {
    //            auto e{boolean.getEdge(boolean.value(x), x)};
    //            std::cout << boolean.getLiteral(boolean.value(x), x) << " <=> " << e << std::endl;
    //
    //            if(e.falsified(*this)) {
    //                std::cout << "bug on b" << x << ": " << numeric.upper(e.to) << " - " << numeric.lower(e.from) << " > " << e.distance << std::endl;
    //                exit(1);
    //            } else {
    //                std::cout << "ok: " << numeric.upper(e.to) << " - " << numeric.lower(e.from) << " <= " << e.distance << std::endl;
    //            }
    //        }
    //    }
    
    ++num_solutions;
    SolutionFound.trigger(*this);
}

template <typename T> BooleanVar<T> Solver<T>::newBoolean() {
    auto x{boolean.newVar()};
    boolean_constraints.resize(std::max(numConstraint(), 2 * boolean.size()));
    clauses.newBooleanVar();
    ++boolean_size;
    return x;
}

template <typename T>
BooleanVar<T> Solver<T>::newDisjunct(const DistanceConstraint<T> &d1,
                                     const DistanceConstraint<T> &d2) {
    auto x{boolean.newDisjunct(d1, d2)};
    boolean_constraints.resize(std::max(numConstraint(), 2 * boolean.size()));
    clauses.newBooleanVar();
    if (d1 != Constant::NoEdge<T>)
        post(new EdgeConstraint<T>(*this, boolean.getLiteral(true, x.id())));
    
    if (d2 != Constant::NoEdge<T>)
        post(new EdgeConstraint<T>(*this, boolean.getLiteral(false, x.id())));
    ++boolean_size;
    return x;
}

template <typename T> NumericVar<T> Solver<T>::newConstant(const T k) {
    return NumericVar(Constant::K, k);
}

template <typename T>
NumericVar<T> Solver<T>::newOffset(NumericVar<T> &x, const T k) {
    return NumericVar<T>(x.id(), k);
}

template <typename T>
NumericVar<T> Solver<T>::_newNumeric_(const T lb, const T ub) {
    auto x{numeric.newVar(Constant::Infinity<T>)};
    
    changed.reserve(numeric.size());
    numeric_constraints.resize(std::max(numConstraint(), 2 * numeric.size()));
    core.newVertex(x.id());
    
    clauses.newNumericVar();
    
    set(geq<T>(x.id(), lb));
    set(leq<T>(x.id(), ub));
    propagate();
    
    
    ++numeric_size;
    return x;
}

template <typename T>
NumericVar<T> Solver<T>::newNumeric(const T lb, const T ub) {
    if (lb > ub) {
        throw Failure<T>({&bound_exp, Constant::NoHint});
    }
    else if (lb == ub) {
        return NumericVar(Constant::K, lb);
    }
    else {
        return _newNumeric_(lb,ub);
    }
}

//template <typename T>
//Interval<T> Solver<T>::newInterval(const T mindur, const T maxdur,
//                                   const T earliest_start, const T latest_start,
//                                   const T earliest_end, const T latest_end,
//                                   const BooleanVar<T> opt) {
//    return Interval<T>(*this, mindur, maxdur, earliest_start, latest_start, earliest_end, latest_end, opt);
//}


template <typename T>
Interval<T> Solver<T>::between(const NumericVar<T> s, const NumericVar<T> e) {
    Interval<T> i(s, e, e - s, BooleanVar<T>(Constant::True));
    post(i.duration >= 0);
    return i;
}

template <typename T>
Interval<T> Solver<T>::continuefor(const NumericVar<T> s, const NumericVar<T> d) {
    Interval<T> i(s, s+d, d, BooleanVar<T>(Constant::True));
    return i;
}

template <typename T>
Interval<T> Solver<T>::maybe_between(const NumericVar<T> s, const NumericVar<T> e) {
    Interval<T> i(s, e, e - s, newBoolean());
    post(i.duration >= 0);
    return i;
}

template <typename T>
Interval<T> Solver<T>::maybe_continuefor(const NumericVar<T> s, const NumericVar<T> d) {
    Interval<T> i(s, s+d, d, newBoolean());
    return i;
}

template <typename T>
Interval<T> Solver<T>::between(const NumericVar<T> s, const NumericVar<T> e, const BooleanVar<T> optional) {
    Interval<T> i(s, e, e - s, optional);
    post(i.duration >= 0);
    return i;
}

template <typename T>
Interval<T> Solver<T>::continuefor(const NumericVar<T> s, const NumericVar<T> d, const BooleanVar<T> optional) {
    Interval<T> i(s, s+d, d, optional);
    return i;
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
size_t Solver<T>::numDecision() const {
    return decisions.size();
}

template <typename T>
Explanation<T> Solver<T>::getReason(const Literal<T> l) const {
    return getReason(propagationStamp(l));
}

template <typename T> Literal<T> Solver<T>::getLiteral(const index_t i) const {
    return trail[i];
}

template <typename T>
Explanation<T> Solver<T>::getReason(const index_t i) const {
    return trail[i].getExplanation();
}

template <typename T> int Solver<T>::getLevel(const index_t i) const {
  if (i >= trail.size())
    return level() + 1;
  return trail[i].level();
}

template <typename T>
void Solver<T>::set(const DistanceConstraint<T> &c, const index_t r) {
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & PROPAGATION)) {
        std::cout << "set constraint " << c << " @" << numLiteral() << std::endl;
    }
#endif
    
    core.emplace_edge(c.from, c.to, c.distance, r);
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & PROPAGATION)) {
        if(numLiteral() == 70) {
            std::cout << "--> precedences:\n"; //<< core << std::endl;
            displayPrecedences(std::cout);
        }
    }
#endif
    
    boundClosure(c.from, c.to, c.distance, r);
}

template <typename T>
void Solver<T>::set(Literal<T> l, const Explanation<T> &e) {
    
    if (l.isNumeric()) {
        setNumeric(l, e);
    } else {
        setBoolean(l, e);
    }
}

template <typename T>
void Solver<T>::setNumeric(Literal<T> l, const Explanation<T> &e,
                           const bool do_update) {
    
    if (not numeric.satisfied(l)) {
        
#ifdef DBG_TRACE
        auto stamp{numLiteral()};
        if (DBG_BOUND and (DBG_TRACE & PROPAGATION)) {
            std::cout << "set " << pretty(l) << " @" << stamp << " b/c "
            << e  << (do_update ? "" : "*") << std::endl;
        }
#endif
        
        trail.emplace_back(l, level(), e);
        numeric.set(l);
        
        if (numeric.falsified(l)) {
            
#ifdef DBG_TRACE
            if (DBG_BOUND and (DBG_TRACE & PROPAGATION)) {
                std::cout << "failure* on " << pretty(l) << std::endl;
            }
#endif
            
#ifdef DBG_FAIL
            if (DBG_FAIL) {
                std::cout << "failure* on " << pretty(l) << std::endl;
            }
#endif
            
            throw Failure<T>({&bound_exp, static_cast<hint>(l.variable())});
        }
        
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
void Solver<T>::setBoolean(Literal<T> l, const Explanation<T> &e) {
    
    if(boolean.satisfied(l))
        return;
    
    if (boolean.falsified(l)) {
        
#ifdef DBG_TRACE
        if (DBG_BOUND and (DBG_TRACE & PROPAGATION)) {
            std::cout << "failure on " << pretty(l) << " @" << numLiteral()
            << " b/c " << e << std::endl;
        }
#endif
        
#ifdef DBG_FAIL
        if (DBG_FAIL) {
            std::cout << "failure on " << pretty(l) << " @" << numLiteral()
            << " b/c " << e << std::endl;
        }
#endif
        
        throw Failure<T>(e);
    }
    
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & PROPAGATION)) {
        if(numLiteral() == 69) {
            std::cout << "before set " << pretty(l) << " @" << numLiteral() << " b/c " << e
            << std::endl;
            
            std::cout << "--> precedences:\n"; //<< core << std::endl;
            displayPrecedences(std::cout);
        }
    }
#endif
    
    trail.emplace_back(l, level(), e);
    boolean.set(l);
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & PROPAGATION)) {
        std::cout << "set " << pretty(l) << " @" << numLiteral() << " b/c " << e
        << std::endl;
        if(numLiteral() == 70) {
            std::cout << "--> precedences:\n"; //<< core << std::endl;
            displayPrecedences(std::cout);
        }
    }
#endif
    
    if (boolean_search_vars.has(l.variable()))
        boolean_search_vars.remove_front(l.variable());
}

template <typename T>
void Solver<T>::boundClosure(const var_t x, const var_t y, const T d,
                             const index_t r) {
    // closure w.r.t. 0 (0 -> x -(d)-> y -> 0)
    
    Explanation<T> e{&graph_exp, static_cast<hint>(r)};
    if (r == Constant::NoIndex)
        e = Constant::NoReason<T>;
    
    if (numeric.lower(y) != -Constant::Infinity<T>) {
        setNumeric(geq<T>(x, numeric.lower(y) - d), e);
    }
    
    if (numeric.upper(x) != Constant::Infinity<T>) {
        setNumeric(leq<T>(y, numeric.upper(x) + d), e);
    }
}

template <typename T> void Solver<T>::restart(const bool on_solution) {
    
    ++num_restarts;
    env.restore(init_level);
//    decisions.clear();
    
    assert(decisions.empty());
    
    if (on_solution) {
        restartPolicy.initialize();
    } else {
        restartPolicy.reset();
    }
    
    SearchRestarted.trigger(on_solution);
    
    if (options.verbosity > Options::NORMAL) {
        std::cout << std::setw(8) << "restart(" << (on_solution ? "s" : "l") << ")";
        displayProgress(std::cout);
    }
    
}

template <typename T> void Solver<T>::backtrack(Explanation<T> &e) {
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
        std::cout << "failure @level " << env.level() << "/" << init_level
        << " b/c " << e << ":\n";
    }
#endif
    
    ConflictEncountered.trigger(e);
    propagation_queue.clear();
    
    try {
        if (options.learning)
            learnConflict(e);
        else
            branchRight();
    } catch (Failure<T> &f) {
        backtrack(f.reason);
    }
}

template <typename T>
index_t Solver<T>::propagationStamp(const Literal<T> l) const {
    if (l.isNumeric()) {
        return (l.variable() == Constant::K ? 0 : numeric.litIndex(l));
    } else {
        return boolean.litIndex(l);
    }
}

template <typename T>
index_t Solver<T>::decisionLevel(const Literal<T> p) const {
    return getLevel(propagationStamp(p));
}

template <typename T> void Solver<T>::minimizeClause() {
    if (cut.empty())
        return;
    cut.sort();
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
        std::cout << std::endl << "minimize " << cut << "\n";
    }
#endif
    
    cut.uncache(cut.begin()->first);
    minimizeSlice(cut.begin(), cut.end());
    //    minimizeSlice(cut.rbegin(), cut.rend());
}

template <typename T>
template <typename Iter>
void Solver<T>::minimizeSlice(Iter beg, Iter stop) {
    
    for (auto lit{beg}; lit != stop; ++lit) {
        
#ifdef DBG_TRACE
        if (DBG_BOUND and (DBG_TRACE & MINIMIZATION)) {
            std::cout << "* Try to minimize " << lit->second << " away\n";
        }
#endif
        
        index_t rstamp{getRelevance(*lit)};
        
#ifdef DBG_TRACE
        if (DBG_BOUND and (DBG_TRACE & MINIMIZATION)) {
            std::cout << "* relevant stamp = " << rstamp;
            if (rstamp < ground_stamp) {
                std::cout << " (remove)";
            } else if (rstamp < lit->first) {
                std::cout << " (change " << lit->second << " to " << getLiteral(rstamp)
                << ")";
            } else {
                std::cout << " (keep)";
            }
            std::cout << "\n";
        }
#endif
        
        if (rstamp < ground_stamp) {
            cut.remove(lit);
        } else if (rstamp < lit->first) {
            assert(lit->second.isNumeric());
            cut.change(getLiteral(rstamp), rstamp);
        }
    }
}

template <typename T>
index_t Solver<T>::getRelevance(const std::pair<index_t, Literal<T>> &lit) {
    return getRelevantBoundRec(lit.first, lit.second, lit.first,
                               options.minimization, 0);
}

template <typename T>
index_t Solver<T>::getRelevantBoundRec(const index_t l_stamp,
                                       const Literal<T> p,
                                       const index_t p_stamp, const int depth,
                                       const index_t stamp) {
    
    if (cut.cached(p_stamp)) {
#ifdef DBG_TRACE
        if (DBG_BOUND and (DBG_TRACE & MINIMIZATION)) {
            if (depth == options.minimization)
                std::cout << "cache @depth 0!!\n";
            for (auto d{0}; d < (options.minimization - depth); ++d)
                std::cout << "  ";
            std::cout << "** " << p << " is in redundant cache\n";
        }
#endif
        return 0;
    }
    
    index_t relevant_stamp{l_stamp};
    
    if (depth > 0) {
        auto r{getReason(p_stamp)};
        if (r != Constant::NoReason<T>) {
            relevant_stamp = stamp;
            
            auto beg_lit{lit_buffer.size()};
            r.explain(p, lit_buffer);
            auto end_lit{lit_buffer.size()};
            
            for (auto cur_lit{beg_lit}; cur_lit < end_lit; ++cur_lit) {
                auto q{lit_buffer[cur_lit]};
                auto q_lvl{propagationStamp(q)};
                
#ifdef DBG_TRACE
                if (DBG_BOUND and (DBG_TRACE & MINIMIZATION)) {
                    for (auto d{0}; d < (options.minimization - depth + 1); ++d)
                        std::cout << "  ";
                    std::cout << "is " << pretty(q) << " relevant?\n";
                }
#endif
                
                auto l{getLiteral(l_stamp)};
                if (q.sameVariable(l) and q.sign() == l.sign()) {
                    relevant_stamp = std::max(relevant_stamp, q_lvl);
                    if (relevant_stamp >= l_stamp) {
                        break;
                    }
#ifdef DBG_TRACE
                    else if (DBG_BOUND and (DBG_TRACE & MINIMIZATION)) {
                        for (auto d{0}; d < (options.minimization - depth + 1); ++d)
                            std::cout << "  ";
                        std::cout << "=> relevance = " << relevant_stamp << "\n";
                    }
#endif
                } else {
                    
                    if (q_lvl >= ground_stamp and not cut.entailed(q)) {
                        relevant_stamp = std::max(
                                                  relevant_stamp,
                                                  getRelevantBoundRec(l_stamp, q, q_lvl, depth - 1,
                                                                      relevant_stamp));
                        if (relevant_stamp >= l_stamp) {
                            break;
                        }
#ifdef DBG_TRACE
                        else if (DBG_BOUND and (DBG_TRACE & MINIMIZATION)) {
                            for (auto d{0}; d < (options.minimization - depth + 1);
                                 ++d)
                                std::cout << "  ";
                            std::cout << "=> relevance = " << relevant_stamp << "\n";
                        }
#endif
                    }
#ifdef DBG_TRACE
                    else if (DBG_BOUND and (DBG_TRACE & MINIMIZATION)) {
                        for (auto d{0}; d < (options.minimization - depth + 1); ++d)
                            std::cout << "  ";
                        std::cout << q << " is redundant";
                        if (q_lvl < ground_stamp)
                            std::cout << " (fact)";
                        if (cut.entailed(q))
                            std::cout << " (entailed)";
                        std::cout << "\n";
                    }
#endif
                }
            }
            lit_buffer.resize(beg_lit);
            
#ifdef DBG_TRACE
            if (DBG_BOUND and (DBG_TRACE & MINIMIZATION)) {
                for (auto d{0}; d < (options.minimization - depth); ++d)
                    std::cout << "  ";
                std::cout << p << " is ";
                if (relevant_stamp == 0)
                    std::cout << "completely redundant\n";
                else if (relevant_stamp < l_stamp) {
                    std::cout << "redundant if " << getLiteral(relevant_stamp) << "\n";
                } else {
                    std::cout << "relevant\n";
                }
            }
#endif
            
            if (relevant_stamp == 0) {
                
#ifdef DBG_TRACE
                if (DBG_BOUND and (DBG_TRACE & MINIMIZATION)) {
                    for (auto d{0}; d < (options.minimization - depth); ++d)
                        std::cout << "  ";
                    std::cout << "add " << p << " to redundant cache\n";
                }
#endif
                
                cut.cache(p_stamp);
            }
            
        }
#ifdef DBG_TRACE
        else if (DBG_BOUND and (DBG_TRACE & MINIMIZATION)) {
            for (auto d{0}; d < (options.minimization - depth); ++d)
                std::cout << "  ";
            std::cout << "decision\n";
        }
#endif
    }
#ifdef DBG_TRACE
    else if (DBG_BOUND and (DBG_TRACE & MINIMIZATION)) {
        for (auto d{0}; d < (options.minimization - depth); ++d)
            std::cout << "  ";
        std::cout << "max depth was reached\n";
    }
#endif
    
    return relevant_stamp;
}

template <typename T> void Solver<T>::decisionCut(Explanation<T> &) {
    
    std::cout << "decision cut!\n";
}


template <typename T>
void Solver<T>::shrinkClause() {
    
    if (cut.empty())
        return;
    cut.sort();
    
    //    std::cout << "hello\n";
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & SHRINKING)) {
        std::cout << std::endl << "++shrink " << cut << "\n";
        
        for (auto p : cut) {
            auto l{p.second};
            std::cout << std::setw(4) << decisionLevel(l) << " " << pretty(l)
            << std::endl;
        }
        std::cout << std::endl;
    }
#endif
    
    cut.uncache(cut.begin()->first);
    
    auto beg_slice{cut.rbegin()};
    auto lvl{getLevel(beg_slice->first)};
    auto num_lit{1};
    //  for (auto it{beg_slice + 1}; it != cut.rend(); ++it) {
    //    auto l{getLevel(it->first)};
    //    if (l != lvl) {
    //      if (num_lit > 2 and not shrinkSlice(beg_slice, it)) {
    //        minimizeSlice(beg_slice, it);
    //      }
    //      beg_slice = it;
    //      lvl = l;
    //        num_lit = 1;
    //    } else {
    //      ++num_lit;
    //    }
    //  }
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & SHRINKING)) {
        std::cout << "slice @lvl " << lvl << ":\n-- " << std::setw(4) << lvl << " "
        << std::setw(4) << beg_slice->first << ": " << beg_slice->second
        << std::endl;
    }
#endif
    
    //    auto tit{beg_slice+1};
    auto n{cut.size()};
    for (size_t i{n - 1}; i > 0; --i) {
        auto it{cut.rend() - i};
        auto l{getLevel(it->first)};
        
        if (l != lvl) {
            
#ifdef DBG_TRACE
            if (num_lit == 1 and DBG_BOUND and (DBG_TRACE & SHRINKING)) {
                std::cout << "skip " << beg_slice->second << " at level "
                << getLevel(beg_slice->first) << std::endl;
            }
#endif
            
            if (num_lit > 1 and not shrinkSlice(beg_slice, it)) {
                minimizeSlice(beg_slice, it);
            }
            beg_slice = (cut.rend() - i);
            lvl = l;
            num_lit = 1;
            
#ifdef DBG_TRACE
            if (DBG_BOUND and (DBG_TRACE & SHRINKING)) {
                std::cout << "slice @lvl " << l << ":\n-- " << std::setw(4)
                << getLevel(beg_slice->first) << " " << std::setw(4)
                << beg_slice->first << ": " << beg_slice->second << std::endl;
            }
#endif
            
        } else {
            ++num_lit;
            
#ifdef DBG_TRACE
            if (DBG_BOUND and (DBG_TRACE & SHRINKING)) {
                std::cout << "-- " << std::setw(4) << getLevel(it->first) << " "
                << std::setw(4) << it->first << ": " << it->second
                << std::endl;
            }
#endif
        }
    }
    
    //    auto tit{beg_slice+1};
    
    //    auto n{cut.size()};
    //  for (auto i{n - 1}; i-- > 0;) {
    //    auto l{getLevel(cut[i].first)};
    //
    // #ifdef DBG_TRACE
    //    if (DBG_BOUND and (DBG_TRACE & SHRINKING)) {
    //      std::cout << "-- " << std::setw(4) << l << " " << std::setw(4)
    //                << cut[i].first << ": " << cut[i].second << std::endl;
    //    }
    // #endif
    //
    //    if (l != lvl) {
    //        auto offset{(n - i - 1)};
    //      auto it{cut.rbegin() + offset};
    //
    // #ifdef DBG_TRACE
    //      if (num_lit == 1 and DBG_BOUND and (DBG_TRACE & SHRINKING)) {
    //        std::cout << "skip " << beg_slice->second << " at level "
    //                  << getLevel(beg_slice->first) << std::endl;
    //      }
    // #endif
    //
    //      if (num_lit > 1 and not shrinkSlice(beg_slice, it)) {
    //        minimizeSlice(beg_slice, it);
    //      }
    //      beg_slice = (cut.rbegin() + offset);
    //      lvl = l;
    //      num_lit = 1;
    //    } else {
    //      ++num_lit;
    //    }
    //  }
}

template <typename T>
template <typename Iter>
bool Solver<T>::shrinkSlice(Iter beg, Iter stop) {
    auto lvl{getLevel(beg->first)};
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & SHRINKING)) {
        std::cout << "shrink slice [" << beg->second << " -> " << (stop - 1)->second
        << "] @lvl" << lvl << std::endl;
    }
#endif
    
    int num_lit{0};
    
    auto waypoint{cut.saveCache()};
    for (auto it{beg}; it != stop; ++it) {
        
#ifdef DBG_TRACE
        if (DBG_BOUND and (DBG_TRACE & SHRINKING)) {
            std::cout << std::setw(4) << getLevel(it->first) << " " << std::setw(4)
            << it->first << ": " << it->second << std::endl;
        }
#endif
        
        cut.cache(it->first);
        ++num_lit;
    }
    
    auto lit_pointer = (stop - 1)->first;
    
    while (num_lit > 1) {
        
        auto l{getLiteral(lit_pointer)};
        auto exp{getReason(lit_pointer)};
        
#ifdef DBG_TRACE
        if (DBG_BOUND and (DBG_TRACE & SHRINKING)) {
            std::cout << " - resolve " << pretty(l) << " (#lits=" << num_lit
            << ") b/c " << exp << std::endl;
        }
#endif
        
        lit_buffer.clear();
        exp.explain(l, lit_buffer);
        
        //      std::cout << "lit_buffer.size() = " << lit_buffer.size() <<
        //      std::endl;
        
        --num_lit;
        
        for (auto p : lit_buffer) {
            
#ifdef DBG_TRACE
            if (DBG_BOUND and (DBG_TRACE & SHRINKING)) {
                std::cout << "  ** " << pretty(p) << ":";
            }
#endif
            
            auto p_stamp{propagationStamp(p)};
            auto p_lvl{getLevel(p_stamp)};
            
            if (p_stamp < ground_stamp) {
                // if p is a ground fact, we ignore it
#ifdef DBG_TRACE
                if (DBG_BOUND and (DBG_TRACE & SHRINKING)) {
                    std::cout << " => ground fact\n";
                }
#endif
                continue;
            } else if (p_lvl < lvl) {
                // p is not from the current decision level
                
                if (cut.entailed(p)) {
                    // p is entailed by the current conflict
#ifdef DBG_TRACE
                    if (DBG_BOUND and (DBG_TRACE & SHRINKING)) {
                        std::cout << " => entailed by cut\n";
                    }
#endif
                    continue;
                } else {
                    
#ifdef DBG_TRACE
                    if (DBG_BOUND and (DBG_TRACE & SHRINKING)) {
                        std::cout << " => from lvl " << p_lvl << ">1 (FAIL!)\n";
                    }
#endif
                    
                    cut.restoreCache(waypoint);
                    
#ifdef DBG_TRACE
                    if (DBG_BOUND and (DBG_TRACE & SHRINKING)) {
                        std::cout << "----\n";
                        for (auto p : cut) {
                            if (cut.has(p.second)) {
                                std::cout << std::setw(4) << p.first << " " << pretty(p.second)
                                << std::endl;
                            }
                        }
                        std::cout << "----\n";
                    }
#endif
                    
                    return false;
                }
                // p_stamp >= decision_stamp
            } else if (cut.cached(p_stamp)) {
                
#ifdef DBG_TRACE
                if (DBG_BOUND and (DBG_TRACE & SHRINKING)) {
                    std::cout << " => already in 'to-explore'\n";
                }
#endif
                
                continue;
                
            } else {
                
                // literal from the current decision level, and since there
                // are now at least two of them, we do not add it to the
                // conflict yet
                cut.cache(p_stamp);
                ++num_lit;
                
#ifdef DBG_TRACE
                if (DBG_BOUND and (DBG_TRACE & SHRINKING)) {
                    std::cout << " => add in 'to-explore':";
                    for (auto ii{lit_pointer}; --ii > 0;) {
                        if (cut.cached(ii)) {
                            std::cout << " [" << getLiteral(ii) << "]";
                        }
                    }
                    std::cout << std::endl;
                }
#endif
            }
        }
        
        while (not cut.cached(--lit_pointer))
            ;
    }
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & SHRINKING)) {
        std::cout << "successful shrinking: " << pretty(getLiteral(lit_pointer))
        << "\n";
    }
#endif
    
    auto novel{true};
    for (auto it{beg}; it != stop; ++it) {
        if (it->first != lit_pointer) {
            
#ifdef DBG_TRACE
            if (DBG_BOUND and (DBG_TRACE & SHRINKING)) {
                std::cout << " rm <" << it->first << ", " << it->second << ">\n";
            }
#endif
            
            cut.remove(it);
        } else {
            novel = false;
        }
    }
    if (novel) {
        cut.add(getLiteral(lit_pointer), lit_pointer);
        
#ifdef DBG_TRACE
        if (DBG_BOUND and (DBG_TRACE & SHRINKING)) {
            std::cout << "add new " << getLiteral(lit_pointer) << " in\n";
        }
#endif
        
    }
#ifdef DBG_TRACE
    else if (DBG_BOUND and (DBG_TRACE & SHRINKING)) {
        std::cout << "leave " << getLiteral(lit_pointer) << " in\n";
    }
    
    if (DBG_BOUND and (DBG_TRACE & SHRINKING)) {
        std::cout << "----\n";
        for (auto p : cut) {
            if (cut.has(p.second)) {
                std::cout << std::setw(4) << getLevel(p.first) << " " << std::setw(4)
                << p.first << " " << pretty(p.second) << std::endl;
            }
        }
        std::cout << "----\n";
    }
#endif
    
    //  for (auto it{beg}; it != stop; ++it) {
    //    cut.remove(it);
    //  }
    //
    //  //    std::cout << "cut.size(): " << cut.size() << std::endl;
    //  cut.add(getLiteral(lit_pointer), lit_pointer);
    cut.uncache(lit_pointer);
    
    return true;
}

template <typename T>
void Solver<T>::analyze(Explanation<T> &e, const bool only_boolean) {
    
    cut.clear(this);
    learnt_clause.clear();
    
    index_t decision_stamp;
    
    auto UIP{1 - only_boolean};
    
    //  /*
    //   We distinguish between standard conflict analysis and conflict analysis
    //   from a global fail (@decision level 0)
    //   */
    if (decisions.empty()) {
        
        if(assumption_stamp <= ground_stamp) {
            decision_stamp = ground_stamp;
            ground_stamp = 1;
        } else {
            decision_stamp = assumption_stamp;
        }
        
        
//        std::cout << "NO DECISIONS\n----------\n";
//        displayTrail(std::cout);
//        std::cout << "----------\n";
//        std::cout << "    ground_stamp = " << ground_stamp << std::endl;
//        std::cout << "assumption_stamp = " << assumption_stamp << std::endl;
//        std::cout << "  decision_stamp = " << decision_stamp << std::endl;
//        
        
        UIP = 0;
    } else {
        // marker to distinguish literal from the current decision level from the
        // rest
        decision_stamp = propagationStamp(decisions.back());
    }
    
    // number of literals of the current decision level in the conflict clause
    int num_lit{0};
    index_t lit_pointer{static_cast<index_t>(numLiteral())};
    Literal<T> l{Contradiction<T>};
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
        std::cout << "analyze conflict (" << e << "):\n";
    }
#endif
    
    // the first reason is the conflicting clause itself
    Explanation<T> &exp = e;
    do {
        
#ifdef DBG_TRACE
        if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
            if(l == Contradiction<T>)
                std::cout << "(#lit=" << num_lit << ") explain contradiction" << std::endl;
            else
                std::cout << "(#lit=" << num_lit << ") explore " << l << " @" << decisionLevel(l) << " b/c " << exp << std::endl;
        }
#endif
        
        if (exp == Constant::NoReason<T>) {
            // can happen when only_boolean = true (pass over the UIP)
            cut.add(l, lit_pointer);
            
#ifdef DBG_TRACE
            if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                std::cout << " => add to cut (decision): " << cut << std::endl;
            }
#endif
            
        } else {
            lit_buffer.clear();
            
            exp.explain(l, lit_buffer);
            
            
#ifdef DBG_CLPLUS
//            if (cl_file != NULL) {
//                *cl_file << "1 " << (lit_buffer.size() + (l != Contradiction<T>));
//                for (auto p : lit_buffer) {
//                    writeLiteral(p);
//                }
//                if (l != Contradiction<T>)
//                    writeLiteral(~l);
//                *cl_file << std::endl;
//            }
            writeExplanation(l);
            checkClauseLimit();
#endif
            
            for (auto p : lit_buffer) {
                
#ifdef DBG_TRACE
                if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                    std::cout << "   ** " << pretty(p) << " @" << decisionLevel(p);
                    std::cout.flush();
                }
#endif
                
                // rank of p in the trail
                auto p_stamp{propagationStamp(p)};
                
                if (p_stamp < ground_stamp) {
                    // if p is a ground fact, we ignore it
#ifdef DBG_TRACE
                    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                        std::cout << " => ground fact\n";
                    }
#endif
                    continue;
                } else if ((not(only_boolean and p.isNumeric())) and
                           p_stamp < decision_stamp) {
                    // p is not from the current decision level
                    
                    if (cut.entailed(p)) {
                        // p is entailed by the current conflict
#ifdef DBG_TRACE
                        if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                            std::cout << " => entailed by cut\n";
                        }
#endif
                        continue;
                    } else {
                        
                        cut.add(p, p_stamp);
                        
#ifdef DBG_TRACE
                        if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                            std::cout << " => add to cut: " << cut << std::endl;
                        }
#endif
                    }
                    // p_stamp >= decision_stamp
                } else if (cut.cached(p_stamp)) {
                    
#ifdef DBG_TRACE
                    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                        std::cout << " => already in 'to-explore'\n";
                    }
#endif
                    
                    continue;
                    
                } else {
                    
                    // literal from the current decision level, and since there
                    // are now at least two of them, we do not add it to the
                    // conflict yet
                    cut.cache(p_stamp);
                    ++num_lit;
                    
#ifdef DBG_TRACE
                    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                        std::cout << " => add in 'to-explore':";
                        for (auto ii{lit_pointer}; --ii > 0;) {
                            if (cut.cached(ii)) {
                                std::cout << " [" << getLiteral(ii) << "]";
                            }
                        }
                        std::cout << std::endl;
                    }
#endif
                }
            }
        }
        
        if (num_lit >= 1) {
            while (not cut.cached(--lit_pointer))
                ;
            l = getLiteral(lit_pointer);
            exp = getReason(lit_pointer);
        }
        
        --num_lit;
    } while (num_lit >= UIP and (only_boolean or lit_pointer >= decision_stamp));
    
    // l i the UIP, but to add it to the conflict, we must do the numeric test
    if (l != Contradiction<T> and not only_boolean and cut.add(l, lit_pointer)) {
        cut.uncache(lit_pointer);
    }
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
        std::cout << " => add UIP to cut: " << cut << std::endl;
    }
#endif
    
#ifdef DBG_TRACE
    auto size_before{cut.size()};
#endif
    
#ifdef DBG_CL
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
        std::vector<Literal<T>> buffer;
        for(auto p : cut) {
            if(p.first != 0) {
                buffer.push_back(~p.second);
            }
        }
        std::cout << "\nlearn clause of size " << buffer.size() << std::endl;
        for (auto p : buffer) {
            std::cout << std::setw(4) << decisionLevel(~p) << " " << pretty(p)
            << std::endl;
        }
        std::cout << std::endl;
    }
#endif
    writeConflict();
    checkClauseLimit();
#endif
    
    if (options.shrinking) {
        shrinkClause();
    } else if (options.minimization > 0) {
        minimizeClause();
    }
    
    cut.get(learnt_clause);
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
        std::cout << "\nlearn clause of size " << learnt_clause.size();
        if (size_before > learnt_clause.size())
            std::cout << " (minimized from " << size_before << ")";
        std::cout << std::endl;
        
        for (auto l : learnt_clause) {
            std::cout << std::setw(4) << decisionLevel(~l) << " " << pretty(l)
            << std::endl;
        }
        std::cout << std::endl;
    }
#endif
    
    //    shrinkClause();
    
#ifdef DBG_CL
    writeClause();
    checkClauseLimit();
#endif
    
    assert(static_cast<int>(decisionLevel(~(learnt_clause[0]))) == level());
}

template <typename T> void Solver<T>::analyzeDecisions(Explanation<T> &e) {
    
    cut.clear(this);
    learnt_clause.clear();
    
    //  /*
    //   We distinguish between standard conflict analysis and conflict analysis
    //   from a global fail (@decision level 0)
    //   */
    if (decisions.empty()) {
        return;
    }
    
    // number of literals of the current decision level in the conflict clause
    index_t lit_pointer{static_cast<index_t>(numLiteral())};
    Literal<T> l{Contradiction<T>};
    auto num_lit{0};
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
        std::cout << "analyze conflict (" << e << "):\n";
    }
#endif
    
    // the first reason is the conflicting clause itself
    Explanation<T> &exp = e;
    do {
        
        if (exp == Constant::NoReason<T>) {
            cut.add(l, lit_pointer);
        } else {
            
            lit_buffer.clear();
            exp.explain(l, lit_buffer);
            
#ifdef DBG_CLPLUS
            writeExplanation(l);
            checkClauseLimit();
#endif
            
            for (auto p : lit_buffer) {
                
#ifdef DBG_TRACE
                if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                    std::cout << "   ** " << pretty(p) << " @" << decisionLevel(p);
                    std::cout.flush();
                }
#endif
                
                // rank of p in the trail
                auto p_stamp{propagationStamp(p)};
                
                if (p_stamp < ground_stamp) {
                    // if p is a ground fact, we ignore it
#ifdef DBG_TRACE
                    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                        std::cout << " => ground fact\n";
                    }
#endif
                    continue;
                } else if (not cut.cached(p_stamp)) {
                    // p is not a ground fact nor a decision, and is not marked 'to
                    // explore' yet
                    cut.cache(p_stamp);
                    ++num_lit;
                    
#ifdef DBG_TRACE
                    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                        std::cout << " => add in 'to-explore'\n";
                    }
#endif
                }
            }
        }
        
#ifdef DBG_TRACE
        if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
            std::cout << num_lit << " open literals:";
            for (auto ii{lit_pointer}; --ii > 0;) {
                if (cut.cached(ii)) {
                    std::cout << " [" << getLiteral(ii) << "]";
                }
            }
            std::cout << std::endl;
        }
#endif
        
        if (num_lit == 0)
            break;
        
        while (not cut.cached(--lit_pointer))
            ;
        
        l = getLiteral(lit_pointer);
        exp = getReason(lit_pointer);
        
#ifdef DBG_TRACE
        if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
            std::cout << " explore " << l << " b/c " << exp << std::endl;
        }
#endif
        
        --num_lit;
        
    } while (true);
    
    cut.get(learnt_clause);
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
        std::cout << "\nlearn clause of size " << learnt_clause.size() << std::endl;
        
        for (auto l : learnt_clause) {
            std::cout << std::setw(4) << decisionLevel(~l) << " " << pretty(l)
            << std::endl;
        }
        std::cout << std::endl;
    }
#endif
    
#ifdef DBG_CL
    if (cl_file != NULL) {
        *cl_file << "0 " << (learnt_clause.size() + 1) << " 0 1 "
        << numeric.upper(1);
        for (auto p : learnt_clause) {
            writeLiteral(~p);
        }
        *cl_file << std::endl;
    }
    if (++num_clauses > DBG_CL) {
        std::cout << "exit because of dbg clause limit (#fails = " << num_fails
        << ", #cpts = " << num_choicepoints << ")\n";
        exit(1);
    }
#endif
    
    assert(static_cast<int>(decisionLevel(~(learnt_clause[0]))) == level());
}

template <typename T> bool Solver<T>::isAssertive(std::vector<Literal<T>> &conf) const {
    
    if(decisions.empty())
        return true;
    
    if (conf.empty()) {
        std::cout << "learnt empty clause @lvl " << level()
        << "! (#fails=" << num_fails << ")\n";
        return true;
    }
    
    if(conf[0].isNumeric() and numeric.falsified(conf[0])) {
        std::cout << "learned";
        for(auto l : conf)
            std::cout << " " << l;
        std::cout << " @" << level() << std::endl;
        
        std::cout << conf[0] << " is falsified! (#fails=" << num_fails << ")\n";
        return false;
    }
    if(conf[0].isNumeric() and numeric.satisfied(conf[0])) {
        std::cout << "learned";
        for(auto l : conf)
            std::cout << " " << l;
        std::cout << " @" << level() << std::endl;
        
        std::cout << conf[0] << " is satisfied! (#fails=" << num_fails << ")\n";
        return false;
    }
    if(not conf[0].isNumeric() and boolean.falsified(conf[0])) {
        std::cout << "learned";
        for(auto l : conf)
            std::cout << " " << l;
        std::cout << " @" << level() << std::endl;
        
        std::cout << conf[0] << " is falsified! (#fails=" << num_fails << ")\n";
        return false;
    }
    if(not conf[0].isNumeric() and boolean.satisfied(conf[0])) {
        std::cout << "learned";
        for(auto l : conf)
            std::cout << " " << l;
        std::cout << " @" << level() << std::endl;
        
        std::cout << conf[0] << " is satisfied! (#fails=" << num_fails << ")\n";
        return false;
    }
    
    if (conf.size() < 2)
        return true;
    
    for (size_t i{1}; i < conf.size(); ++i) {
        if(not conf[i].isNumeric() and not boolean.falsified(conf[i])) {
            
            std::cout << "learned";
            for(auto l : conf)
                std::cout << " " << l;
            std::cout << " @" << level() << std::endl;
            
            std::cout << conf[i] << " is not falsified! (#fails=" << num_fails
            << ")\n";
            return false;
        }
        if(conf[i].isNumeric() and not numeric.falsified(conf[i])) {
            
            std::cout << "learned";
            for(auto l : conf)
                std::cout << " " << l;
            std::cout << " @" << level() << std::endl;
            
            std::cout << conf[i] << " is not falsified! (#fails=" << num_fails
            << ")\n";
            return false;
        }
    }
    return true;
}

template <typename T> void Solver<T>::learnConflict(Explanation<T> &e) {
    
    //    std::cout << (int)(options.cut_type) << std::endl;
    
    if(options.cut_type == Options::Cut::Decisions) {
        analyzeDecisions(e);
    } else {
        analyze(e, options.cut_type == Options::Cut::Booleans);
    }
    
    assert(decisions.size() == static_cast<size_t>(level() - init_level));
    
    if (decisions.empty()) {
        throw SearchExhausted();
    }
    
    auto bt_level{std::max(init_level,cut.backtrackLevel(*this))};
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
        if (not learnt_clause.empty())
            std::cout << "learn clause of size " << learnt_clause.size() << " @lvl"
            << level() << ", backtrack to level " << bt_level
            << " and deduce " << pretty(learnt_clause[0]);
        else
            std::cout << "learn empty clause!";
    }
#endif
    
//    if (bt_level < init_level)
//        throw SearchExhausted();

    //    lit_buffer.clear();
    //    for (auto l : learnt_clause) {
    //        lit_buffer.push_back(~l);
    //    }
    //    for (auto i : cut.cached_) {
    //        lit_buffer.push_back(getLiteral(i));
    //    }
    ClauseAdded.trigger(*this);

    restoreState(bt_level);
        
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
        std::cout << " @lvl " << level() << std::endl;
    }
#endif
    
//    for (auto i{learnt_clause.begin()}; i != learnt_clause.end(); ++i) {
//        for (auto j{learnt_clause.begin()}; j != learnt_clause.end(); ++j) {
//            if (i != j) {
//                if (i->isNumeric() == j->isNumeric() and
//                    i->variable() == j->variable() and i->sign() == j->sign()) {
//                    std::cout << "duplicates!!! (" << num_fails << ")\n";
//                    
//                    std::cout << *i << " AND " << *j << std::endl;
//                    
//                    exit(1);
//                }
//            }
//        }
//    }
    
    assert(isAssertive(learnt_clause));

    clauses.add(learnt_clause.begin(), learnt_clause.end(), true, cut.glueScore(*this));

}

template <typename T> void Solver<T>::branchRight() {
    
//    std::cout << "branch right @" << level() << " (" << init_level << ")\n";
    
    assert(decisions.size() == static_cast<size_t>(level() - init_level));
    
    if (env.level() <= init_level)
        throw SearchExhausted();
    
    assert(not decisions.empty());
    
    auto deduction{~decisions.back()};
    DeductionMade.trigger(deduction);
    
    restoreState(env.level() - 1);
//    decisions.pop_back();
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
        std::cout << "-- backtrack to lvl " << env.level() << " & deduce "
        << deduction << ": " << pretty(deduction)
        << " [i=" << num_choicepoints << "] --\n";
        printTrace();
    }
#endif
    
    set(deduction);
}

template <typename T> void Solver<T>::initializeSearch() {
    
    start_time = cpu_time();
    stopWatch.start();
    
    if(not initialized) {
        
//        booleanActivityMap.resize(boolean.size(), heuristics::impl::ActivityMap::baseIncrement);
//        numericActivityMap.resize(numeric.size(), heuristics::impl::ActivityMap::baseIncrement);
        
        
//        std::cout << "here: " << numConstraint() << " / " << constraint_size << std::endl;
        
//        post(&clauses);
        
//        std::cout << "here: " << numConstraint() << " / " << constraint_size << std::endl;
        
        
        clauses.post(clauses.id());
        
        restartPolicy.initialize();
        if (not heuristic.isValid()) {
            heuristic = heuristics::make_heuristic(*this);
        }
        
        propag_pointer = 1;
        
        propagate();
        
        ground_stamp = assumption_stamp = numLiteral();
        
        initialized = true;
        if (options.verbosity >= Options::QUIET) {
            displayHeader(std::cout);
        }
    }
    
    
    //    std::cout << "INIT:\n";
    //    std::cout << boolean_size << " / " << boolean.size() << std::endl;
    //    std::cout << numeric_size << " / " << numeric.size() << std::endl;
    //    std::cout << constraint_size << " / " << constraints.size() << std::endl;
    
    assert(boolean.size() >= boolean_size);
    assert(numeric.size() >= numeric_size);
    assert(constraints.size() >= constraint_size);
    
//    boolean_size = boolean.size();
//    numeric_size = numeric.size();
//    constraint_size = constraints.size();
    
    //    heuristic.notifyStartSearch(*this);
}

template<typename T>
template<heuristics::heuristic<T> H>
void Solver<T>::setBranchingHeuristic(H &&h) {
    this->heuristic = std::forward<H>(h);
}

template <typename T> boolean_state Solver<T>::satisfiable() {
    
    auto satisfiability{UnknownState};
    try {
        initializeSearch();
    } catch(Failure<T>& f) {
        satisfiability = FalseState;
    }
    
    if(satisfiability == UnknownState) {
        satisfiability = search();
        if (satisfiability == TrueState) {
            saveSolution();
        }
    }
    if(options.verbosity >= Options::QUIET)
        displaySummary(
                       std::cout,
                       (satisfiability == TrueState
                        ? "sat "
                        : (satisfiability == FalseState ? "unsat " : "unknown ")));
    return satisfiability;
}

template <typename T> void Solver<T>::check_clause(const index_t i) {
    lit_buffer.clear();
    auto cl{clauses[i]};
    if (cl != NULL) {
        unsigned lsat{0};
        unsigned lunsat{-0};
        unsigned lundef{0};
        for (auto l : *cl) {
            if (l.isNumeric()) {
                if (numeric.satisfied(l)) {
                    ++lsat;
                } else if (numeric.falsified(l)) {
                    ++lunsat;
                } else {
                    ++lundef;
                    lit_buffer.push_back(l);
                }
            } else {
                if (boolean.satisfied(l)) {
                    ++lsat;
                } else if (boolean.falsified(l)) {
                    ++lunsat;
                } else {
                    ++lundef;
                    lit_buffer.push_back(l);
                }
            }
        }
        if (lsat == 0) {
            std::cout << "cl" << i << ": " << lundef << "/" << lunsat << "/" << lsat
            << ":";
            for (auto l : lit_buffer) {
                std::cout << " " << pretty(l);
            }
            std::cout << " watched: " << pretty(cl->watched(0)) << " & "
            << pretty(cl->watched(1)) << std::endl;
        } else {
            std::cout << "cl" << i << ": sat\n";
        }
    } else {
        std::cout << "cl" << i << ": does not exist\n";
    }
}

template <typename T> void Solver<T>::check_clauses(const char* msg) {
    auto cl{clauses.consistent()};
    if (cl != NULL) {
        std::cout << "\n" << msg << "\n@" << num_choicepoints << ": " << *cl << std::endl;
        
        std::cout << " watched:" << pretty(cl->watched(0)) << " & "
        << pretty(cl->watched(1)) << std::endl;
        
        for (auto p : *cl) {
            std::cout << " * " << pretty(p);
            if (p.isNumeric()) {
                if (numeric.satisfied(p))
                    std::cout << " (sat)\n";
                else if (numeric.falsified(p))
                    std::cout << " (unsat)\n";
                else
                    std::cout << " (undef)\n";
            } else {
                if (boolean.satisfied(p))
                    std::cout << " (sat)\n";
                else if (boolean.falsified(p))
                    std::cout << " (unsat)\n";
                else
                    std::cout << " (undef)\n";
            }
        }
        
        exit(1);
    }
}

template <typename T> void Solver<T>::minimize(const NumericVar<T> x) {
    MinimizationObjective<T> objective(x);
    optimize(objective);
}

template <typename T> void Solver<T>::maximize(const NumericVar<T> x) {
    MaximizationObjective<T> objective(x);
    optimize(objective);
}

template <typename T>
template <typename S>
void Solver<T>::optimize(S &objective) {
    objective.X.extract(*this);
    objective_var = objective.X.id();

    T lower_bound{0};

    try {
        initializeSearch();
        
        lower_bound = objective.getDual(*this);
        if (options.verbosity >= Options::NORMAL) {
            std::cout << "-- Lower bound = " << lower_bound << std::endl;
        }
        
    } catch (Failure<T> &f) {
        objective.setDual(objective.primalBound());
    }
    
    while (objective.gap() and not KillHandler::instance().signalReceived() and
           not searchCancelled and
           not(num_fails >= options.search_limit)) {
        auto satisfiability = search();
        if (satisfiability == TrueState) {
            
            
            
            auto best{objective.value(*this)};
            if (options.verbosity >= Options::NORMAL) {
                std::cout << std::setw(10) << best;
                displayProgress(std::cout);
            }
            
//            std::cout << "apply best\n";
            
//            if(best == 1328) {
////                std::cout << *this << std::endl;
//                debug_flag = true;
//            }
            
            objective.apply(best, *this);
            saveSolution();
            
//            std::cout << "restart\n";
            
            restart(true);
            //            assumption_stamp = numLiteral();
            try {
                objective.setPrimal(best, *this);
                propagate();
                
                auto lb{objective.getDual(*this)};
                if (options.verbosity >= Options::NORMAL and lb > lower_bound) {
                    lower_bound = lb;
                    std::cout << "-- Lower bound = " << lower_bound << std::endl;
                }
            } catch (Failure<T> &f) {
                satisfiability = FalseState;
            }
            
//            std::cout << "ok\n";
        }
        
        if (satisfiability == FalseState) {
            
//            std::cout << "set dual to primals\n";
            
            objective.setDual(objective.primalBound());
            
            if (options.verbosity >= Options::NORMAL) {
                lower_bound = objective.primalBound();
                std::cout << "-- Lower bound = " << lower_bound << std::endl;
            }
        }
    }
    
//    std::cout << "finished\n";
    
    //    std::cout << "hi " << objective.getDual(*this) << "\n";
    
    
    if (options.verbosity >= Options::QUIET) {
        auto msg{"optimal"};
        if(objective.gap() > 0) {
            if(KillHandler::instance().signalReceived() or searchCancelled)
                msg = "killed";
            else
                msg = "limit";
        }
        displaySummary(std::cout, msg);
        
        //        displaySummary(std::cout, (objective.gap() > 0 ? "killed" : "optimal"));
    }
}

template <typename T>
template <typename S, lns::relaxation_policy P>
void Solver<T>::largeNeighborhoodSearch(S &objective, P &&relaxationPolicy) {
    objective.X.extract(*this);
    objective_var = objective.X.id();
    initializeSearch();
    
    if (not boolean.hasSolution()) {
        auto satisfiability{search()};
        if (satisfiability == TrueState) {
            auto best{objective.value(*this)};
            if (options.verbosity >= Options::NORMAL) {
                std::cout << std::setw(10) << best;
                displayProgress(std::cout);
            }
            saveSolution();
            restart(true);
            try {
                objective.setPrimal(best, *this);
            } catch (Failure<T> &f) {
                objective.setDual(objective.primalBound());
            }
        } else {
            if (options.verbosity >= Options::QUIET)
                displaySummary(std::cout, "unsatisfiable");
            return;
        }
    }
    
    
    while (objective.gap() and not KillHandler::instance().signalReceived() and not searchCancelled) {
        lns::AssumptionProxy surrogate = *this;
        std::forward<P>(relaxationPolicy).relax(surrogate);
        
        // return to level 0 if there is no relaxation
        if (surrogate.getState() == lns::AssumptionState::Empty) {
            restoreState(0);
        }
        
        auto satisfiability = UnknownState;
        if (surrogate.getState() != lns::AssumptionState::Fail) {
            
            //            std::cout << "FIXED " << (assumption_stamp -
            //            ground_stamp) << " LITERALS\n";
            
            satisfiability = search();
        }
        
        if (satisfiability == TrueState) {
            auto best{objective.value(*this)};
            if (options.verbosity >= Options::NORMAL) {
                std::cout << std::setw(10) << best;
                displayProgress(std::cout);
            }
            saveSolution();
            std::forward<P>(relaxationPolicy).notifySuccess(num_fails);
            restoreState(0);
            assert(decisions.empty());
            
            try {
                objective.setPrimal(best, *this);
                propagate();
            } catch (Failure<T> &f) {
                objective.setDual(objective.primalBound());
            }
        } else {
            
            // no assumptions made and still failure => no improving solution exists
            if (surrogate.getState() == lns::AssumptionState::Empty and
                not KillHandler::instance().signalReceived()) {
                objective.setDual(objective.primalBound());
            } else {
                std::forward<P>(relaxationPolicy).notifyFailure(num_fails);
                restoreState(0);
                assert(decisions.empty());
            }
        }
    }
    
    if (options.verbosity >= Options::QUIET)
        displaySummary(std::cout, (objective.gap() > 0 ? "killed" : "optimal"));
}


template <typename T>
template <concepts::typed_range<Literal<T>> L>
void Solver<T>::makeAssumptions(const L &literals) {
    
    initializeSearch();
    saveState();
    
    for (auto lit : literals) {
        set(lit);
    }
    
    propagate();
    
    assumption_stamp = numLiteral();
    
    
//    std::cout << "make assumptions : " << ground_stamp << " " << assumption_stamp << std::endl;
    
}

template <typename T> void Solver<T>::makeAssumption(const Literal<T> lit) {
    
    initializeSearch();
    saveState();
    
    //  if(assumption_stamp == 0)
    //     assumption_stamp = numLiteral();
    
    set(lit);
    
    propagate();
    
    assumption_stamp = numLiteral();
    
//    std::cout << "make assumption : " << ground_stamp << " " << assumption_stamp << std::endl;
}

template <typename T> boolean_state Solver<T>::search() {
    
    //
    
    //    assumption_stamp =
    init_level = env.level();
    
    boolean_state satisfiability{UnknownState};
    while (satisfiability == UnknownState and not searchCancelled and
           not KillHandler::instance().signalReceived() and
           not(num_fails >= options.search_limit)) {
        
        try {
            
            PropagationInitiated.trigger(*this);
            propagate();
            PropagationCompleted.trigger(*this);
            
#ifdef DBG_TRACE
            if (DBG_BOUND) {
                std::cout << "--- propag [i=" << num_choicepoints << "] ---\n";
                printTrace();
            }
#endif
            
            // make a checkpoint
            saveState();
            
            // all resource constraints are accounted for => a solution has been found
            if (boolean_search_vars.empty() /* and numeric_search_vars.empty()*/) {
                satisfiability = TrueState;
                
#ifdef DBG_TRACE
                if (DBG_BOUND) {
                    std::cout << "--- new solution [i=" << num_choicepoints << "] ---\n";
                    //          printTrace();
                }
#endif
                
            } else {
                ++num_choicepoints;
                
                Literal<T> d = heuristic.branch(*this);
                decisions.push_back(d);
                ChoicePoint.trigger(d);
                
#ifdef DBG_TRACE
                if (DBG_BOUND) {
                    std::cout << "--- search node (lvl=" << env.level() << "/"
                    << init_level << ") [i=" << num_choicepoints
                    << "] (|trail|=" << numLiteral() << ")---\n";
                    std::cout << " ** take decision " << d << ": " << pretty(d)
                    << " **\n";
                    printTrace();
                }
#endif
                
                set(d);
            }
        } catch (Failure<T> &f) {
            
            avg_fail_level = (avg_fail_level * num_fails + env.level()) / (num_fails + 1);
            ++num_fails;
            
#ifdef DBG_FAIL
            f.reason.explain(Solver<T>::Contradiction, lit_buffer);
            for (auto l : lit_buffer) {
                std::cout << pretty(l) << std::endl;
            }
            lit_buffer.clear();
#endif
            
            try {
                backtrack(f.reason);
                BackTrackCompleted.trigger();
                if (restartPolicy.limit()) {
                    restart();
                }
            } catch (const SearchExhausted &f) {
                satisfiability = FalseState;
            }
            
#ifdef DBG_TRACE
            if (DBG_BOUND) {
                std::cout << "--- after fail [i=" << num_choicepoints << "] ---\n";
                printTrace();
            }
#endif
        }
    }
    
    return satisfiability;
}

template <typename T> void tempo::Solver<T>::propagate() {
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & PROPAGATION) and
        static_cast<index_t>(propag_pointer) < numLiteral()) {
        std::cout << "propagate\n";
        for (index_t i{propag_pointer}; i < numLiteral(); ++i)
            std::cout << " * " << pretty(trail[i]) << std::endl;
        std::cout << std::endl;
    }
#endif
    
    index_t p_index{static_cast<index_t>(propag_pointer)};
    
    clauses.clearTriggers();
    
    while (not propagation_queue.empty() or numLiteral() > p_index) {
        
        while (numLiteral() > p_index) {
            ++num_literals;
            //      Literal<T> l{trail[p_index]};
            Literal<T> l{getLiteral(p_index)};
            //      auto culprit{reason[p_index].expl};
            auto culprit{getReason(p_index).expl};
            
#ifdef DBG_TRACE
            if (DBG_BOUND and (DBG_TRACE & QUEUE)) {
                std::cout << "triggers for (" << l << ") b/c " << culprit->id() << "/"
                << getReason(p_index) << std::endl;
            }
#endif
            
            // TODO: not sure why it is better to do it like this than with the
            // standard constraint queue system (PRIORITY?)
            if (not l.isNumeric()) {
                
                clauses.unit_propagate_boolean(l);
                ++num_unit_propagations;
                
            } else {
                if (numeric.lastLitIndex(l.sign(), l.variable()) > p_index) {
                    ++p_index;
                    continue;
                    //                    std::cout << "subsumed\n";
                }
            }
            
            const std::vector<int> &cid =
            (l.isNumeric() ? numeric_constraints[l] : boolean_constraints[l]);
            const std::vector<unsigned> &rank =
            (l.isNumeric() ? numeric_constraints.rank(l)
             : boolean_constraints.rank(l));
            
            for (auto i{cid.size()}; i-- > 0;) {
                auto cons{constraints[cid[i]]};
                if (cons->idempotent and culprit == cons) {
                    continue;
                }
                
#ifdef DBG_TRACE
                if (DBG_BOUND and (DBG_TRACE & QUEUE)) {
                    std::cout << " -" << *cons << " (" << cons->id() << ")" << std::endl;
                }
#endif
                
                propagation_queue.triggers(l, rank[i], cons);
            }
            
            ++p_index;
        }
        
        if (not propagation_queue.empty()) {
            auto cons{propagation_queue.pop_front()};
            
#ifdef DBG_TRACE
            if (DBG_BOUND and (DBG_TRACE & QUEUE)) {
                std::cout << "propagate " << *cons << std::endl;
            }
#endif
            
            ++num_cons_propagations;
            cons->propagate();
        }
    }
    
    propag_pointer = p_index;
    if(options.ground_update)
        updateGroundTruth();
}

template <typename T> int Solver<T>::saveState() {
    assert(propag_pointer == static_cast<index_t>(numLiteral()));
    
    int lvl{env.level()};
    env.save();
    ReversibleObject::save();
    return lvl;
}

template <typename T> void Solver<T>::restoreState(const int l) {
    env.restore(l);
    
    if(boolean_size < boolean.size()) {
        
//        std::cout << "remove " << (boolean.size() - boolean_size) << " Boolean variables when backtracking to lvl " << l << std::endl;
        
        boolean.resize(boolean_size);
    }
    if(numeric_size < numeric.size()) {
        
//        std::cout << "remove " << (numeric.size() - numeric_size) << " numeric variables when backtracking to lvl " << l << std::endl;
        
        numeric.resize(numeric_size);
    }
    const size_t csize{constraint_size};
    
    assert(csize >= 1);
    
    
//    if(csize < constraints.size())
//        std::cout << "remove " << (constraints.size() - csize) << " constraints when backtracking to lvl " << l << std::endl;
    
    while(csize < constraints.size()) {
        delete constraints.back();
        constraints.pop_back();
    }
    //
    //    std::cout << "RESTORESTATE:\n";
    //    std::cout << boolean_size << " / " << boolean.size() << std::endl;
    //    std::cout << numeric_size << " / " << numeric.size() << std::endl;
    //    std::cout << constraint_size << " / " << constraints.size() << std::endl;
}

template <typename T> void Solver<T>::undo() {
    size_t n{propag_pointer};
    while (numLiteral() > n) {
        Literal<T> l{trail.back()};
        
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
                    << (bounds == bound::lower ? "from " : "to ") << v << "("
                    << (shortest_path[u] + w) << "/" << shortest_path[v] << ")"
                    << std::endl;
                }
#endif
                
                if (v == s) {
                    
#ifdef DBG_FAIL
                    if (DBG_FAIL) {
                        std::cout << " negative cyle\n";
                        
                        Explanation exp{&graph_exp,
                            static_cast<hint>(Literal<T>::index(bounds, s))};
                        exp.explain(Solver<T>::Contradiction, lit_buffer);
                        for (auto l : lit_buffer) {
                            std::cout << pretty(l) << std::endl;
                        }
                        lit_buffer.clear();
                    }
#endif
                    
                    throw Failure<T>(
                                     {&graph_exp, static_cast<hint>(Literal<T>::index(bounds, s))});
                }
                setNumeric(makeNumericLiteral<T>(bounds, v, shortest_path[u] + w),
                           {&graph_exp,
                    static_cast<hint>((edge.stamp() != Constant::NoIndex
                                       ? edge.stamp()
                                       : numeric.lastLitIndex(bounds, u)))},
                           false);
                
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

template <typename T> void Solver<T>::post(const DistanceConstraint<T> &c) {
    if (c.from == Constant::K) {
        // to - from <= dist -> to <= dist
        set(makeNumericLiteral(bound::upper, c.to, c.distance));
    } else if (c.to == Constant::K) {
        // to - from <= dist -> from >= -dist
        set(makeNumericLiteral(bound::lower, c.from, c.distance));
    } else {
        set(c);
    }
    propagate();
}

template <typename T> void Solver<T>::post(BooleanVar<T> con) {
    con.post(*this);
    propagate();
}

template <typename T> void Solver<T>::post(Constraint<T> *con) {
    
    constraints.push_back(con);
    propagation_queue.resize(constraints.size());
    
    boolean_constraints.resize(std::max(2 * boolean.size(), numConstraint()));
    numeric_constraints.resize(std::max(2 * numeric.size(), numConstraint()));
    
    con->post(numConstraint() - 1);
    
    propagate();
    
    constraint_size = numConstraint();
}

template <typename T> void Solver<T>::relax(Constraint<T> *con) {
    if (boolean_constraints.indegree(con->id()) > 0) {
        boolean_constraints.remove(con->id(), IN);
    }
    if (numeric_constraints.indegree(con->id()) > 0) {
        numeric_constraints.remove(con->id(), IN);
    }
}

template <typename T> void Solver<T>::addToSearch(const NumericVar<T> &x) {
    var_t var_id{x.id()};
    if (numeric.lower(var_id) < numeric.upper(var_id)) {
        numeric_search_vars.reserve(var_id + 1);
        if (not numeric_search_vars.has(var_id))
            numeric_search_vars.add(var_id);
    }
}

template <typename T> void Solver<T>::addToSearch(const BooleanVar<T> &x) {
    var_t var_id{x.id()};
    if (boolean.isUndefined(var_id)) {
        boolean_search_vars.reserve(var_id + 1);
        if (not boolean_search_vars.has(var_id))
            boolean_search_vars.add(var_id);
    }
}

template <typename T>
size_t Solver<T>::wake_me_on(const Literal<T> l, const int c) {
    if (l.isNumeric()) {
        if(numeric_constraints[l].empty() or numeric_constraints[l].back() != c)
            numeric_constraints.add(l, c);
        return numeric_constraints.rank(l).back();
    } else {
        if(boolean_constraints[l].empty() or boolean_constraints[l].back() != c)
            boolean_constraints.add(l, c);
        return boolean_constraints.rank(l).back();
    }
}

template <typename T>
template <typename ItLit>
void Solver<T>::postCardinality(const ItLit beg_lit, const ItLit end_lit,
                                const T bound) {
    post(new CardinalityConst<T>(*this, beg_lit, end_lit, bound));
}

template <typename T>
template <typename ItVar>
void Solver<T>::postCardinality(ItVar beg_var, ItVar end_var, const bool sign,
                                const T bound) {
    post(new CardinalityConst<T>(*this, beg_var, end_var, sign, bound));
}

template <typename T>
template <typename ItLit, typename ItW>
void Solver<T>::postPseudoBoolean(const ItLit beg_lit, const ItLit end_lit,
                                  ItW w, const T bound) {
    post(new PseudoBooleanConst<T>(*this, beg_lit, end_lit, w, bound));
}

template <typename T>
template <concepts::typed_range<Interval<T>> Tasks>
void Solver<T>::postEdgeFinding(Interval<T> &schedule, Tasks &&taskRange, Matrix<Literal<T>> lits) {
    post(new DisjunctiveEdgeFinding<T>(*this, schedule, std::forward<Tasks>(taskRange), std::move(lits)));
}

template <typename T>
template <concepts::typed_range<Interval<T>> Tasks>
void Solver<T>::postTransitivity(Interval<T> &schedule, Tasks &&taskRange, Matrix<Literal<T>> lits) {
    post(new Transitivity<T>(*this, schedule, std::forward<Tasks>(taskRange), std::move(lits)));
}

template <typename T>
template <typename ItRes>
FullTransitivity<T> *Solver<T>::postFullTransitivity(const ItRes beg_res,
                                                     const ItRes end_res) {
    auto c = new FullTransitivity<T>(*this);
    for (auto res{beg_res}; res != end_res; ++res) {
        using namespace std::views;
        auto disjunctive = res->getDisjunctiveLiterals().rawData() |
        filter([](auto lit) { return lit != Contradiction<T>; }) | common;
        c->addResource(disjunctive);
    }
    post(c);
    return c;
}

template <typename T>
template <concepts::typed_range<Interval<T>> Tasks, concepts::typed_range<NumericVar<T>> Demands>
void Solver<T>::postCumulative(const NumericVar<T> c, Tasks &&tasks, Demands &&demands, Matrix<Literal<T>> lits) {
    post(new CumulativeCheck<T>(*this, c, std::forward<Tasks>(tasks), std::forward<Demands>(demands), std::move(lits)));
}

template <typename T>
template <concepts::typed_range<Interval<T>> Tasks>
void Solver<T>::postCumulativeIncrementality(Tasks &&tasks,
                                             Matrix<Literal<T>> lits) {
    post(new Incrementality<T>(*this, std::forward<Tasks>(tasks),
                               std::move(lits)));
}

template <typename T>
template <typename ItTask, typename ItNVar>
void Solver<T>::postTimetabling(const NumericVar<T> c, const ItTask beg_task,
                                const ItTask end_task, const ItNVar beg_dem) {
    post(new CumulativeTimetabling<T>(*this, c, beg_task, end_task, beg_dem));
}

template <typename T>
template <typename ItTask, typename ItNVar>
void Solver<T>::postStrongEdgeFinding(const Interval<T> s,
                                      const NumericVar<T> c,
                                      const ItTask beg_task,
                                      const ItTask end_task,
                                      const ItNVar beg_dem, const bool tt,
                                      Incrementality<T> *b, const int approx) {
    post(new CumulativeEdgeFinding<T>(*this, s, c, beg_task, end_task, beg_dem,
                                      tt, b, approx));
}

template <typename T>
template <typename ItTask, typename ItNVar>
void Solver<T>::postOverlapFinding(
                                   const Interval<T> s, const NumericVar<T> c, const ItTask beg_task,
                                   const ItTask end_task, const ItNVar beg_dem, Matrix<Literal<T>> lits) {
    post(new CumulativeOverlapFinding<T>(*this, s, c, beg_task, end_task, beg_dem, lits));
}

template <typename T> double Solver<T>::looseness(const Literal<T> &l) const {
    if (l.isNumeric()) {
        
        if(numeric.falsified(l))
            return 0;
        
        auto lb{numeric.lower(l.variable())};
        auto ub{numeric.upper(l.variable())};
        
        if (l.sign() == bound::lower) {
            auto b{-l.value_unsafe()};
            double ls{static_cast<double>(ub - b)};
            assert(ls >= 0);
            return ls / static_cast<double>(ub - lb);
        } else {
            auto b{l.value_unsafe()};
            double ls{static_cast<double>(b - lb)};
            assert(ls >= 0);
            return ls / static_cast<double>(ub - lb);
        }
    } else if (l.hasSemantic()) {
        auto c{boolean.getEdge(l)};
        if(c != Constant::NoEdge<T>)
            return static_cast<double>(numeric.upper(c.from) - numeric.lower(c.to) +
                                       c.distance) /
            static_cast<double>(numeric.upper(c.to) - numeric.lower(c.from) -
                                c.distance);
    }
    return .5;
}

template <typename T>
std::ostream &Solver<T>::displayHeader(std::ostream &os,
                                       const int width) const {
    os << std::right << std::setw(width)
    << " objective   failures   branches    nds/s    lvl   clauses  size";
    os << "   cpu         wall-time\n";
    return os;
}

template <typename T>
std::ostream &Solver<T>::displaySummary(std::ostream &os,
                                        std::string msg) const {
    os << std::setw(10) << msg;
    displayProgress(os);
    return os;
}

template <typename T>
std::ostream &Solver<T>::displayProgress(std::ostream &os) const {
    
    auto cpu{(cpu_time() - start_time)};
    const auto wallTime = stopWatch.elapsed<std::chrono::milliseconds>();
    
    os << "  " << std::setw(9) << num_fails << "  " << std::setw(9)
    << num_choicepoints << "  " << std::setw(7)
    << static_cast<unsigned long>(static_cast<double>(num_choicepoints) / cpu)
    << "  " << std::setw(5) << static_cast<unsigned>(avg_fail_level) << "  "
    << std::setw(8) << clauses.size();
    if (clauses.size() == 0)
        os << "   n/a";
    else
        os << "  " << std::setw(4)
        << static_cast<unsigned>(static_cast<double>(clauses.volume()) /
                                 static_cast<double>(clauses.size()));
    
    
    os << "   " << std::left << std::setw(9) << cpu << "   " << wallTime << std::right << std::endl;
    
    return os;
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
    os << " " << boolean_search_vars.size() + boolean_search_vars.frontsize() << " Boolean search vars";
    if(boolean.size() > boolean_search_vars.capacity())
        os << " out of " << boolean.size() ;
    if(boolean_search_vars.frontsize() > 0)
        os << ", (" << boolean_search_vars.frontsize() << ") units:";
    for (auto b{boolean_search_vars.fbegin()}; b!=boolean_search_vars.fend(); ++b) {
        os << " " << pretty(boolean.getLiteral(boolean.isTrue(*b), *b)) ;
    }
    os << " " << numLiteral() << " literals" ;
    
    return os;
}

template <typename T>
std::ostream &Solver<T>::displayPrecedences(std::ostream &os) const {
    
    for (auto a : core) {
        for (auto e : core[a]) {
            int b{e};
            os << "x" << b;
            if (e.label() > 0) {
                os << " - " << e.label();
            } else if (e.label() < 0) {
                os << " + " << -e.label();
            }
            os << " <= x" << a << std::endl;
        }
    }
    
    return os;
}

template <typename T>
std::ostream &Solver<T>::displayTrail(std::ostream &os) const {
    
    size_t i{0};
    size_t j{1};
    while (j < numLiteral()) {
        auto l{getLiteral(j)};
        if (j == ground_stamp) {
            os << "g =[" << pretty(l) << "]";
        } else if (j == assumption_stamp) {
            os << "a =[" << pretty(l) << "]";
        } else if (i < decisions.size() and decisions[i] == l) {
            os << "d" << (i + 1) << "=[" << pretty(l) << "]";
            ++i;
        } else {
            os //<< (j > 1 ? ", " : " ")
            << pretty(l);
        }
        ++j;
        os << std::endl;
    }
    
    return os;
}

template <typename T>
std::ostream &Solver<T>::displayVariables(std::ostream &os) const {
    for (auto x : boolean_search_vars) {
        os << " b" << x;
        if (boolean.hasSemantic(x))
            os << ":[" << boolean.getEdge(true, x) << " or "
            << boolean.getEdge(false, x) << "]";
    }
    os << std::endl;
    return os;
}

template <typename T>
std::ostream &Solver<T>::displayConstraints(std::ostream &os) const {
    for (auto c : constraints) {
        os << *c << ": ";
        for (auto x : boolean_constraints.backward()[c->id()]) {
            os << (Literal<T>::sgn(x) == true ? " b" : " ¬b") << Literal<T>::var(x) ;
            //      os << " ¬b" << x;
        }
        os << " /";
        for (auto x : numeric_constraints.backward()[c->id()]) {
            os << (Literal<T>::sgn(x) == bound::lower ? " lb(x" : " ub(x") << Literal<T>::var(x) << ")";
        }
        os << std::endl;
    }
    return os;
}

//template <typename T>
//std::ostream &Solver<T>::displayPrecedences(std::ostream &os) const {
//  os << core;
//  return os;
//}

#ifdef LEARNING_RATE_STUFF
template <typename T> void Solver<T>::updateActivity(const Literal<T> l) {
    if (l.isNumeric()) {
        numeric.updateActivity(l.variable());
    } else {
        boolean.updateActivity(l.variable());
    }
}
#endif

template <typename T>
std::string Solver<T>::pretty(const Literal<T> l) const {
    std::stringstream ss;
    if(not l.isNumeric() and l.hasSemantic()) {
        ss << l << "(" << boolean.getEdge(l) << ")";
    } else {
        ss << l;
    }
    return ss.str();
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
        os << std::endl;
    }
    if (dom) {
        os << "domains:\n";
        displayDomains(os);
        os << std::endl;
    }
    if (sva) {
        os << "search vars:\n";
        displayVariables(os);
        os << std::endl;
    }
    if (pre) {
        os << "precedences:\n"; //<< core << std::endl;
        displayPrecedences(os);
        //      os << std::endl;
    }
    if (cla) {
        os << "clauses:\n" ;
        clauses.displayClauses(os);
        //      os << std::endl;
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
            os << l.lit << " b/c " << getReason(i++) << std::endl;
        }
    }
    if (con) {
        os << "constraints:\n";
        displayConstraints(os);
    }
    
    return os;
}

#ifdef DBG_TRACE
template <typename T> void Solver<T>::printTrace() const {
    if (DBG_TRACE & SEARCH) {
        display(std::cout, (DBG_TRACE & DOMAINS), (DBG_TRACE & BRANCH), false,
                false, (DBG_TRACE & CLAUSES), false, false, false,
                (DBG_TRACE & TRAIL));
        
        //    std::cout << "bsv: " << boolean_search_vars << std::endl;
    }
}
#endif

#ifdef DBG_CL
template <typename T> void Solver<T>::writeLiteral(const Literal<T> l) const {
    if (l.isNumeric()) {
        if (l.sign() == bound::lower) {
            *cl_file << " " << l.variable() << " 0 " << l.value();
        } else {
            *cl_file << " 0 " << l.variable() << " " << l.value();
        }
    } else {
        assert(l.hasSemantic());
        auto c{boolean.getEdge(l)};
        
        *cl_file << " " << c.from << " " << c.to << " " << c.distance;
    }
}


template <typename T> void Solver<T>::writeConflict() const {
    if (cl_file != NULL) {
        auto pb{options.ground_update and numeric.upper(1) != Constant::Infinity<T>};
        std::vector<Literal<T>> buffer;
        for (auto p : cut) {
            if(p.first != 0)
                buffer.push_back(p.second);
        }
        *cl_file << 2 << " " << (buffer.size() + pb);
        if(pb)
            *cl_file << " 0 1 " << numeric.upper(1);
        for (auto p : buffer) {
            writeLiteral(p);
        }
        *cl_file << std::endl;
    }
}

template <typename T> void Solver<T>::writeClause() const {
    if (cl_file != NULL) {
        
//        std::cout << "write clause " << learnt_clause.size() << " dec=" << decisions.size() << std::endl;
        
        auto pb{options.ground_update and numeric.upper(1) != Constant::Infinity<T>};
        *cl_file << 2 << " " << (learnt_clause.size() + pb);
        if(pb)
            *cl_file << " 0 1 " << numeric.upper(1);
        for (auto p : learnt_clause) {
            writeLiteral(~p);
        }
        *cl_file << std::endl;
    }
}

template <typename T> void Solver<T>::writeExplanation(const Literal<T> l) const {
    if (cl_file != NULL) {
        *cl_file << "1 " << (lit_buffer.size() + (l != Contradiction<T>));
        for (auto p : lit_buffer) {
            writeLiteral(p);
        }
        if (l != Contradiction<T>)
            writeLiteral(~l);
        *cl_file << std::endl;
    }
}


template <typename T> void Solver<T>::checkClauseLimit() {
    if (++num_clauses > DBG_CL) {
        std::cout << "exit because of dbg clause limit (#fails = " << num_fails
        << ", #cpts = " << num_choicepoints << ")\n";
        exit(1);
    }
    
//    std::cout << "learn clause #" <<num_clauses << std::endl;
}



#endif

template <typename T>
std::ostream &operator<<(std::ostream &os, const Solver<T> &x) {
  return x.display(os);
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const ConflictSet<T> &x) {
  return x.display(os);
}
}

#endif

