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
#include "constraints/CumulativeCheck.hpp"
#include "constraints/CumulativeEdgeFinding.hpp"
#include "constraints/DisjunctiveEdgeFinding.hpp"
#include "constraints/EdgeConstraint.hpp"
//#include "constraints/OptionalEdgeConstraint.hpp"
#include "constraints/PseudoBoolean.hpp"
#include "constraints/Transitivity.hpp"
#include "constraints/FullTransitivity.hpp"
#include "heuristics/heuristic_factories.hpp"
#include "heuristics/impl/DecayingEventActivityMap.hpp"
#include "RelaxationPolicy.hpp"
#include "util/KillHandler.hpp"
#include "util/Options.hpp"
#include "util/SubscribableEvent.hpp"


namespace tempo {

//! T is the numeric variable domain type
template<typename T> class Solver;

template <typename T> class Solution {
    
public:
    Solution();
    Solution(const Solver<T> &);
    void load(const Solver<T> &);
    
    bool value(const BooleanVar<T> x) const;
    bool value(const var_t x) const;
    
    T lower_bound(const NumericVar<T> x) const;
    T lower_bound(const var_t x) const;
    
    T upper_bound(const NumericVar<T> x) const;
    T upper_bound(const var_t x) const;
    
private:
    std::vector<bool> boolean_variable;
    std::vector<T> numeric_lower_bound;
    std::vector<T> numeric_upper_bound;
};

//! Boolean variables and literals manager
/*!
 Responsible for:
 Storage (memory)
 Read/write access
 */
template<typename T>
class BooleanStore {
    
public:
    /**
     * @name constructors
     */
    //@{
    BooleanStore(Solver<T> &s);
    ~BooleanStore() = default;
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
    //  DisjunctVar<T> newDisjunct(const DistanceConstraint<T> &d1,
    //                             const DistanceConstraint<T> &d2);
    BooleanVar<T> newDisjunct(const DistanceConstraint<T> &d1,
                              const DistanceConstraint<T> &d2);
    //@}
    
    /**
     * @name utils
     */
    //@{
    // number of Boolean variables
    size_t size() const;
    
    // set literal l to true
    void set(Literal<T> l);
    
    // unset literal l
    void undo(Literal<T> l);
    
    // returns a literal from a sign and a variable
    Literal<T> getLiteral(const bool s, const var_t x) const;
    
    // returns the rank of literal l in the trail
    index_t litIndex(const Literal<T> l) const;
    
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
    
    // saves the current solution
    void saveSolution() { best_solution = polarity; }
    bool hasSolution() const { return not best_solution.empty(); }
    auto bestSolution() const noexcept -> const std::vector<bool> & { return best_solution; }
    //@}
    
protected:
    Solver<T> &solver;
    
    // [for each literal] the current polarity (x is undefined if x and ~x are
    // both false)
    std::vector<bool> polarity;
    
    // [for each literal] the polarity in the best solution
    std::vector<bool> best_solution;
    
    // the rank of the difference logic constraint in "edges" for each Boolean
    // variable
    std::vector<info_t> edge_index;
    
    // the list of difference logic constraints
    std::vector<DistanceConstraint<T>> edges;
    
    // the rank of each literal in the trail (Constant::NoIndex if the literal
    // is not on the trail)
    std::vector<index_t> propagation_level;
};

//! Numeric variables and literals manager
/*!
 Responsible for:
 Storage (memory)
 Read/write access
 */
template<typename T>
class NumericStore {
    
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
    T upper(const NumericVar<T> x) const;
    T lower(const NumericVar<T> x) const;
    
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
    
    size_t size() const;
    
    void set(Literal<T> l);
    void undo(Literal<T> l);
    
    const std::vector<T> &get(const int b) const;
    //@}
    
    /**
     * @name learning helpers
     */
    //@{
    index_t getConflictIndex(const Literal<T> l) const;
    void setConflictIndex(const Literal<T> l, T v);
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
        best_solution[bound::lower] = bound[bound::lower];
        best_solution[bound::upper] = bound[bound::upper];
    }
    bool hasSolution() const { return not best_solution[bound::lower].empty(); }
    std::vector<T> &bestSolution(const int b) const { return best_solution[b]; }
    //@}
    
private:
    Solver<T> &solver;
    
    // [for each numeric signed_var] the current bounds
    std::vector<T> bound[2];
    
    // [for each numeric signed_var] the bound in the best solution
    std::vector<T> best_solution[2];
    
    // [for each numeric signed_var] the current index in the 'propagation_events'
    // stack
    std::vector<std::vector<index_t>> bound_index[2];
    
    // used for learning
    std::vector<index_t> conflict_index[2];
};

//! Explainer for literals from difference logic
template <typename T = int> class GraphExplainer : public Explainer<T> {
    
public:
    GraphExplainer(Solver<T> &s);
    void xplain(const Literal<T>, const hint, std::vector<Literal<T>> &) override;
    std::ostream &print_reason(std::ostream &os, const hint) const override;
    //    int getType() const override;
    
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
    mutable SubscribableEvent<const std::vector<Literal<T>> &>
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
    //    // create an internal Boolean variable with a difference logic
    //    semantic,
    //    // post the channelling constraints, and return a model object
    //    pointing to
    //    // it
    //    DisjunctVar<T> newDisjunct(const DistanceConstraint<T> &,
    //                               const DistanceConstraint<T> &);
    //    // create an internal Boolean variable with a difference logic
    //    semantic,
    //    // post the channelling constraints (including the optionality), and
    //    return
    //    // a model object pointing to it
    //    DisjunctVar<T> newDisjunct(const DistanceConstraint<T> &d1,
    //                               const DistanceConstraint<T> &d2,
    //                               const var_t opt);
    // create an internal Boolean variable with a difference logic semantic,
    // post the channelling constraints, and return a model object pointing to
    // it
    BooleanVar<T> newDisjunct(const DistanceConstraint<T> &,
                              const DistanceConstraint<T> &);
    //    // create an internal Boolean variable with a difference logic semantic,
    //    // post the channelling constraints (including the optionality), and return
    //    // a model object pointing to it
    //    BooleanVar<T> newDisjunct(const DistanceConstraint<T> &d1,
    //                              const DistanceConstraint<T> &d2, const var_t opt);
    // returns the constant 0
    NumericVar<T> zero() { return NumericVar<T>(Constant::K,0); }
    // returns the constant true
    BooleanVar<T> truism() { return BooleanVar<T>(0); }
    // create an internal numeric variable and return a model object pointing to it
    NumericVar<T> newNumeric(const T lb = -Constant::Infinity<T>,
                             const T ub = Constant::Infinity<T>);
    
    NumericVar<T> newConstant(const T k);
    NumericVar<T> newOffset(NumericVar<T> &x, const T k);
    // create an internal temporal variable and return a model object pointing
    // to it
    //    TemporalVar<T> newTemporal(const T offset = 0);
    //    NumericVar<T> newTemporal(const T offset = 0);
    // create the internal variables (depending on the type of Interval) and
    // return a model object pointing to them
    Interval<T> newInterval(const T mindur = 0,
                            const T maxdur = Constant::Infinity<T>,
                            const T earliest_start = -Constant::Infinity<T>,
                            const T latest_start = Constant::Infinity<T>,
                            const T earliest_end = -Constant::Infinity<T>,
                            const T latest_end = Constant::Infinity<T>,
                            const BooleanVar<T> opt = Constant::True
                            );
    
    //    Interval<T> between(const NumericVar<T> s, const NumericVar<T> e, const bool opt=false);
    //    Interval<T> continuefor(const NumericVar<T> s, const NumericVar<T> d, const bool opt=false);
    //
    
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
    
    // get the bound literal that precedes l
    //    Literal<T> previousBound(const Literal<T> l) const;
    //    bool subsumed(const Literal<T> l) const;
    //    Literal<T> previousLiteral(const Literal<T> l) const;
    
    //    // get the most recent Literal that entails l
    //    Literal<T> getImplicant(const Literal<T> l) const;
    
    // get the explanation for the i-th literal
    Explanation<T> getReason(const index_t i) const;
    
    // get the explanation for the i-th literal
    Explanation<T> getReason(const Literal<T> l) const;
    
    // get the index in the propagation queue of the last Literal involving
    // variable x (to be used parcimoniously, not so efficient)
    index_t propagationLevel(const Literal<T> l) const;
    
    // for debugging purpose only, very inefficient
    index_t decisionLevel(const Literal<T> l) const;
    
    // whether the current conflict implies this literal (its level is passed as argument for efficiency purpose)
    bool entailedByConflict(Literal<T>, const index_t) const;
    
    // set literal l true with explanation e
    void set(Literal<T> l, const Explanation<T> &e = Constant::Decision<T>);
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
    template <typename ItTask, typename ItVar>
    void postEdgeFinding(Interval<T> &schedule, const ItTask beg_task,
                         const ItTask end_task, const ItVar beg_var);
    
    // create and post a new precedence reasoning propagator
    template <typename ItTask, typename ItVar>
    void postTransitivity(Interval<T> &schedule, const ItTask beg_task,
                          const ItTask end_task, const ItVar beg_var);
    
    // create and post a new full transitivity propagator
    template <typename ItRes>
    FullTransitivity<T> *postFullTransitivity(const ItRes beg_res,
                                              const ItRes end_res);
    
    // create and post a new precedence reasoning propagator
    template <typename ItTask, typename ItNVar, typename ItBVar>
    void postCumulative(const NumericVar<T> c, const ItTask beg_task,
                        const ItTask end_task, const ItNVar beg_dem,
                        const ItBVar beg_disj);
    template <typename ItTask, typename ItNVar, typename ItBVar>
    void postStrongEdgeFinding(const Interval<T> s, const NumericVar<T> c,
                               const ItTask beg_task, const ItTask end_task,
                               const ItNVar beg_dem, const ItBVar beg_disj);
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
    void setHeuristic(heuristics::MovableHeuristic<heuristics::PolymorphicHeuristic<T>> h);
    
    template<heuristics::heuristic<T> H>
    void setBranchingHeuristic(H &&h);
    // record the backtrack-environment level when calling 'initializeSearch()'
    int init_level{0};
    // init_level + 1 if there are assumptions, init_level otherwise
    //    int assumption_level{0};
    // index of the first decision in the trail
    index_t ground_level{0};
    index_t assumption_level{0};
    // apply the assumptions
    template <typename IterLit>
    void makeAssumptions(IterLit beg_a, IterLit end_a);
    //@}
    
    /**
     * @name search
     */
    //@{
    boolean_state search();
    
    template <typename S> void optimize(S &objective);
    
    template <typename S, typename A>
    void largeNeighborhoodSearch(S &objective, A &relaxationPolicy);
    
    boolean_state satisfiable();
    void minimize(const NumericVar<T> x);
    void maximize(const NumericVar<T> x);
    
    void restart(const bool on_solution = false);
    void backtrack(Explanation<T> &e);
    void branchRight();
    void learnConflict(Explanation<T> &e);
    void analyze(Explanation<T> &e);
    void minimize_clause();
    
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
    
private:
    // reversible strutures
    BacktrackEnvironment env;
    
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
    
    //    // reversible strutures
    //    BacktrackEnvironment env;
    
    // decision stack
    std::vector<Literal<T>> decisions;
    
    
    // the stack of Literals reprensenting all the changes so far
    std::vector<Literal<T>> trail;
    
    // the reason for each propagation event
    std::vector<Explanation<T>> reason;
    
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
    
    /**
     * @name helpers for conflict-analysis
     */
    //@{
    std::vector<bool> explored;
    
    std::vector<Literal<T>> lit_buffer;
    std::vector<Literal<T>> conflict;
    std::vector<index_t> literal_lvl;
    std::vector<Literal<T>> learnt_clause;
    std::vector<Literal<T>> minimal_clause;
    
public:
    std::vector<Literal<T>>::iterator begin_learnt() {
        return learnt_clause.begin();
    }
    std::vector<Literal<T>>::iterator end_learnt() {
        return learnt_clause.end();
    }
    // @}
    
private:
    // specialisation to numeric literals
    void setNumeric(Literal<T> l,
                    const Explanation<T> &e = Constant::Decision<T>,
                    const bool do_update = true);
    
    // specialisation to Boolean literals
    void setBoolean(Literal<T> l,
                    const Explanation<T> &e = Constant::Decision<T>);
    
    
#ifdef DBG_TRACE
    void printTrace() const;
#endif
    
    heuristics::impl::EventActivityMap *activityMap{nullptr};
    
    // to restack assumptions after initial propagation
    //    SubscriberHandle assumptionHandler;
    
public:
    void setActivityMap(heuristics::impl::EventActivityMap *map) {
        activityMap = map;
    }
    heuristics::impl::EventActivityMap *getActivityMap() {
        return activityMap;
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
    // average depth of the search tree
    double avg_fail_level{0};
    
    //@}
    
    // Global constant for failures
    static constexpr Literal<T> Contradiction =
    makeBooleanLiteral<T>(false, Constant::NoVar, 0);
    
    /**
     * @name to collect the modeling construct in order to free memory up
     */
    //@{
    std::vector<ExpressionImpl<T>*> trash_bin;
    //@}
    
private:
    /**
     * @name debug
     */
    //@{
    bool isAssertive(std::vector<Literal<T>> &conf) const;
    
#ifdef DBG_CL
    std::ofstream *cl_file{NULL};
    int num_clauses{0};
    
    void writeLiteral(const Literal<T> l) const;
#endif
    
    //@}
};

template <typename T> Solution<T>::Solution() {}

template <typename T> Solution<T>::Solution(const Solver<T> &solver) {
    load(solver);
}

template <typename T> void Solution<T>::load(const Solver<T> &solver) {
    boolean_variable = solver.boolean.bestSolution();
    numeric_lower_bound = solver.numeric.bestSolution(bound::lower);
    numeric_upper_bound = solver.numeric.bestSolution(bound::upper);
}

template <typename T> bool Solution<T>::value(const BooleanVar<T> x) const {
    return boolean_variable[x.id()];
}

template <typename T> bool Solution<T>::value(const var_t x) const {
    return boolean_variable[x];
}

template <typename T> T Solution<T>::lower_bound(const NumericVar<T> x) const {
    return numeric_lower_bound[x.id()];
}

template <typename T> T Solution<T>::lower_bound(const var_t x) const {
    return numeric_lower_bound[x];
}

template <typename T> T Solution<T>::upper_bound(const NumericVar<T> x) const {
    return numeric_upper_bound[x.id()];
}

template <typename T> T Solution<T>::upper_bound(const var_t x) const {
    return numeric_upper_bound[x];
}

/*!
 GraphExplainer implementation
 */
template <typename T>
void GraphExplainer<T>::xplain(const Literal<T> l, const hint h,
                               std::vector<Literal<T>> &Cl) {
    
    if (l == Solver<T>::Contradiction) {
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

//template <typename T> int GraphExplainer<T>::getType() const {
//    return CYCLEEXPL;
//}

template <typename T>
GraphExplainer<T>::GraphExplainer(Solver<T> &s) : solver(s) {}

template <typename T>
void BoundExplainer<T>::xplain(const Literal<T> l, const hint h,
                               std::vector<Literal<T>> &Cl) {
    
    if (l == Solver<T>::Contradiction) {
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
            std::cout << solver.getLiteral(lidx) << " AND " << solver.getLiteral(uidx)
            << std::endl;
        }
#endif
        
        Literal<T> le;
        Explanation<T> exp;
        
        if (lidx < uidx) {
            le = solver.getLiteral(lidx);
            exp = solver.getReason(uidx);
        } else {
            le = solver.getLiteral(uidx);
            exp = solver.getReason(lidx);
        }
        
#ifdef DBG_TRACE
        if (DBG_CBOUND and (DBG_TRACE & LEARNING)) {
            std::cout << " => " << le << " AND " << exp << std::endl;
        }
#endif
        
        Cl.push_back(le);
        exp.explain(~le, Cl);
        
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

//template <typename T> int BoundExplainer<T>::getType() const {
//    return BOUNDEXPL;
//}

template <typename T>
BoundExplainer<T>::BoundExplainer(Solver<T> &s) : solver(s) {}

/*!
 BooleanStore implementation
 */
template <typename T>
Literal<T> BooleanStore<T>::getLiteral(const bool s, const var_t x) const {
    return makeBooleanLiteral<T>(s, x, edge_index[x]);
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

template <typename T>
const DistanceConstraint<T> &BooleanStore<T>::getEdge(const index_t i) const {
    return edges[i];
}

template <typename T> bool BooleanStore<T>::hasSemantic(const var_t x) const {
    return edge_index[x] >= Constant::SomeSemantic;
}

template <typename T> BooleanStore<T>::BooleanStore(Solver<T> &s) : solver(s) {
    // I don't remember what's that for
    edges.push_back(Constant::NoEdge<T>);
    edges.push_back(Constant::NoEdge<T>);
    
    // this is to have the constants true/false
    auto x{newVar()};
    set(makeBooleanLiteral<T>(true, x.id()));
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

// template <typename T>
// DisjunctVar<T> BooleanStore<T>::newDisjunct(const DistanceConstraint<T> &d1,
//                                             const DistanceConstraint<T> &d2)
//                                             {
//     info_t d{static_cast<info_t>(edges.size())};
//     DisjunctVar<T> x{newVar(d).id(), d};
//     edges.push_back(d1);
//     edges.push_back(d2);
//     return x;
// }

template <typename T>
BooleanVar<T> BooleanStore<T>::newDisjunct(const DistanceConstraint<T> &d1,
                                           const DistanceConstraint<T> &d2) {
    info_t d{static_cast<info_t>(edges.size())};
    BooleanVar<T> x{newVar(d).id(), d};
    edges.push_back(d1);
    edges.push_back(d2);
    return x;
}

template <typename T> void BooleanStore<T>::set(Literal<T> l) {
    propagation_level[l.variable()] = (solver.numLiteral() - 1);
    polarity[l] = true;
    if (l.hasSemantic()) {
        assert(l.constraint() == (edge_index[l.variable()] + l.sign()));
        auto e{edges[l.constraint()]};
        //        if (e != DistanceConstraint<T>::none)
        if (e != Constant::NoEdge<T>)
            solver.set(e, static_cast<index_t>(solver.numLiteral() - 1));
    }
    assert(l.hasSemantic() or edge_index[l.variable()] == Constant::NoSemantic);
}

template <typename T> void BooleanStore<T>::undo(Literal<T> l) {
    polarity[l] = false;
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

// template <typename T> bool BooleanStore<T>::visited(const Literal<T> l) const
// {
//   return explored[l.variable()];
// }
//
// template <typename T>
// void BooleanStore<T>::markVisited(const Literal<T> l)  {
//   explored[l.variable()] = true;
// }
//
// template <typename T>
// void BooleanStore<T>::unmarkVisited(const Literal<T> l)  {
//   explored[l.variable()] = false;
// }

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

/*!
 NumericStore implementation
 */
template <typename T> NumericStore<T>::NumericStore(Solver<T> &s) : solver(s) {
    // useful to have the constant 0 as a regular variable
    // other constants can be TemporalVars pointing to the constant 0
    //  newVar(0,0);
    newVar(0);
}

template <typename T> size_t NumericStore<T>::size() const {
    return bound[bound::lower].size();
}

template <typename T>
const std::vector<T> &NumericStore<T>::get(const int b) const {
    return bound[b];
}

// template <typename T> NumericVar<T> NumericStore<T>::newVar(const T lb, const
// T ub) {
//   NumericVar<T> x{static_cast<var_t>(size())};
//
//   bound[bound::lower].push_back(lb);
//   bound[bound::upper].push_back(ub);
//
//   conflict_index[bound::lower].push_back(Constant::NoIndex);
//   conflict_index[bound::upper].push_back(Constant::NoIndex);
//
//   bound_index[bound::lower].resize(size());
//   bound_index[bound::upper].resize(size());
//   bound_index[bound::lower].back().push_back(Constant::InfIndex);
//   bound_index[bound::upper].back().push_back(Constant::InfIndex);
//
//   return x;
// }

template <typename T> NumericVar<T> NumericStore<T>::newVar(const T b) {
    NumericVar<T> x{static_cast<var_t>(size())};
    
    bound[bound::lower].push_back(b);
    bound[bound::upper].push_back(b);
    
    conflict_index[bound::lower].push_back(Constant::NoIndex);
    conflict_index[bound::upper].push_back(Constant::NoIndex);
    
    bound_index[bound::lower].resize(size());
    bound_index[bound::upper].resize(size());
    bound_index[bound::lower].back().push_back(Constant::InfIndex);
    bound_index[bound::upper].back().push_back(Constant::InfIndex);
    
    return x;
}

template <typename T> void NumericStore<T>::set(Literal<T> l) {
    auto s{l.sign()};
    auto v{l.variable()};
    
    assert(bound[s][v] > l.value());
    
    bound[s][v] = l.value();
    bound_index[s][v].push_back(static_cast<index_t>(solver.numLiteral() - 1));
}

template <typename T> void NumericStore<T>::undo(Literal<T> l) {
    auto s{l.sign()};
    auto v{l.variable()};
    bound_index[s][v].pop_back();
    bound[s][v] = solver.getLiteral(bound_index[s][v].back()).value();
}

// template <typename T> bool NumericStore<T>::visited(const Literal<T> l) const
// {
//   return explored_bound[l.sign()][l.variable()] >= l.value();
// }
//
// template <typename T>
// void NumericStore<T>::markVisited(const Literal<T> l)  {
//   explored_bound[l.sign()][l.variable()] = l.value();
// }
//
// template <typename T>
// void NumericStore<T>::unmarkVisited(const Literal<T> l)  {
//   explored_bound[l.sign()][l.variable()] = -Constant::Infinity<T>;
// }

template <typename T>
index_t NumericStore<T>::getConflictIndex(const Literal<T> l) const {
    return conflict_index[l.sign()][l.variable()];
}

template <typename T>
void NumericStore<T>::setConflictIndex(const Literal<T> l, T v) {
    conflict_index[l.sign()][l.variable()] = v;
}

template <typename T> T NumericStore<T>::upper(const NumericVar<T> x) const {
    
    assert(best_solution[bound::upper].size() == bound[bound::upper].size());
    
    assert(best_solution[bound::upper].size() > x.id());
    
    return best_solution[bound::upper][x.id()] + x.offset();
}
template <typename T> T NumericStore<T>::lower(const NumericVar<T> x) const {
    
    assert(best_solution[bound::lower].size() == bound[bound::lower].size());
    
    assert(best_solution[bound::lower].size() > x.id());
    
    return -best_solution[bound::lower][x.id()] + x.offset();
}

template <typename T> T NumericStore<T>::upper(const var_t x) const {
    return bound[bound::upper][x];
}

template <typename T> T NumericStore<T>::lower(const var_t x) const {
    return -bound[bound::lower][x];
}

template <typename T>
Literal<T> NumericStore<T>::getLiteral(const bool s, const var_t x) const {
    return solver.getLiteral(bound_index[s][x].back());
}

template <typename T>
Literal<T> NumericStore<T>::previousBound(const Literal<T> l) const {
    //    auto i{numeric.bound_index[l.sign()][l.variable()]};
    
    auto i{bound_index[l.sign()][l.variable()].rbegin()};
    
    assert(solver.getLiteral(*i) == l);
    
    //    std::cout << "prev bound of " << l << " (top of stack = " <<
    //    solver.getLiteral(*i) << ")\n";
    //
    //    while (solver.getLiteral(*i) != l) {
    //        std::cout << " -- " <<  solver.getLiteral(*i) << "? (" <<
    //        std::distance(i, bound_index[l.sign()][l.variable()].rend()) <<
    //        ")\n";
    //        --i;
    //    }
    
    --i;
    
    return solver.getLiteral(*i);
}

template <typename T>
index_t NumericStore<T>::lastLitIndex(const bool s, const var_t x) const {
    return bound_index[s][x].back();
}

template <typename T>
index_t NumericStore<T>::litIndex(const Literal<T> l) const {
    auto i{bound_index[l.sign()][l.variable()].rbegin()};
    while (solver.getLiteral(*(i + 1)).value() <= l.value())
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

/*!
 Solver implementation
 */
template <typename T>
Solver<T>::Solver()
: ReversibleObject(&env), boolean(*this), numeric(*this), clauses(*this),
core(&env), boolean_search_vars(0, &env), numeric_search_vars(0, &env),
propag_pointer(1, &env), propagation_queue(constraints),
boolean_constraints(&env), numeric_constraints(&env),
restartPolicy(*this), graph_exp(*this), bound_exp(*this)
//, assumptionHandler(PropagationCompleted.subscribe_handled([this](const auto
//&) { this->makeAssumptions();}))
{
    trail.emplace_back(Constant::NoVar, Constant::Infinity<T>, detail::Numeric{});
    reason.push_back(Constant::GroundFact<T>);
    core.newVertex(0);
    seed(options.seed);
}

template <typename T> Solver<T>::~Solver() {
    for (auto exp : trash_bin) {
        delete exp;
    }
    
    for (auto c : constraints) {
        if (c != &clauses) {
            //          std::cout << "delete " << *c << std::endl;
            delete c;
        }
    }
}

/*!
 Solver implementation
 */
template <typename T>
Solver<T>::Solver(Options opt)
: ReversibleObject(&env), boolean(*this), numeric(*this), clauses(*this),
core(&env), boolean_search_vars(0, &env), numeric_search_vars(0, &env),
options(std::move(opt)), propag_pointer(1, &env),
propagation_queue(constraints), boolean_constraints(&env),
numeric_constraints(&env), restartPolicy(*this), graph_exp(*this),
bound_exp(*this)
//, assumptionHandler(PropagationCompleted.subscribe_handled([this](const auto
//&) { this->makeAssumptions(); }))
{
    trail.emplace_back(Constant::NoVar, Constant::Infinity<T>, detail::Numeric{});
    reason.push_back(Constant::Decision<T>);
    core.newVertex(0);
    seed(options.seed);
    
#ifdef DBG_CL
    if (options.dbg_file != "")
        cl_file = new std::ofstream(options.dbg_file, std::ofstream::out);
#endif
}

template <typename T> BooleanVar<T> Solver<T>::newBoolean() {
    auto x{boolean.newVar()};
    clauses.newBooleanVar(x.id());
    boolean_constraints.resize(std::max(numConstraint(), 2 * boolean.size()));
    return x;
}

// template <typename T>
// DisjunctVar<T> Solver<T>::newImplication(const DistanceConstraint<T> &d) {
//     auto x{boolean.newDisjunct(d, DistanceConstraint<T>::none)};
//     clauses.newBooleanVar(x.id());
//     boolean_constraints.resize(std::max(numConstraint(), 2 *
//     boolean.size()));
//
//     post(new EdgeConstraint<T>(*this, boolean.getLiteral(true, x.id())));
//
//     return x;
// }

// template <typename T>
// DisjunctVar<T> Solver<T>::newDisjunct(const DistanceConstraint<T> &d1,
//                                       const DistanceConstraint<T> &d2) {
//     auto x{boolean.newDisjunct(d1, d2)};
//     clauses.newBooleanVar(x.id());
//     boolean_constraints.resize(std::max(numConstraint(), 2 *
//     boolean.size()));
//
//     if (d1 != DistanceConstraint<T>::none)
//       post(new EdgeConstraint<T>(*this, boolean.getLiteral(true, x.id())));
//
//     if (d2 != DistanceConstraint<T>::none)
//       post(new EdgeConstraint<T>(*this, boolean.getLiteral(false, x.id())));
//
//     return x;
// }
//
// template <typename T>
// DisjunctVar<T> Solver<T>::newDisjunct(const DistanceConstraint<T> &d1,
//                                       const DistanceConstraint<T> &d2,
//                                       const var_t opt) {
//   auto x{boolean.newDisjunct(d1, d2)};
//   clauses.newBooleanVar(x.id());
//   boolean_constraints.resize(std::max(numConstraint(), 2 * boolean.size()));
//
//   post(new OptionalEdgeConstraint<T>(*this, x.id(), opt));
//
//   return x;
// }

template <typename T>
BooleanVar<T> Solver<T>::newDisjunct(const DistanceConstraint<T> &d1,
                                     const DistanceConstraint<T> &d2) {
    auto x{boolean.newDisjunct(d1, d2)};
    clauses.newBooleanVar(x.id());
    boolean_constraints.resize(std::max(numConstraint(), 2 * boolean.size()));
    
    if (d1 != Constant::NoEdge<T>)
        post(new EdgeConstraint<T>(*this, boolean.getLiteral(true, x.id())));
    
    if (d2 != Constant::NoEdge<T>)
        post(new EdgeConstraint<T>(*this, boolean.getLiteral(false, x.id())));
    
    return x;
}

//template <typename T>
//BooleanVar<T> Solver<T>::newDisjunct(const DistanceConstraint<T> &d1,
//                                     const DistanceConstraint<T> &d2,
//                                     const var_t opt) {
//  auto x{boolean.newDisjunct(d1, d2)};
//  clauses.newBooleanVar(x.id());
//  boolean_constraints.resize(std::max(numConstraint(), 2 * boolean.size()));
//
//  post(new OptionalEdgeConstraint<T>(*this, x.id(), opt));
//
//  return x;
//}

template <typename T> NumericVar<T> Solver<T>::newConstant(const T k) {
    return newNumeric(k, k);
}

template <typename T>
NumericVar<T> Solver<T>::newOffset(NumericVar<T> &x, const T k) {
    return NumericVar<T>(x.id(), k);
}

template <typename T>
NumericVar<T> Solver<T>::newNumeric(const T lb, const T ub) {
    //    auto x{numeric.newVar(-lb, ub)};
    
    if (lb > ub) {
        throw Failure<T>({&bound_exp, Constant::FactHint});
    } else if (lb == ub) {
        return NumericVar(Constant::K, lb);
    } else {
        
        auto x{numeric.newVar(Constant::Infinity<T>)};
        
        changed.reserve(numeric.size());
        clauses.newNumericVar(x.id());
        numeric_constraints.resize(std::max(numConstraint(), 2 * numeric.size()));
        core.newVertex(x.id());
        
        set(geq<T>(x.id(), lb));
        set(leq<T>(x.id(), ub));
        propagate();
        
        return x;
    }
}

// template <typename T> TemporalVar<T> Solver<T>::newTemporal(const T offset) {
//     TemporalVar<T> x{newNumeric().id(), offset};
//     core.newVertex(x);
//     return x;
// }

// template <typename T> NumericVar<T> Solver<T>::newTemporal(const T offset) {
//   NumericVar<T> x{newNumeric().id(), offset};
//   //  core.newVertex(x);
//   return x;
// }

template <typename T>
Interval<T> Solver<T>::newInterval(const T mindur, const T maxdur,
                                   const T earliest_start, const T latest_start,
                                   const T earliest_end, const T latest_end,
                                   const BooleanVar<T> opt) {
    return Interval<T>(*this, mindur, maxdur, earliest_start, latest_start, earliest_end, latest_end, opt);
}


template <typename T>
Interval<T> Solver<T>::between(const NumericVar<T> s, const NumericVar<T> e) {
    Interval<T> i(*this, s, e, e - s, BooleanVar<T>(Constant::True));
    post(i.duration >= 0);
    return i;
}

template <typename T>
Interval<T> Solver<T>::continuefor(const NumericVar<T> s, const NumericVar<T> d) {
    
    //    BooleanVar<T> o{};
    Interval<T> i(*this, s, s+d, d, BooleanVar<T>(Constant::True));
    return i;
}

template <typename T>
Interval<T> Solver<T>::maybe_between(const NumericVar<T> s, const NumericVar<T> e) {
    Interval<T> i(*this, s, e, e - s, newBoolean());
    post(i.duration >= 0);
    return i;
}

template <typename T>
Interval<T> Solver<T>::maybe_continuefor(const NumericVar<T> s, const NumericVar<T> d) {
    
    //    BooleanVar<T> o{};
    Interval<T> i(*this, s, s+d, d, newBoolean());
    return i;
}

template <typename T>
Interval<T> Solver<T>::between(const NumericVar<T> s, const NumericVar<T> e, const BooleanVar<T> optional) {
    Interval<T> i(*this, s, e, e - s, optional);
    post(i.duration >= 0);
    return i;
}

template <typename T>
Interval<T> Solver<T>::continuefor(const NumericVar<T> s, const NumericVar<T> d, const BooleanVar<T> optional) {
    
    //    BooleanVar<T> o{};
    Interval<T> i(*this, s, s+d, d, optional);
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
Explanation<T> Solver<T>::getReason(const index_t i) const {
    return reason[i];
}

template <typename T>
Explanation<T> Solver<T>::getReason(const Literal<T> l) const {
    return reason[getReason(propagationLevel(l))];
}

template <typename T> Literal<T> Solver<T>::getLiteral(const index_t i) const {
    return trail[i];
}

// template <typename T>
// bool Solver<T>::subsumed(const Literal<T> l) const {
//     assert(l.isNumeric());
//
//     return(l.value() > numeric.getLiteral(l.sign(), l.variable()).value());
// }
//
// template <typename T>
// Literal<T> Solver<T>::previousBound(const Literal<T> l) const {
//     return numeric.previousBound(l);
////    std::cout << "get previous bound of " << l << std::endl;
////
////
////    auto i{numeric.lastLitIndex(l.sign(), l.variable())};
////
////
////    std::cout << i << std::endl;
////
////    std::cout << trail[i] << std::endl;
////    std::cout << trail[0] << std::endl;
////
////    assert(i > 0);
////
////    while (getLiteral(i--) != l);
////
////    assert(i >= 0);
////
////    return getLiteral(i);
//}

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
        if (DBG_BOUND and (DBG_TRACE & PROPAGATION)) {
            std::cout << "set " << pretty(l) << " @" << trail.size() << " b/c " << e << std::endl;
        }
#endif
        
        reason.emplace_back(e);
        trail.push_back(l);
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
            std::cout << "failure on " << pretty(l) << " @" << trail.size() << " b/c " << e << std::endl;
        }
#endif
        
#ifdef DBG_FAIL
        if (DBG_FAIL) {
            std::cout << "failure on " << pretty(l) << " @" << trail.size()
            << " b/c " << e << std::endl;
        }
#endif
        
        throw Failure<T>(e);
    }
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & PROPAGATION)) {
        std::cout << "set " << pretty(l) << " @" << trail.size() << " b/c " << e << std::endl;
    }
#endif
    
    reason.emplace_back(e);
    trail.push_back(l);
    boolean.set(l);
    
    if (boolean_search_vars.has(l.variable()))
        boolean_search_vars.remove_back(l.variable());
}

template <typename T>
void Solver<T>::boundClosure(const var_t x, const var_t y, const T d,
                             const index_t r) {
    // closure w.r.t. 0 (0 -> x -(d)-> y -> 0)
    
    Explanation<T> e{&graph_exp, static_cast<hint>(r)};
    if (r == Constant::NoIndex)
        e = Constant::Decision<T>;
    
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
    decisions.clear();
    
    //    auto sb{trail.size()};
    //
    //    std::cout << env.level() << std::endl;
    //    std::cout << trail.size() << std::endl;
    //  undo();
    //    std::cout << trail.size() << std::endl;
    //
    //    if(trail.size() != sb) {
    //        std::cout << "here " << trail.size() << "/" << sb << "\n";
    //        exit(1);
    //    }
    
    if (on_solution) {
        restartPolicy.initialize();
    } else {
        restartPolicy.reset();
    }
    
    SearchRestarted.trigger(on_solution);
    
    if (options.verbosity > Options::NORMAL) {
        std::cout << std::setw(10) << "restart ";
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
    
    //    if (env.level() == init_level) {
    //        throw SearchExhausted();
    //    }
    
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
index_t Solver<T>::propagationLevel(const Literal<T> l) const {
    if (l.isNumeric()) {
        return numeric.litIndex(l);
    } else {
        return boolean.litIndex(l);
    }
}

template <typename T>
index_t Solver<T>::decisionLevel(const Literal<T> p) const {
    auto p_lvl{propagationLevel(p)};
    index_t jump{0};
    for (auto d{decisions.rbegin() + jump};
         d != decisions.rend() and propagationLevel(*d) > p_lvl; ++d) {
        ++jump;
    }
    return decisions.size()-jump;
}

template <typename T>
bool Solver<T>::entailedByConflict(Literal<T> p, const index_t p_lvl) const {
    if (p.isNumeric()) {
        auto p_idx{numeric.getConflictIndex(p)};
        if (p_idx != Constant::NoIndex) {
            return conflict[p_idx].value() <= p.value_unsafe();
        }
    }
    return explored[p_lvl]; // boolean.visited(p);
}

template <typename T> void Solver<T>::minimize_clause() {
    //  auto fact_lvl{propagationLevel(decisions[0])};
    //
    //    if(fact_lvl != ground_level)
    //    {
    //        std::cout << ground_level << "/" << fact_lvl << std::endl;
    //        exit(1);
    //    }
    
    if (learnt_clause.empty())
        return;
    
    minimal_clause.clear();
    
#ifdef DBG_MINIMIZATION
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
        std::cout << std::endl << ~learnt_clause[0] << " (UIP:" << decisionLevel(~learnt_clause[0]) << "/" << propagationLevel(~learnt_clause[0]) << "\n";
    }
#endif
    
    index_t i{0};
    
    if(not decisions.empty()) {
        // in this case the UIP is always in the minimized clause
        i = 1;
        minimal_clause.push_back(learnt_clause[0]);
        if(learnt_clause[0].isNumeric())
            numeric.setConflictIndex(~learnt_clause[0], Constant::NoIndex);
        explored[literal_lvl[0]] = false;
    }
    
    for(; i<learnt_clause.size(); ++i) {
        auto p_lvl{literal_lvl[i]};
        
        auto r{reason[p_lvl]};
        bool relevant{true};
        
#ifdef DBG_MINIMIZATION
        if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
            std::cout << ~learnt_clause[i] << "(" << decisionLevel(~learnt_clause[i]) << "/" << propagationLevel(~learnt_clause[i]) << ")" << ":";
        }
#endif
        
        if(r != Constant::Decision<T> and r != Constant::Assumption<T>) {
            auto p{~learnt_clause[i]};
            
            if(learnt_clause[i].isNumeric())
                numeric.setConflictIndex(p, Constant::NoIndex);
            explored[p_lvl] = false;
            
            
            relevant = false;
            
            lit_buffer.clear();
            r.explain(p, lit_buffer);
            
            for(auto q : lit_buffer) {
                
#ifdef DBG_MINIMIZATION
                if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                    std::cout << " " << q;
                }
#endif
                
                auto q_lvl{propagationLevel(q)};
                
                if (q_lvl >= ground_level and
                    not entailedByConflict(q, q_lvl)) {
                    relevant = true;
                    break;
                }
#ifdef DBG_MINIMIZATION
                else if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                    std::cout << " (";
                    if (q_lvl < ground_level)
                        std::cout << "f";
                    if(entailedByConflict(q, q_lvl))
                        std::cout << "e";
                    std::cout << ")";
                }
#endif
            }
        }
        
#ifdef DBG_MINIMIZATION
        else if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
            std::cout << " (decision)";
        }
#endif
        
        if(relevant) {
            assert(p_lvl == literal_lvl[i]);
            
            literal_lvl[minimal_clause.size()] = p_lvl; //literal_lvl[i];
            minimal_clause.push_back(learnt_clause[i]);
        }
        
#ifdef DBG_MINIMIZATION
        if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
            std::cout << (relevant ? " ==> relevant\n" : " ==> useless\n");
        }
#endif
    }
    
#ifdef DBG_MINIMIZATION
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
        if(minimal_clause.size() < learnt_clause.size()) {
            std::cout << "minimized away " << (learnt_clause.size() - minimal_clause.size()) << " literals\n" ;
        }
    }
#endif
    
    std::swap(minimal_clause, learnt_clause);
    literal_lvl.resize(learnt_clause.size());
    
#ifdef DBG_MINIMIZATION
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
        for(size_t i{0}; i<learnt_clause.size(); ++i) {
            assert(propagationLevel(~(learnt_clause[i])) == literal_lvl[i]);
        }
    }
#endif
}

template <typename T> void Solver<T>::analyze(Explanation<T> &e) {
    explored.resize(trail.size(), false);
    conflict.clear();
    literal_lvl.clear();
    learnt_clause.clear();
    auto decision_lvl{(decisions.empty() ? assumption_level
                       : propagationLevel(decisions.back()))};
    //    auto fact_lvl{propagationLevel(decisions[0])};
    
    //    if(fact_lvl != ground_level)
    //    {
    //        std::cout << ground_level << "/" << fact_lvl << std::endl;
    //        exit(1);
    //    }
    
    int num_lit{0};
    index_t li{static_cast<index_t>(trail.size() - 1)};
    Literal<T> l{Contradiction};
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
        std::cout << "analyze conflict:\n";
        //        displayBranches(std::cout);
        for(auto i{ground_level}; i<trail.size(); ++i) {
            if(i == assumption_level)
                std::cout << "** ";
            std::cout << std::setw(4) << i << " " << pretty(trail[i]) << std::endl;
        }
    }
#endif
    
    Explanation<T> &exp = e;
    do {
        int csize{static_cast<int>(conflict.size())};
        
#ifdef DBG_TRACE
        if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
            std::cout << "resolve ";
            if (l == Contradiction) {
                std::cout << "contradiction";
            } else {
                std::cout << pretty(l) << " @" << propagationLevel(l);
            }
            std::cout << " by " << exp << std::endl;
        }
#endif
        
        exp.explain(l, conflict);
        
        
#ifdef DBG_CLPLUS
        if (num_clauses == DBG_CL and cl_file != NULL) {
            *cl_file << "0 " << (conflict.size() + num_lit + 1) << " 0 1 " << numeric.upper(1);
            for(index_t i{0}; i<conflict.size(); ++i) {
                writeLiteral(conflict[i]);
            }
            for(index_t i{li}; i>=decision_lvl; --i) {
                if(explored[i])
                    writeLiteral(trail[i]);
            }
            *cl_file << std::endl;
        }
#endif
        
        for (int i{static_cast<int>(conflict.size()) - 1}; i >= csize;) {
            
            auto p{conflict[i]};
            
            //@TODO: remove
            auto p_lvl{propagationLevel(p)};
            
            //            if (l != Contradiction and p_lvl > propagationLevel(l)) {
            //              std::cout << "bug here " << num_choicepoints << "\n";
            //
            //              std::cout << env.level() << std::endl;
            //
            //              std::cout << " explain " f<< pretty(l) << " @"
            //                        << propagationLevel(l) << " with " << pretty(p) << " @"
            //                        << p_lvl << std::endl;
            //
            //              std::cout << *(clauses[exp.the_hint]) << std::endl;
            //              //                << ": " << *(exp.expl)
            //              //                << std::endl;
            //
            //              exit(1);
            //            }
            
#ifdef DBG_TRACE
            if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                std::cout << " ** " << pretty(p) << " (" << ground_level << "/"
                << p_lvl << "/" << decision_lvl << ") i.e., "
                << getLiteral(p_lvl);
                std::cout.flush();
            }
#endif
            
            //@TODO: need to separate decisions from ground facts in Explanation
            if (p_lvl < ground_level) {
#ifdef DBG_TRACE
                if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                    std::cout << " => ground fact\n";
                }
#endif
                --i;
            } else if (p_lvl < decision_lvl) {
                if (entailedByConflict(p, p_lvl)) {
#ifdef DBG_TRACE
                    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                        std::cout << " => entailed by cut\n";
                    }
#endif
                    --i;
                } else {
                    bool need_add{true};
                    if (p.isNumeric()) {
                        auto idx_p{numeric.getConflictIndex(p)};
                        if (idx_p != Constant::NoIndex) {
                            if (conflict[idx_p].value() > p.value()) {
                                conflict[idx_p].setValue(p.value());
                                literal_lvl[idx_p] = p_lvl;
                            }
                            need_add = false;
                        } else {
                            numeric.setConflictIndex(p, csize);
                        }
                    }
                    if (need_add) {
                        
#ifdef DBG_TRACE
                        if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                            std::cout << " => add to confict [";
                            for (int z{0}; z < csize; ++z) {
                                std::cout << " " << z << " " << conflict[z];
                                std::cout.flush();
                            }
                            std::cout << " ]\n";
                        }
#endif
                        
                        if(not p.isNumeric())
                            explored[p_lvl] = true;
                        std::swap(conflict[csize], conflict[i]);
                        ++csize;
                        literal_lvl.push_back(p_lvl);
                    } else {
                        --i;
#ifdef DBG_TRACE
                        if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                            std::cout << " => update confict [";
                            for (int z{0}; z < csize; ++z) {
                                std::cout << " " << z << " " << conflict[z];
                            }
                            std::cout << " ]\n";
                        }
#endif
                    }
                }
            } else if (explored[p_lvl]) {
#ifdef DBG_TRACE
                if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                    std::cout << " => already explored\n";
                }
#endif
                --i;
            } else {
                
                explored[p_lvl] = true;
                ++num_lit;
                --i;
                
#ifdef DBG_TRACE
                if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                    auto count{0};
                    std::cout << " => to explore [ ";
                    for (index_t z{li}; z >= decision_lvl; --z) {
                        if (explored[z] or z == p_lvl) {
                            std::cout << " " << trail[z];
                            std::cout.flush();
                            ++count;
                        }
                    }
                    std::cout << "]\n";
                    if (count != num_lit) {
                        std::cout << "bug " << count << "/" << num_lit << "\n";
                        //                      exit(1);
                    }
                }
#endif
            }
        }
        
        conflict.resize(csize);
        
#ifdef DBG_TRACE
        if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
            std::cout << "num_lit: " << num_lit << std::endl;
        }
#endif
        
        if(num_lit > 0) {
            while (not explored[li--])
                ;
            
            l = trail[li + 1];
            exp = reason[li + 1];
            
            explored[li + 1] = false;
            
#ifdef DBG_TRACE
            if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                std::cout << (li + 1) << "//" << decision_lvl << std::endl;
            }
#endif
        }
        
        --num_lit;
        
    } while (num_lit > 0 and li > ground_level); // or l.isNumeric());
    
    //    std::cout << "conflict.size() = " << conflict.size() << " numlit = " << num_lit << std::endl;
    
    
    if(num_lit >= 0) {
        bool need_add{true};
        if (l.isNumeric()) {
            auto idx_l{numeric.getConflictIndex(l)};
            if (idx_l != Constant::NoIndex) {
                if (conflict[idx_l].value() > l.value_unsafe()) {
                    conflict[idx_l].setValue(l.value_unsafe());
                    literal_lvl[idx_l] = li + 1;
                }
                need_add = false;
                
#ifdef DBG_TRACE
                if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                    std::cout << " UIP (updated): " << conflict[idx_l] << std::endl;
                }
#endif
            }
        }
        
        if (need_add) {
            conflict.push_back(l);
            literal_lvl.push_back(li + 1);
            
#ifdef DBG_TRACE
            if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                std::cout << " UIP (added): " << l << std::endl;
            }
#endif
        }
    }
    
    //    else {
    //        std::cout << "conflict.size() = " << conflict.size() << std::endl;
    //    }
    
    assert(conflict.size() == literal_lvl.size());
    
    std::sort(literal_lvl.begin(), literal_lvl.end(),
              [&](const index_t a, const index_t b) { return a > b; });
    
    for (size_t i{0}; i < literal_lvl.size(); ++i) {
        auto p{trail[literal_lvl[i]]};
        if (p.isNumeric()) {
            auto idx_p{numeric.getConflictIndex(p)};
            if (idx_p != Constant::NoIndex)
                p = conflict[idx_p];
        }
        learnt_clause.push_back(~p);
    }
    
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
        
        displayDomains(std::cout);
        
        std::cout << "\nlearn clause\n";
        for (auto l : learnt_clause) {
            std::cout << " " << pretty(l) << " @" << propagationLevel(~l) << ": "
            << (l.isNumeric() ? numeric.falsified(l)
                : boolean.falsified(l));
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
#endif
    
    
    if (options.minimization > 0) {
        minimize_clause();
    }
    
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
        std::cout << "\nminimized clause\n";
        for (auto l : learnt_clause) {
            std::cout << " " << pretty(l) << " @" << propagationLevel(~l) << ": "
            << (l.isNumeric() ? numeric.falsified(l)
                : boolean.falsified(l));
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
#endif
    
    for (auto p : learnt_clause)
        if (p.isNumeric()) {
            numeric.setConflictIndex(~p, Constant::NoIndex);
        }
    
    for (auto i : literal_lvl) {
        explored[i] = false;
    }
    
#ifdef DBG_CL
    if (cl_file != NULL) {
        *cl_file << "0 " << (learnt_clause.size() + 1) << " 0 1 " << numeric.upper(1);
        for (auto p : learnt_clause) {
            writeLiteral(~p);
        }
        *cl_file << std::endl;
    }
#endif
}

template <typename T> bool Solver<T>::isAssertive(std::vector<Literal<T>> &conf) const {
    if(conf[0].isNumeric() and numeric.falsified(conf[0])) {
        return false;
    }
    if(conf[0].isNumeric() and numeric.satisfied(conf[0])) {
        return false;
    }
    if(not conf[0].isNumeric() and boolean.falsified(conf[0])) {
        return false;
    }
    if(not conf[0].isNumeric() and boolean.satisfied(conf[0])) {
        return false;
    }
    
    for (size_t i{1}; i < conf.size(); ++i) {
        if(not conf[i].isNumeric() and not boolean.falsified(conf[i])) {
            return false;
        }
        if(conf[i].isNumeric() and not numeric.falsified(conf[i])) {
            return false;
        }
    }
    return true;
}

template <typename T> void Solver<T>::learnConflict(Explanation<T> &e) {
    
    analyze(e);
    
    //ClauseAdded.trigger(conflict);
    
    int jump{1};
    
    if (learnt_clause.size() == 1 or decisions.size() <= 1) {
        jump = std::max(1, static_cast<int>(decisions.size()));
    } else {
        auto uip_lvl{literal_lvl[1]};
        for (auto d{decisions.rbegin() + jump};
             d != decisions.rend() and propagationLevel(*d) > uip_lvl; ++d) {
            ++jump;
        }
    }
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
        std::cout << "learn clause @lvl" << level() << ":";
        for (auto l : learnt_clause) {
            std::cout << " " << pretty(l) << " (" << propagationLevel(l) << ")";
        }
        std::cout << std::endl; //<< *this << std::endl;
    }
#endif
    //        std::cout << "jump to " << (env.level() - jump) << std::endl;
    //    if(jump == 0)
    //        exit(1);
    
    if ((env.level() - jump) < init_level)
        throw SearchExhausted();
    restoreState(env.level() - jump);
    
    
    decisions.resize(decisions.size() - jump);
    
    
    //    ClauseAdded.trigger(conflict);
    ClauseAdded.trigger(learnt_clause);
    
    
    
    
    //    assert(isAssertive(learnt_clause));
    
    //#ifdef DBG_TRACE
    //    auto cl =
    //#endif
    clauses.add(learnt_clause.begin(), learnt_clause.end(), true);
    
    //#ifdef DBG_TRACE
    //    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
    //        if (clauses.size() > 0 and cl != NULL) {
    //            std::cout << "learn conflict" << *cl << std::endl;
    //        }
    //    }
    //#endif
    
#ifdef DBG_CL
    if (++num_clauses > DBG_CL) {
        std::cout << "exit because of dbg clause limit (#fails = " << num_fails << ", #cpts = " << num_choicepoints << ")\n";
        exit(1);
    }
#endif
}

template <typename T> void Solver<T>::branchRight() {
    
    auto deduction{~decisions.back()};
    
    DeductionMade.trigger(deduction);
    
    if (env.level() <= init_level)
        throw SearchExhausted();
    restoreState(env.level() - 1);
    decisions.pop_back();
    
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

//template <typename T> void Solver<T>::setHeuristic(heuristics::MovableHeuristic<heuristics::PolymorphicHeuristic<T>> h) {
//    heuristic = h;
//}

template <typename T> void Solver<T>::initializeSearch() {
    if(not initialized) {
        start_time = cpu_time();
        post(&clauses);
        
        restartPolicy.initialize();
        if (not heuristic.isValid()) {
            heuristic = heuristics::make_heuristic(*this);
        }
        
        if (/*not initialized and*/ options.verbosity >= Options::QUIET) {
            displayHeader(std::cout);
        }
        
        propag_pointer = 1;
        
        propagate();
        
        assumption_level = ground_level = trail.size();
        
        initialized = true;
    }
}

template<typename T>
template<heuristics::heuristic<T> H>
void Solver<T>::setBranchingHeuristic(H &&h) {
    this->heuristic = std::forward<H>(h);
}

template <typename T> boolean_state Solver<T>::satisfiable() {
    initializeSearch();
    auto satisfiability{search()};
    if (satisfiability == TrueState) {
        boolean.saveSolution();
        numeric.saveSolution();
        
        ++num_solutions;
        SolutionFound.trigger(*this);
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
    conflict.clear();
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
                    conflict.push_back(l);
                }
            } else {
                if (boolean.satisfied(l)) {
                    ++lsat;
                } else if (boolean.falsified(l)) {
                    ++lunsat;
                } else {
                    ++lundef;
                    conflict.push_back(l);
                }
            }
        }
        if (lsat == 0) {
            std::cout << "cl" << i << ": " << lundef << "/" << lunsat << "/" << lsat
            << ":";
            for (auto l : conflict) {
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
    initializeSearch();
    //    if(options.verbosity >= Options::QUIET)
    //    displayHeader(std::cout);
    
    while (objective.gap() and not KillHandler::instance().signalReceived() and
           not(num_fails >= options.search_limit)) {
        auto satisfiability = search();
        if (satisfiability == TrueState) {
            auto best{objective.value(*this)};
            if (options.verbosity >= Options::NORMAL) {
                std::cout << std::setw(10) << best;
                displayProgress(std::cout);
            }
            
            objective.apply(best, *this);
            boolean.saveSolution();
            numeric.saveSolution();
            ++num_solutions;
            SolutionFound.trigger(*this);
            restart(true);
            try {
                objective.setPrimal(best, *this);
            } catch (Failure<T> &f) {
                objective.setDual(objective.primalBound());
            }
        } else if (satisfiability == FalseState) {
            objective.setDual(objective.primalBound());
        }
    }
    
    if (options.verbosity >= Options::QUIET)
        displaySummary(std::cout, (objective.gap() > 0 ? "killed" : "optimal"));
}

template <typename T>
template <typename S, typename A>
void Solver<T>::largeNeighborhoodSearch(S &objective, A &relaxationPolicy) {
    objective.X.extract(*this);
    objective_var = objective.X.id();
    initializeSearch();
    
    //    Solution<T> current_best;
    
    if (not boolean.hasSolution()) {
        auto satisfiability{search()};
        if (satisfiability == TrueState) {
            
//            std::cout << "INIT SOL\n";
            
            auto best{objective.value(*this)};
            if (options.verbosity >= Options::NORMAL) {
                std::cout << std::setw(10) << best;
                displayProgress(std::cout);
            }
            boolean.saveSolution();
            numeric.saveSolution();
            
            ++num_solutions;
            SolutionFound.trigger(*this);
            //        current_best.load(*this);
            restart(true);
            try {
                objective.setPrimal(best, *this);
            } catch (Failure<T> &f) {
                objective.setDual(objective.primalBound());
            }
        } else {
            if(options.verbosity >= Options::QUIET)
                displaySummary(std::cout, "unsatisfiable");
            return;
        }
    }
    
    std::vector<Literal<T>> assumptions;
    
    while (objective.gap() and not KillHandler::instance().signalReceived()) {
        relaxationPolicy.select(assumptions);
        
//        std::cout << " ASSUMPTIONS:";
//        for(auto l : assumptions)
//            std::cout << " " << pretty(l);
//        std::cout << std::endl;
        
        auto satisfiability{FalseState};
        
        try {
            makeAssumptions(assumptions.begin(), assumptions.end());
            satisfiability = UnknownState;
        }
        catch(Failure<T>& f) {
//            std::cout << "assumptions are inconsistent under the current UB\n";
        }
        
        if(satisfiability == UnknownState)
            satisfiability = search();
            
        if (satisfiability == TrueState) {
            auto best{objective.value(*this)};
            if (options.verbosity >= Options::NORMAL) {
                std::cout << std::setw(10) << best;
                displayProgress(std::cout);
            }
            boolean.saveSolution();
            numeric.saveSolution();
            
            ++num_solutions;
            SolutionFound.trigger(*this);
            
            relaxationPolicy.notifySuccess();
            restoreState(0);
        
            try {
                objective.setPrimal(best, *this);
                propagate();
            } catch (Failure<T> &f) {
                objective.setDual(objective.primalBound());
            }
            ground_level = trail.size();
        } else {
            
//            std::cout << "failed to improve the current best\n";
            
//            std::cout << "learn clause:\n";
//            for(auto l : learnt_clause) {
//                std::cout << pretty(l) << std::endl;
//            }
//            std::cout << std::endl;
            
//            std::cout << learnt_clause.size() << "/" << assumptions.size() << std::endl;
            
            relaxationPolicy.notifyFailure();
            restoreState(0);
        }
    }
    
    if (options.verbosity >= Options::QUIET)
        displaySummary(std::cout, (objective.gap() > 0 ? "killed" : "optimal"));
}

//
//template <typename T>
//template <typename S>
//void Solver<T>::optimize(S &objective) {
//    objective.X.extract(*this);
//  initializeSearch();
//    
//    // first, find an initial
//    
//
//  while (objective.gap() and not KillHandler::instance().signalReceived()) {
//    auto satisfiability = search();
//    if (satisfiability == TrueState) {
//      auto best{objective.value(*this)};
//      if (options.verbosity >= Options::NORMAL) {
//        std::cout << std::setw(10) << best;
//        displayProgress(std::cout);
//      }
//        objective.apply(best, *this);
//      boolean.saveSolution();
//      numeric.saveSolution();
//      restart(true);
//      try {
//        objective.setPrimal(best, *this);
//      } catch (Failure<T> &f) {
//        objective.setDual(objective.primalBound());
//      }
//    } else if (satisfiability == FalseState) {
//      objective.setDual(objective.primalBound());
//    }
//  }
//
//    if(options.verbosity >= Options::QUIET)
//  displaySummary(std::cout, (objective.gap() > 0 ? "killed" : "optimal"));
//}

template <typename T>
template <typename IterLit>
 void Solver<T>::makeAssumptions(IterLit beg_a, IterLit end_a) {
    initializeSearch();
    saveState();
    for (auto a{beg_a}; a!=end_a; ++a)
      set(*a);
    propagate();
    assumption_level = numLiteral()-1;
}

template <typename T> boolean_state Solver<T>::search() {

  //    assumption_level =
  init_level = env.level();
  boolean_state satisfiability{UnknownState};
  while (satisfiability == UnknownState and
         not KillHandler::instance().signalReceived() and
         not(num_fails >= options.search_limit)) {
      
    try {
#ifdef DBG_TRACE
      if (DBG_BOUND) {
        std::cout << "--- propag [i=" << num_choicepoints << "] ---\n";
        printTrace();
      }
#endif

      PropagationInitiated.trigger(*this);
      propagate();
      PropagationCompleted.trigger(*this);

      // make a checkpoint
      saveState();

      // all resource constraints are accounted for => a solution has been found
      if (boolean_search_vars.empty() /* and numeric_search_vars.empty()*/) {
        satisfiability = TrueState;

#ifdef DBG_TRACE
        if (DBG_BOUND) {
          std::cout << "--- new solution [i=" << num_choicepoints << "] ---\n";
          printTrace();
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
                    << init_level << ") [i=" << num_choicepoints << "] (|trail|=" << trail.size() << ")---\n";
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
        f.reason.explain(Solver<T>::Contradiction, conflict);
        for (auto l : conflict) {
          std::cout << pretty(l) << std::endl;
        }
        conflict.clear();
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
    }
  }

  return satisfiability;
}

template <typename T> void tempo::Solver<T>::propagate() {

#ifdef DBG_TRACE
  if (DBG_BOUND and (DBG_TRACE & PROPAGATION) and
      static_cast<index_t>(propag_pointer) < trail.size()) {
    std::cout << "propagate\n";
    for (index_t i{propag_pointer}; i < trail.size(); ++i)
      std::cout << " * " << pretty(trail[i]) << std::endl;
    std::cout << std::endl;
  }
#endif

  index_t p_index{static_cast<index_t>(propag_pointer)};

  clauses.clearTriggers();

  while (not propagation_queue.empty() or trail.size() > p_index) {

    while (trail.size() > p_index) {
      ++num_literals;
      Literal<T> l{trail[p_index]};
      auto culprit{reason[p_index].expl};

#ifdef DBG_TRACE
      if (DBG_BOUND and (DBG_TRACE & QUEUE)) {
        std::cout << "triggers for (" << l << ") b/c " << culprit->id() << "/"
                  << reason[p_index] << std::endl;
      }
#endif

        //TODO: not sure why it is better to do it like this than with the standard constraint queue system (PRIORITY?)
            if (not l.isNumeric())
              clauses.unit_propagate_boolean(l);
            else {
              if (numeric.lastLitIndex(l.sign(), l.variable()) > p_index) {
                ++p_index;
                continue;
                //                    std::cout << "subsumed\n";
              }
            }

//            else if(options.full_up and ) {
//
//            }
            //            else if (2*env.level() <
            //            static_cast<int>(avg_fail_level) and clauses.notify(l,
            //            0)) {
            //
            //#ifdef DBG_TRACE
            //              if (DBG_BOUND and (DBG_TRACE & QUEUE)) {

            //                std::cout << " -" << clauses << " (" <<
            //                clauses.id() << ")"
            //                          << std::endl;
            //              }
            //#endif
            //
            ////                need_up_num = true;
            //              propagation_queue.activate(&clauses);
            //            }

            const std::vector<int> &cid =
                (l.isNumeric() ? numeric_constraints[l]
                               : boolean_constraints[l]);
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
    //      if(need_up_num and not up_done and propagation_queue.empty() and
    //      trail.size() == p_index) {
    //          clauses.propagate();
    //          up_done = true;
    //      }
  }

  propag_pointer = p_index;
}

template <typename T> int Solver<T>::saveState() {
  assert(propag_pointer == static_cast<index_t>(trail.size()));

  int lvl{env.level()};
  env.save();
  ReversibleObject::save();
  return lvl;
}

template <typename T> void Solver<T>::restoreState(const int l) {
//  if (l < init_level) {
//    throw SearchExhausted();
//  } else {
    env.restore(l);
//  }
}

template <typename T> void Solver<T>::undo() {
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

//for(auto x : core) {
//    for(auto e : core[x]) {
//        auto y{int(e)};
//        auto d{e.label()};
////        map[x] -- d --> map[y]
//    }
//}

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
            exp.explain(Solver<T>::Contradiction, conflict);
            for (auto l : conflict) {
              std::cout << pretty(l) << std::endl;
            }
            conflict.clear();
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
}

template <typename T> void Solver<T>::relax(Constraint<T> *con) {
  if (boolean_constraints.indegree(con->id()) > 0) {
    boolean_constraints.remove(con->id(), IN);
  }
  if (numeric_constraints.indegree(con->id()) > 0) {
    numeric_constraints.remove(con->id(), IN);
  }
    
//    exit(1);
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
//      return numeric_constraints.indegree(c)-1;
      return numeric_constraints.rank(l).back();
  } else {
      if(boolean_constraints[l].empty() or boolean_constraints[l].back() != c)
        boolean_constraints.add(l, c);
//      return boolean_constraints.indegree(c)-1;
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
template <typename ItTask, typename ItVar>
void Solver<T>::postEdgeFinding(Interval<T> &schedule, ItTask beg_task,
                                ItTask end_task, const ItVar beg_var) {
  post(new DisjunctiveEdgeFinding<T>(*this, schedule, beg_task, end_task,
                                  beg_var
                                  ));
}

template <typename T>
template <typename ItTask, typename ItVar>
void Solver<T>::postTransitivity(Interval<T> &schedule, ItTask beg_task,
                                 ItTask end_task, const ItVar beg_var) {
  post(new Transitivity<T>(*this, schedule, beg_task, end_task,
                                  beg_var
                                  ));
}

template <typename T>
template <typename ItRes>
FullTransitivity<T> *Solver<T>::postFullTransitivity(const ItRes beg_res,
                                                     const ItRes end_res) {
    auto c = new FullTransitivity<T>(*this);
  for (auto res{beg_res}; res != end_res; ++res) {
    c->addResource(res->begDisjunct(), res->endDisjunct());
  }
  post(c);
  return c;
}

template <typename T>
template <typename ItTask, typename ItNVar, typename ItBVar>
void Solver<T>::postCumulative(const NumericVar<T> c, const ItTask beg_task,
                               const ItTask end_task, const ItNVar beg_dem,
                               const ItBVar beg_disj) {
  post(new CumulativeCheck<T>(*this, c, beg_task, end_task, beg_dem, beg_disj));
}

template <typename T>
template <typename ItTask, typename ItNVar, typename ItBVar>
void Solver<T>::postStrongEdgeFinding(
    const Interval<T> s, const NumericVar<T> c, const ItTask beg_task,
    const ItTask end_task, const ItNVar beg_dem, const ItBVar beg_disj) {
  post(new CumulativeEdgeFinding<T>(*this, s, c, beg_task, end_task, beg_dem,
                                    beg_disj));
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
  os << "   cpu\n";
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


  os << "   " << std::left << cpu << std::right << std::endl;

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
//  for (size_t i{0}; i < decisions.size(); ++i) {
//    os << std::setw(3) << i << ": " << decisions[i] << std::endl;
//  }
//  return os;
    for (auto b{boolean_search_vars.bbegin()}; b!=boolean_search_vars.bend(); ++b) {
      os << " " << pretty(boolean.getLiteral(boolean.isTrue(*b), *b)) ;
    }
    
    return os;
}

template <typename T>
std::ostream &Solver<T>::displayTrail(std::ostream &os) const {
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
      os << "\nd" << (i + 1) << "=[" << pretty(l) << "]";
      ++i;
    } else {
      os << (j > 1 ? ", " : " ") << pretty(l);
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

template <typename T>
std::ostream &Solver<T>::displayPrecedences(std::ostream &os) const {
  os << core;
  return os;
}

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
      os << l << " b/c " << reason[i++] << std::endl;
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
            false, (DBG_TRACE & CLAUSES), false, false, false, false);
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
#endif

template <typename T>
std::ostream &operator<<(std::ostream &os, const Solver<T> &x) {
  return x.display(os);
}
}

#endif

