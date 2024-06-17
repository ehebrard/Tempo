
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
#include "Objective.hpp"
#include "Restart.hpp"
//#include "TemporalNetwork.hpp"
#include "constraints/DisjunctiveEdgeFinding.hpp"
#include "constraints/EdgeConstraint.hpp"
//#include "constraints/Transitivity.hpp"
#include "heuristics/HeuristicManager.hpp"
#include "heuristics/ValueHeuristicsManager.hpp"
#include "heuristics/impl/DecayingEventActivityMap.hpp"
//#include "util/Heap.hpp"
#include "util/KillHandler.hpp"
#include "util/Options.hpp"
#include "util/SubscribableEvent.hpp"


//#define DBG_MINIMIZATION

namespace tempo {

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
    const DistanceConstraint<T> &getEdge(const index_t l) const;
    
    bool hasSemantic(const var_t x) const;
    
    void saveBestSolution() { best_solution = polarity; }
    
    //    bool visited(const Literal<T> l) const;
    //    void markVisited(const Literal<T> l) ;
    //    void unmarkVisited(const Literal<T> l) ;
    
protected:
    Solver<T> &solver;
    
    std::vector<bool> polarity;
    std::vector<bool> best_solution;
    
    std::vector<info_t> edge_index;
    
    std::vector<DistanceConstraint<T>> edges;
    
    std::vector<index_t> propagation_level;
    
    //  std::vector<bool> explored;
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
    
    std::ostream &displayLiteralTrail(std::ostream &os, const bool s,
                                      const var_t x) const {
        for (auto i : bound_index[s][x]) {
            os << " @" << i << ": " << solver.getLiteral(i);
        }
        return os;
    }
    
    //  bool visited(const Literal<T> l) const;
    //  void markVisited(const Literal<T> l) ;
    //  void unmarkVisited(const Literal<T> l) ;
    
    index_t getConflictIndex(const Literal<T> l) const;
    void setConflictIndex(const Literal<T> l, T v);
    
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
    
    // used for learning
    //  std::vector<T> explored_bound[2];
    std::vector<index_t> conflict_index[2];
};

template <typename T = int> class GraphExplainer : public NewExplainer<T> {
    
public:
    GraphExplainer(Solver<T> &s);
    void xplain(const Literal<T>, const hint, std::vector<Literal<T>> &) override;
    std::ostream &print_reason(std::ostream &os, const hint) const override;
    int getType() const override;
    
private:
    Solver<T> &solver;
};

template <typename T = int> class BoundExplainer : public NewExplainer<T> {
    
public:
    BoundExplainer(Solver<T> &s);
    void xplain(const Literal<T>, const hint, std::vector<Literal<T>> &) override;
    std::ostream &print_reason(std::ostream &os, const hint) const override;
    int getType() const override;
    
private:
    Solver<T> &solver;
};

template <typename T>
void GraphExplainer<T>::xplain(const Literal<T> l, const hint h,
                               std::vector<Literal<T>> &Cl) {
    
    if (l == Solver<T>::Contradiction) {
        auto s{Literal<T>::sgn(h)};
        auto x{Literal<T>::var(h)};
        //    //        auto lcycle{solver.getLiteral(static_cast<index_t>(h))};
        //    std::cout << "explain contradiction: negative cycle "
        //              << (s == bound::lower ? "lb(x" : "ub(x") << x << ")" <<
        //              std::endl;
        
        auto end_cycle{x};
        
        //        do {
        //
        //
        //        }
        
        //        exit(1);
        
        do {
            //      auto r_idx{solver.numeric.lastLitIndex(s, x)};
            //            auto e{solver.getReason(r_idx)};
            //
            ////            std::cout << e << std::endl;
            //
            //            auto l_idx{e.the_hint};
            auto le{solver.getLiteral(
                                      solver.getReason(solver.numeric.lastLitIndex(s, x)).the_hint)};
            
            //      std::cout << " ** " << le << std::endl;
            
            if (le.isNumeric()) {
                x = le.variable();
                assert(s == le.sign());
            } else {
                
                Cl.push_back(le);
                //            Cl.push_back(le);
                
                auto c{solver.boolean.getEdge(le)};
                
                //                    std::cout << c << std::endl;
                
                //        if (s == bound::lower) {
                //          assert(c.from == x);
                //        } else {
                //          assert(c.to == x);
                //        }
                
                x = (s == bound::lower ? c.to : c.from);
            }
            
            //        std::cout << " x <- " << x << std::endl;
            
            //            exit(1);
            //
            ////            auto p{bounds.reason[bounds.getIndex(LIT(x,s))].the_hint};
            ////
            ////            if(LTYPE(p) == EDGE_LIT) {
            ////                auto el{sched.getEdgeLiteral(FROM_GEN(p))};
            ////                Cl.push_back(EDGE(el));
            ////                x = (s == LOWER ? sched.getEdge(el).to
            ////                                : sched.getEdge(el).from);
            ////            } else {
            ////                x = EVENT(bounds.getConstraint(FROM_GEN(p)).l);
            ////            }
            //
        } while (x != end_cycle);
        ////        while (EVENT(h) != x);
        
        //    exit(1);
        
    } else {
        
        auto r_idx{static_cast<index_t>(h)};
        auto le{solver.getLiteral(r_idx)};
        
        Cl.push_back(le);
        
        //      std::cout << " ** " << le << std::endl;
        
        //      std::cout << "explain lit " << l << " due to lit " << le ;
        
        if (not le.isNumeric()) {
            auto c{solver.boolean.getEdge(le)};
            //          std::cout << " <-> " << c << " and " << Literal<T>(l.sign(),
            //          (l.sign() == bound::lower ? c.to : c.from), l.value() -
            //          c.distance) << ": " ; if(l.sign() == bound::lower) {
            //              std::cout << "-x" << c.to << " <= " << l.value() -
            //              c.distance << std::endl;
            //          } else {
            //              std::cout << "x" << c.from << " <= " << l.value() -
            //              c.distance << std::endl;
            //          }
            //
            //          std::cout << solver.numeric.strongestLiteral(l.sign(),
            //          (l.sign() == bound::lower ? c.to : c.from)) << std::endl;

            Cl.emplace_back(l.sign(), (l.sign() == bound::lower ? c.to : c.from),
                            l.value() - c.distance, detail::Numeric{});
        }
        
        //      std::cout << std::endl;
        //
        //      exit(1);
        
        //      assert(false);
    }
}

template <typename T>
std::ostream &GraphExplainer<T>::print_reason(std::ostream &os,
                                              const hint) const {
    os << "precedence graph";
    return os;
}

template <typename T> int GraphExplainer<T>::getType() const {
    return CYCLEEXPL;
}

template <typename T>
GraphExplainer<T>::GraphExplainer(Solver<T> &s) : solver(s) {}

template <typename T>
void BoundExplainer<T>::xplain(const Literal<T> l, const hint h,
                               std::vector<Literal<T>> &Cl) {
    
    if (l == Solver<T>::Contradiction) {
        
        //      auto bt{SIGN(h)};
        //      auto x{EVENT(h)};
        
        //      auto r_idx{static_cast<index_t>(h)};
        //      auto le{solver.getLiteral(r_idx)};
        //
        //
        //      Cl.push_back(le);
        
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
        NewExplanation<T> exp;
        
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
        
        //        std::cout << le << " and " << ~le << std::endl;
        
        Cl.push_back(le);
        exp.explain(~le, Cl);
        
        //    auto le{solver.getLiteral(std::min(lidx, uidx))};
        //    Cl.push_back(le);
        //    Cl.push_back(~le);
        
    } else {
        std::cout << "explain lit " << l << " due to constraint "
        << solver.boolean.getEdge(static_cast<index_t>(h)) << std::endl;
        
        exit(1);
    }
}

template <typename T>
std::ostream &BoundExplainer<T>::print_reason(std::ostream &os,
                                              const hint) const {
    os << "collapsed bounds";
    return os;
}

template <typename T> int BoundExplainer<T>::getType() const {
    return BOUNDEXPL;
}

template <typename T>
BoundExplainer<T>::BoundExplainer(Solver<T> &s) : solver(s) {}

template <typename T = int> class Solver : public ReversibleObject {
    
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
    mutable SubscribableEvent<NewExplanation<T> &>
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
    
    // get the explanation for the i-th literal
    NewExplanation<T> getReason(const index_t i) const;
    
    // get the explanation for the i-th literal
    NewExplanation<T> getReason(const Literal<T> l) const;
    
    // get the index in the propagation queue of the last Literal involving
    // variable x
    index_t propagationLevel(const Literal<T> l) const;
    index_t decisionLevel(const Literal<T> l) const;
    bool entailedByConflict(Literal<T>, const index_t) const;
    
    void set(Literal<T> l, const NewExplanation<T> &e = Constant::Decision<T>);
    //  void setNumericNoUpdate(Literal<T> l,
    //                  const NewExplanation<T> &e);
    void setNumeric(Literal<T> l,
                    const NewExplanation<T> &e = Constant::Decision<T>,
                    const bool do_update = true);
    void setBoolean(Literal<T> l,
                    const NewExplanation<T> &e = Constant::Decision<T>);
    void set(const DistanceConstraint<T> &c, const index_t r = Constant::NoIndex);
    void boundClosure(const var_t x, const var_t y, const T d, const index_t r);
    template <typename G>
    void update(const bool bt, const int s, const G &neighbors);
    //@}
    
    void post(NewConstraint<T> *);
    void relax(NewConstraint<T> *);
    void wake_me_on(const Literal<T>, const int);
    template <typename X> void addToSearch(const X &x);
    
    
    template <typename ItTask, typename ItVar>
    void postEdgeFinding(Job<T>& schedule, const ItTask beg_task, const ItTask end_task,
                         const ItVar beg_var
      );
    
    /**
     * @name search
     */
    //@{
    void trigger(const Literal<T> l);
    void propagate();
    boolean_state search();
    
    boolean_state satisfiable();
    template <typename S> void optimize(S &objective);
    
    void restart(const bool on_solution = false);
    void backtrack(NewExplanation<T> &e);
    void branchRight();
    void learnConflict(NewExplanation<T> &e);
    void analyze(NewExplanation<T> &e);
    void minimize();
    
    //    bool visited(const Literal<T> l) const {
    //        if(l.isNumeric()) {
    //            return numeric.visited(l);
    //        }
    //        return boolean.visited(l);
    //    }
    //    void markVisited(const Literal<T> l) {
    //        if(l.isNumeric()) {
    //            numeric.markVisited(l);
    //        } else {
    //            boolean.markVisited(l);
    //        }
    //    }
    //    void unmarkVisited(const Literal<T> l) {
    //        if(l.isNumeric()) {
    //            numeric.unmarkVisited(l);
    //        } else {
    //            boolean.unmarkVisited(l);
    //        }
    //    }
    
    //    int takeDecision(const Literal<T> l);
    
    const SparseSet<var_t, Reversible<size_t>> &getBranch() const {
        return boolean_search_vars;
    }
    
    int saveState();
    void restoreState(const int);
    void undo() override;
    //@}
    
    //  void xplain(const Literal<T>, const hint,
    //              std::vector<Literal<T>> &) override; // {}
    //  std::ostream &print_reason(std::ostream &os,
    //                             const hint) const override; // { return os; }
    //  //    int getType() const override;// { return CYCLEEXPL; }
    
    BacktrackEnvironment &getEnv() { return env; }
    const Options &getOptions() const { return options; }
    
    /**
     * @name printing and trace
     */
    //@{
    std::ostream &display(std::ostream &os, const bool dom = true,
                          const bool bra = true, const bool sva = false,
                          const bool pre = false, const bool cla = false,
                          const bool bgr = false, const bool ngr = false,
                          const bool con = false, const bool trl = false) const;
    std::ostream &displayTrail(std::ostream &os) const;
    std::ostream &displayDomains(std::ostream &os) const;
    std::ostream &displayBranches(std::ostream &os) const;
    std::ostream &displayVariables(std::ostream &os) const;
    std::ostream &displayConstraints(std::ostream &os) const;
    
    std::ostream &displayProgress(std::ostream &os) const;
    std::ostream &displayHeader(std::ostream &os, const int width = 59) const;
    std::ostream &displaySummary(std::ostream &os, std::string msg) const;
    
    void check_clauses(const char* msg);
    void check_clause(const index_t i);
    
    std::string pretty(const Literal<T> l) const;
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
    
    std::optional<heuristics::HeuristicManager<T>> heuristic;
    std::optional<heuristics::ValueHeuristicsManager> valueHeuristic;
    //  RestartPolicy *restart_policy = nullptr;
    //  unsigned int restart_limit{static_cast<unsigned int>(-1)};
    RestartManager<Solver<T>> restartPolicy;
    
public:
    SparseSet<var_t, Reversible<size_t>> boolean_search_vars;
    SparseSet<var_t, Reversible<size_t>> numeric_search_vars;
    
private:
    GraphExplainer<T> graph_exp;
    BoundExplainer<T> bound_exp;
    
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
    std::vector<Literal<T>> lit_buffer;
    std::vector<Literal<T>> conflict;
    std::vector<index_t> literal_lvl;
    std::vector<Literal<T>> learnt_clause;
    std::vector<Literal<T>> minimal_clause;
    //    std::vector<Literal<T>> probe;
    
#ifdef DBG_TRACE
    void printTrace() const;
#endif
    
    heuristics::impl::EventActivityMap<T> *activityMap{NULL};
    
public:
    void setActivityMap(heuristics::impl::EventActivityMap<T> *map) {
        activityMap = map;
    }
    heuristics::impl::EventActivityMap<T> *getActivityMap() {
        return activityMap;
    }
    
    double start_time;
    
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
    bool initialized{false};
    void initializeSearch();
    
    static constexpr Literal<T> Contradiction = makeBooleanLiteral<T>(false, Constant::NoVarx, 0);
    
    double looseness(const Literal<T> &l) const;
    
    bool isAssertive(std::vector<Literal<T>> &conf) const;
    
#ifdef DBG_CL
private:
    //  TemporalNetwork<T> cut;
    
    std::ofstream *cl_file{NULL};
    int num_clauses{0};
    
    void writeLiteral(const Literal<T> l) const;
#endif
};


#ifdef DBG_TRACE
template <typename T> void Solver<T>::printTrace() const {
    if (DBG_TRACE & SEARCH) {
        display(std::cout, (DBG_TRACE & DOMAINS), (DBG_TRACE & BRANCH), false,
                false, (DBG_TRACE & CLAUSES), false, false, false, false);
    }
}
#endif

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
    return edge_index[x] != Constant::NoSemantic;
}

template <typename T> BooleanStore<T>::BooleanStore(Solver<T> &s) : solver(s) {
    edges.push_back(Constant::NoEdge<T>);
    edges.push_back(Constant::NoEdge<T>);
}

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

template <typename T> size_t BooleanStore<T>::size() const {
    return polarity.size() / 2;
}

template <typename T> BooleanVar<T> BooleanStore<T>::newVar(const info_t s) {
    BooleanVar<T> x{static_cast<var_t>(size())};
    
    propagation_level.push_back(Constant::NoIndex);
    
    polarity.push_back(false);
    polarity.push_back(false);
    
    //  explored.push_back(false);
    
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
    propagation_level[l.variable()] = (solver.numLiteral() - 1);
    polarity[l] = true;
    if (l.hasSemantic()) {
        assert(l.constraint() == (edge_index[l.variable()] + l.sign()));
        solver.set(edges[l.constraint()],
                   static_cast<index_t>(solver.numLiteral() - 1));
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
    
    //  explored_bound[bound::lower].push_back(-Constant::Infinity<T>);
    //  explored_bound[bound::upper].push_back(-Constant::Infinity<T>);
    //
    conflict_index[bound::lower].push_back(Constant::NoIndex);
    conflict_index[bound::upper].push_back(Constant::NoIndex);
    
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
    bound_index[s][v].push_back(static_cast<index_t>(solver.numLiteral() - 1));
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
    
    //      std::cout << "\ncompute lit index of " << l << std::endl;
    //      displayLiteralTrail(std::cout, l.sign(), l.variable());
    //      std::cout << std::endl;
    
    auto i{bound_index[l.sign()][l.variable()].rbegin()};
    while (solver.getLiteral(*(i + 1)).value() <= l.value())
        ++i;
    
    //    std::cout << " --> " << *i << std::endl;
    
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
numeric_constraints(&env), restartPolicy(*this),
boolean_search_vars(0, &env), numeric_search_vars(0, &env),
graph_exp(*this), bound_exp(*this) {
    trail.emplace_back(Constant::NoVarx, Constant::Infinity<T>, detail::Numeric{});
    reason.push_back(Constant::Decision<T>);
    seed(options.seed);
    
#ifdef DBG_CL
    if (options.dbg_file != "")
        cl_file = new std::ofstream(options.dbg_file, std::ofstream::out);
#endif
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
    clauses.newBooleanVar(x.id());
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

template <typename T>
NewExplanation<T> Solver<T>::getReason(const index_t i) const {
    return reason[i];
}

template <typename T>
NewExplanation<T> Solver<T>::getReason(const Literal<T> l) const {
    return reason[getReason(propagationLevel(l))];
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
    
    if (l.isNumeric()) {
        setNumeric(l, e);
    } else {
        setBoolean(l, e);
    }
}

template <typename T>
void Solver<T>::setNumeric(Literal<T> l, const NewExplanation<T> &e,
                           const bool do_update) {
    
    //#ifdef DBG_TRACE
    //  if (DBG_BOUND and (DBG_TRACE & PROPAGATION)) {
    //    std::cout << "set num literal " << l << " b/c " << e << "?" <<
    //    std::endl;
    //  }
    //#endif
    
    if (not numeric.satisfied(l)) {
        
#ifdef DBG_TRACE
        if (DBG_BOUND and (DBG_TRACE & PROPAGATION)) {
            std::cout << "set " << pretty(l) << " b/c " << e << std::endl;
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
            
            throw NewFailure<T>({&bound_exp, static_cast<hint>(l.variable())});
        }
        
        if (do_update) {
            if (l.sign() == bound::upper) {
                update(bound::upper, l.variable(), core);
            } else {
                update(bound::lower, l.variable(), core.backward());
            }
        }
    }
    
    //    std::cout << "here\n";
}

template <typename T>
void Solver<T>::setBoolean(Literal<T> l, const NewExplanation<T> &e) {
    
    //#ifdef DBG_TRACE
    //  if (DBG_BOUND and (DBG_TRACE & PROPAGATION)) {
    //    std::cout << "set bool literal " << l << " b/c " << e << "?" <<
    //    std::endl;
    //  }
    //#endif
    
    assert(not boolean.satisfied(l));
    
    if (boolean.falsified(l)) {
        
#ifdef DBG_TRACE
        if (DBG_BOUND and (DBG_TRACE & PROPAGATION)) {
            std::cout << "failure on " << pretty(l) << " b/c " << e << std::endl;
        }
#endif
        
        throw NewFailure<T>(e);
    }
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & PROPAGATION)) {
        std::cout << "set " << pretty(l) << " b/c " << e << std::endl;
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
    
    NewExplanation<T> e{&graph_exp, static_cast<hint>(r)};
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
    env.restore(init_level);
    decisions.clear();
    
    //    std::cout << "restart\n" ;
    //    displayBranches(std::cout);
    
    undo();
    //    displayBranches(std::cout);
    
    //    exit(1);
    
    if (on_solution) {
        restartPolicy.initialize();
        //    restart_policy->initialize(restart_limit);
        //    restart_limit += num_fails;
        //    //      std::cout << num_fails << " / " << restart_limit << std::endl;
    } else {
        restartPolicy.reset();
        //    restart_policy->reset(restart_limit);
        //    //      std::cout << num_fails << " / " << restart_limit << std::endl;
        //    //    displayStats(std::cout, "             ");
    }
    
    SearchRestarted.trigger();
    
    if (options.verbosity > Options::NORMAL) {
        std::cout << std::setw(13) << "restart ";
        displayProgress(std::cout);
    }
}

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

// template<typename T>
// void Solver<T>::markExplored(const Literal<T> &l) {
//     explored[propagationLevel(l)] = true;
// }

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

//numeric.setConflictIndex(~p, Constant::NoIndex);
//explored[i] = false;

template <typename T> void Solver<T>::minimize() {
    auto fact_lvl{propagationLevel(decisions[0])};
    minimal_clause.clear();
    minimal_clause.push_back(learnt_clause[0]);
    
#ifdef DBG_MINIMIZATION
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
        std::cout << std::endl << ~learnt_clause[0] << " (UIP:" << decisionLevel(~learnt_clause[0]) << "/" << propagationLevel(~learnt_clause[0]) << "\n";
    }
#endif
    
    
    if(learnt_clause[0].isNumeric())
        numeric.setConflictIndex(~learnt_clause[0], Constant::NoIndex);
    explored[literal_lvl[0]] = false;
    
    for(index_t i{0}; ++i<learnt_clause.size();) {
        auto p_lvl{literal_lvl[i]};
        
 
        
        auto r{reason[p_lvl]};
        bool relevant{true};
        
#ifdef DBG_MINIMIZATION
        if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
            std::cout << ~learnt_clause[i] << "(" << decisionLevel(~learnt_clause[i]) << "/" << propagationLevel(~learnt_clause[i]) << ")" << ":";
        }
#endif
        
        if(r != Constant::Decision<T>) {
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
                

                
                if(q_lvl >= fact_lvl and not entailedByConflict(q, q_lvl))
                {
                    relevant = true;
                    break;
                }
#ifdef DBG_MINIMIZATION
                else if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                    std::cout << " (";
                    if(q_lvl < fact_lvl)
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
//        else {
//
//        }
        
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
    
    
    
    
    //    for(auto l : learnt_clause) {
    //        auto p_lvl{propagationLevel(l)};
    //    }
}

template <typename T> void Solver<T>::analyze(NewExplanation<T> &e) {
    explored.resize(trail.size(), false);
    conflict.clear();
    literal_lvl.clear();
    learnt_clause.clear();
    auto decision_lvl{propagationLevel(decisions.back())};
    auto fact_lvl{propagationLevel(decisions[0])};
    
    int num_lit{0};
    index_t li{static_cast<index_t>(trail.size() - 1)};
    Literal<T> l{Contradiction};
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
        std::cout << "analyze conflict:\n";
        displayBranches(std::cout);
    }
#endif
    
    NewExplanation<T> &exp = e;
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
            
#ifdef DBG_TRACE
            if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                std::cout << " ** " << pretty(p) << " (" << fact_lvl << "/" << p_lvl
                << "/" << decision_lvl << ") i.e., " << getLiteral(p_lvl);
            }
#endif
            
            //@TODO: need to separate decisions from ground facts in Explanation
            if (p_lvl < fact_lvl) {
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
                            std::cout << " => add to confict[";
                        }
#endif
                        
                        if(not p.isNumeric())
                            explored[p_lvl] = true;
                        std::swap(conflict[csize], conflict[i]);
                        ++csize;
                        
                        literal_lvl.push_back(p_lvl);
                    }
#ifdef DBG_TRACE
                    else if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                        std::cout << " => update confict[";
                    }
                    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                        for (int z{0}; z < csize; ++z) {
                            std::cout << " " << conflict[z];
                        }
                        std::cout << " ]\n";
                    }
#endif
                }
            } else if (explored[p_lvl]) {
#ifdef DBG_TRACE
                if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                    std::cout << " => already explored\n";
                }
#endif
                --i;
            } else {
                
#ifdef DBG_TRACE
                if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
                    std::cout << " => to explore [ ";
                    for (index_t z{li - 1}; z >= decision_lvl; --z) {
                        if (explored[z] or z == p_lvl)
                            std::cout << " " << trail[z];
                    }
                    std::cout << "]\n";
                }
#endif
                explored[p_lvl] = true;
                ++num_lit;
                --i;
            }
            
            //      if (p_lvl >= fact_lvl) {
            //        //            markVisited(p);
            //        explored[p_lvl] = true;
            //      }
            
        }
        
        conflict.resize(csize);
        
        while (not explored[li--])
            ;
        
        if (li + 1 < decision_lvl) {
            std::cout << env.level() << " " << num_fails << " " << li << std::endl;
        }
        
        assert(li + 1 >= decision_lvl);
        
        l = trail[li + 1];
        exp = reason[li + 1];
        
        explored[li + 1] = false;
        
        --num_lit;
        
    } while (num_lit > 0); // or l.isNumeric());
    
    //    std::cout << "\n";
    
    
    
    //    for(index_t i{1}; i<trail.size(); ++i) {
    //        if(explored[i])
    //            std::cout << "wrong lit " << trail[i] << std::endl;
    //
    //        assert(not explored[i]);
    //    }
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
        std::cout << "hello\n";
    }
#endif
    
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
        //    std::cout << "levels:";
        //    for (auto d : decisions) {
        //      std::cout << " " << propagationLevel(d);
        //    }
        std::cout << "\nlearn clause\n";
        for (auto l : learnt_clause) {
            std::cout << " " << pretty(l) << " @" << propagationLevel(~l);
            //      if (l.isNumeric()) {
            //        std::cout << " t =";
            //        numeric.displayLiteralTrail(std::cout, l.sign(), l.variable());
            //      }
            std::cout << std::endl;
            //            << "/[" << fact_lvl << ".." << decision_lvl << "]\n";
        }
        std::cout << std::endl; //<< *this << std::endl;
    }
#endif
    
    
    if (options.minimization > 0) {
        minimize();
    }
    
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
        //    std::cout << "levels:";
        //    for (auto d : decisions) {
        //      std::cout << " " << propagationLevel(d);
        //    }
        std::cout << "\nminimized clause\n";
        for (auto l : learnt_clause) {
            std::cout << " " << pretty(l) << " @" << propagationLevel(~l);
            //      if (l.isNumeric()) {
            //        std::cout << " t =";
            //        numeric.displayLiteralTrail(std::cout, l.sign(), l.variable());
            //      }
            std::cout << std::endl;
            //            << "/[" << fact_lvl << ".." << decision_lvl << "]\n";
        }
        std::cout << std::endl; //<< *this << std::endl;
    }
#endif
    
    for (auto p : learnt_clause)
        if (p.isNumeric()) {
            numeric.setConflictIndex(~p, Constant::NoIndex);
        }
    
    for (auto i : literal_lvl) {
        explored[i] = false;
    }
    
    //  for (index_t i{1}; i < trail.size(); ++i) {
    //    assert(explored[i] == false);
    //    assert(entailedByConflict(trail[i], i) == false);
    //  }
    

    
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

template <typename T> void Solver<T>::learnConflict(NewExplanation<T> &e) {
    
    analyze(e);
    
    ClauseAdded.trigger(conflict);
    
    //  assert(propagationLevel(conflict[0]) >= propagationLevel(decisions.back()));
    assert(literal_lvl[0] >= propagationLevel(decisions.back()));
    
    int jump{1};
    
    if (learnt_clause.size() == 1 or decisions.size() == 1) {
        jump = static_cast<int>(decisions.size());
    } else {
        //    auto uip_lvl{propagationLevel(conflict[1])};
        //    for (auto d{decisions.rbegin() + jump};
        //         d != decisions.rend() and propagationLevel(*d) > uip_lvl; ++d) {
        //      ++jump;
        //    }
        
        auto uip_lvl{literal_lvl[1]};
        for (auto d{decisions.rbegin() + jump};
             d != decisions.rend() and propagationLevel(*d) > uip_lvl; ++d) {
            ++jump;
        }
        
//        std::cout << "backjump b/c " << pretty(learnt) <<
    }
    
    restoreState(env.level() - jump);
    
    
    decisions.resize(decisions.size() - jump);
    
    //  for (auto &p : conflict) {
    //    p = ~p;
    //  }
    
    //    for(auto i{0}; i<conflict.size(); ++i) {
    //        if(conflict[i] != learnt_clause[i])
    //        std::cout << num_fails << ": " << conflict[i] << " " << learnt_clause[i] << std::endl;
    //        assert(conflict[i] == learnt_clause[i]);
    //    }
    
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & SEARCH)) {
        std::cout << "learn clause";
        for (auto l : learnt_clause) {
            std::cout << " " << pretty(l) << " (" << propagationLevel(l) << ")";
        }
        std::cout << std::endl; //<< *this << std::endl;
    }
#endif
    
    //  if (conflict[0].isNumeric() and numeric.satisfied(conflict[0])) {
    //    std::cout << "learn tautology " << conflict[0] << " ["
    //              << numeric.lower(conflict[0].variable()) << ".."
    //              << numeric.upper(conflict[0].variable()) << "]\n";
    //  }
    
    assert(isAssertive(learnt_clause));
    
    //  assert(not conflict[0].isNumeric() or not
    //  numeric.falsified(conflict[0]));
    //  assert(not conflict[0].isNumeric() or not
    //  numeric.satisfied(conflict[0])); assert(conflict[0].isNumeric() or not
    //  boolean.falsified(conflict[0])); assert(conflict[0].isNumeric() or not
    //  boolean.satisfied(conflict[0]));
    //
    //  for (size_t i{1}; i < conflict.size(); ++i) {
    //    //      if(conflict[i].isNumeric() and not
    //    numeric.falsified(conflict[i])) {
    //    //          std::cout << conflict[i] << ": [" <<
    //    //          numeric.lower(conflict[i].variable()) << ".." <<
    //    //          numeric.upper(conflict[i].variable()) << "]\n";
    //    //      }
    //
    //    assert(conflict[i].isNumeric() or boolean.falsified(conflict[i]));
    //    assert(not conflict[i].isNumeric() or numeric.falsified(conflict[i]));
    //  }
    
    //    //@TODO: replace
    //    if(conflict[0].isNumeric()) {
    //        if(numeric.satisfied(conflict[0]))
    //
    //    } else {
    //
    //    }
    
#ifdef DBG_TRACE
    auto cl =
#endif
    clauses.add(learnt_clause.begin(), learnt_clause.end(), true);
    
#ifdef DBG_TRACE
    if (DBG_BOUND and (DBG_TRACE & LEARNING)) {
        if (clauses.size() > 0 and cl != NULL) {
            std::cout << "learn conflict" << *cl << std::endl;
        }
    }
#endif
    
#ifdef DBG_CL
    if (++num_clauses > DBG_CL) {
        std::cout << "exit because of dbg clause limit (#fails = " << num_fails << ", #cpts = " << num_choicepoints << ")\n";
        exit(1);
    }
#endif
}

template <typename T> void Solver<T>::branchRight() {
    
    //    if (env.level() == init_level) {
    //      throw SearchExhausted();
    //    }
    
    auto deduction{~decisions.back()};
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

template <typename T> void Solver<T>::initializeSearch() {
    if(not initialized) {
        start_time = cpu_time();
        heuristic.emplace(*this, options);
        post(&clauses);
        initialized = true;
    }
    restartPolicy.initialize();
}

template <typename T> boolean_state Solver<T>::satisfiable() {
    if(options.verbosity >= Options::QUIET)
        displayHeader(std::cout);
    initializeSearch();
    auto satisfiability{search()};
    if(options.verbosity >= Options::QUIET)
        displaySummary(std::cout, (satisfiability == True ? "sat " : (satisfiability == False ? "unsat " : "unknown ")));
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

template <typename T>
template <typename S>
void Solver<T>::optimize(S &objective) {

  initializeSearch();
    if(options.verbosity >= Options::QUIET)
    displayHeader(std::cout);

  while (objective.gap() and not KillHandler::instance().signalReceived()) {
    auto satisfiability = search();
    if (satisfiability == True) {
      auto best{objective.value()};
      if (options.verbosity >= Options::NORMAL) {
        std::cout << std::setw(13) << best;
        displayProgress(std::cout);
      }
        boolean.saveBestSolution();
      restart(true);
      try {
        objective.setPrimal(best);
      } catch (NewFailure<T> &f) {
        objective.setDual(objective.primalBound());
      }
    } else if (satisfiability == False) {
      objective.setDual(objective.primalBound());
    }
  }

    if(options.verbosity >= Options::QUIET)
  displaySummary(std::cout, (objective.gap() > 0 ? "killed" : "optimal"));
}

template <typename T> boolean_state Solver<T>::search() {

  init_level = env.level();
  boolean_state satisfiability{Unknown};
  while (satisfiability == Unknown and
         not KillHandler::instance().signalReceived()) {

    ++num_choicepoints;
    try {
#ifdef DBG_TRACE
      if (DBG_BOUND) {
        std::cout << "--- propag [i=" << num_choicepoints << "] ---\n";
          printTrace();
      }
#endif
//        check_clauses("before propag");
        
        
  

      propagate();
        
//        check_clauses("after propag");



//      assert(clauses.consistent() == NULL);



      // make a checkpoint
      saveState();

      // all resource constraints are accounted for => a solution has been found
      if (boolean_search_vars.empty()) {
        satisfiability = True;

#ifdef DBG_TRACE
        if (DBG_BOUND) {
          std::cout << "--- new solution [i=" << num_choicepoints << "] ---\n";
          printTrace();
        }
#endif

      } else {
        //        ++num_choicepoints;

        var_t x = heuristic->nextChoicePoint(*this);
        //        var_t x = *(boolean_search_vars.begin());
        //        //        lit d = valueHeuristic->choosePolarity(x, *this);
        Literal<T> d{boolean.getLiteral((random() % 2), x)};
        decisions.push_back(d);

#ifdef DBG_TRACE
        if (DBG_BOUND) {
          std::cout << "--- search node (lvl=" << env.level() << "/"
                    << init_level << ") [i=" << num_choicepoints << "] (|trail|=" << trail.size() << ")---\n";
          std::cout << " ** take decision " << d << ": " << pretty(d)
                    << " **\n";
          printTrace();
        }
#endif
          
//          if(numeric.lower(30) > 170)
//              std::cout << "before dec " << num_choicepoints << std::endl;
          
        set(d);
          
//          check_clauses("after decision");
          
//          if(numeric.lower(30) > 170)
//              std::cout << "after dec " << num_choicepoints << std::endl;

//#ifdef DBG_TRACE
//        if (DBG_BOUND) {
//            std::cout << "after dec " << d << std::endl;
//          printTrace();
//        }
//#endif
      }
    } catch (NewFailure<T> &f) {
      try {
          
//          if(numeric.lower(30) > 170)
//              std::cout << "before backtrack " << num_choicepoints << std::endl;
//          
        backtrack(f.reason);
          
//          check_clauses("after backtrack");
          
//          if(numeric.lower(30) > 170)
//              std::cout << "after backtrack " << num_choicepoints << std::endl;
          
        if (restartPolicy.limit()) {
          restart();
        }
      } catch (const SearchExhausted &f) {
        satisfiability = False;
      }
    }
  }

  return satisfiability;
}

template <typename T> void tempo::Solver<T>::propagate() {
    
    
#ifdef DBG_TRACE
      if (DBG_BOUND and (DBG_TRACE)) {
          std::cout << "propagate\n";
          for(index_t i{propag_pointer}; i<trail.size(); ++i)
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
//        if() {
            if (not l.isNumeric())
                clauses.unit_propagate(l);
            else if(clauses.notify(l, 0)) {
                
#ifdef DBG_TRACE
        if (DBG_BOUND and (DBG_TRACE & QUEUE)) {
          std::cout << " -" << clauses << " (" << clauses.id() << ")" << std::endl;
        }
#endif
                
                propagation_queue.activate(&clauses);
            }
//        }


      //              std::cout << "propagate " << info_t(l) << std::endl;
      //
      //        std::cout << boolean_constraints << std::endl;

      const std::vector<int> &cid =
          (l.isNumeric() ? numeric_constraints[l] : boolean_constraints[l]);
      const std::vector<unsigned> &rank =
          (l.isNumeric() ? numeric_constraints.rank(l)
                         : boolean_constraints.rank(l));

      //              std::cout << cons.size() << std::endl;
      //              exit(1);

      for (auto i{cid.size()}; i-- > 0;) {
        auto cons{constraints[cid[i]]};
          if (cons->idempotent and culprit == cons) {
//              std::cout << "idempotency\n";
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
      
//      for(var_t x{0}; x < numeric.size(); ++x) {
//          auto l{numeric.strongestLiteral(bound::lower, x)};
//          auto u{numeric.strongestLiteral(bound::upper, x)};
//          
//  //        std::cout << "prop " << l << std::endl;
//          
//          if(l.value() != Constant::Infinity<T>)
//              clauses.unit_propagate(l);
//          
//  //        std::cout << "prop " << u << std::endl;
//          
//          if(u.value() != Constant::Infinity<T>)
//              clauses.unit_propagate(u);
//      }
      
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
  //    if(decisions.size()
  //    decisions.resize(decisions.size() + l - env.level());
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

          //            std::cout << " negative cycle on edge " << u << "-(" <<
          //            edge.label() << ")-> " << v << "\n"; std::cout << "i.e.,
          //            " << getLiteral(edge.stamp()) << std::endl;

          throw NewFailure<T>(
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
template <typename ItTask, typename ItVar>
void Solver<T>::postEdgeFinding(Job<T>& schedule,  ItTask beg_task,  ItTask end_task,
                                   const ItVar beg_var
                                   ) {
  post(new NewDisjunctiveEdgeFinding<T>(*this, schedule, beg_task, end_task,
                                  beg_var
                                  ));
}

// template <typename T>
// void Solver<T>::xplain(const Literal<T>, const hint,
//                        std::vector<Literal<T>> &) {
//
//
// }

// template <typename T>
// std::ostream &Solver<T>::print_reason(std::ostream &os, const hint) const {
//   os << "shortest-path";
//   return os;
// }

template <typename T> double Solver<T>::looseness(const Literal<T> &l) const {
  if (l.isNumeric()) {
      
      if(numeric.falsified(l))
          return 0;
      
    auto lb{numeric.lower(l.variable())};
    auto ub{numeric.upper(l.variable())};
      
//      std::cout << env.level() << std::endl;
      
//      std::cout << *this << std::endl;
      

    if (l.sign() == bound::lower) {
      auto b{-l.value_unsafe()};
//        std::cout << "lower bound (" << ub << " - " << b << ")\n";
////      assert(b >= lb);
//        std::cout << l << " -> " << (static_cast<double>(ub - b) / static_cast<double>(ub - lb)) << std::endl;
//        assert(b <= ub);
        
        double ls{static_cast<double>(ub - b)};
        assert(ls>=0);
//        {
//            assert(numeric.falsified(l));
//            return 0;
//        }
        
        return ls / static_cast<double>(ub - lb);
    } else {
      auto b{l.value_unsafe()};
//        std::cout << "upper bound (" << b << " - " << lb << ")\n";
////      assert(b <= ub);
//        std::cout << l << " -> " << (static_cast<double>(b - lb) / static_cast<double>(ub - lb)) << std::endl;
//        assert(b >= lb);
        
        double ls{static_cast<double>(b - lb)};
        assert(ls>=0);
//        {
//            assert(numeric.falsified(l));
//            return 0;
//        }
        
        return ls / static_cast<double>(ub - lb);

    }
      
//      exit(1);
      
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
std::ostream &Solver<T>::displayHeader(std::ostream &os,
                                          const int width) const {
  os << std::right << std::setw(width)
    << " objective   failures   branches   clauses    size   cpu" << std::endl;
//     << std::left;
//  os << std::setfill('=') << std::setw(width) << "=" << std::setfill(' ')
//     << std::endl << std::right;
  return os;
}

template <typename T>
std::ostream &Solver<T>::displaySummary(std::ostream &os,
                                           std::string msg) const {

//  os << std::setfill('=') << std::setw(60) << "\n" << std::setfill(' ');
  //    auto offset{msg.size()/2};
  //    os << std::setfill('=') << std::setw(36 + offset) << msg << std::setw(36
  //    - offset)<< "\n" << std::setfill(' ');
  os << std::setw(13) << msg;
  displayProgress(os);
//  os << std::setfill('=') << std::setw(60) << "\n" << std::setfill(' ');
  return os;
}

template <typename T>
std::ostream &Solver<T>::displayProgress(std::ostream &os) const {

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
    for(size_t i{0}; i<decisions.size(); ++i)
    {
        os << std::setw(3) << i << ": " << decisions[i] << std::endl;
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
    os << "\nd" << (i+1) << "=[" << pretty(l) << "]";
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

