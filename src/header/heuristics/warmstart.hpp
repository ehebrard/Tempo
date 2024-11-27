/**
* @author Emmanuel Hebrard
* @brief warmstart heuristics
*/

#ifndef WARMSTART_HPP
#define WARMSTART_HPP

#include <vector>
#include <iostream>

#include "Solver.hpp"
#include "heuristics/Greedy.hpp"

namespace tempo::heuristics {
    /**
     * Solver warm start heuristic for disjunctive scheduling problems
     * @tparam T timing type
     * @param S solver instance
     * @param schedule schedule interval top optimize
     * @param intervals tasks in the problem
     * @param ub upper bound
     */
    template<typename T>
    void warmstartDisjunctive(Solver<T> &S, const Interval<T> &schedule, std::vector<Interval<T>> intervals, T &ub) {
        // try to get a better ub with an initial upper bound insertion heuristic
        Greedy greedy_insertion(S);
        greedy_insertion.addIntervals(std::move(intervals));
        //  for (auto &R : resources) {
        //    greedy_insertion.addResource(R.begDisjunct(), R.endDisjunct());
        //  }
        for (auto x: S.boolean_search_vars) {
            greedy_insertion.addVar(x);
        }

        // the insertion heuristic is randomized so multiple runs can be useful
        S.initializeSearch();
        S.propagate();
        for (auto i{0}; i < S.getOptions().greedy_runs; ++i) {
            auto st{S.saveState()};
            auto sat{greedy_insertion.runEarliestStart()};
            //        auto sat{greedy_insertion.runLex()};
            if (sat) {
                if (schedule.getEarliestEnd(S) <= ub) {
                    S.set(schedule.end.before(schedule.getEarliestEnd(S)));
                    S.boolean.saveSolution();
                    S.numeric.saveSolution();
                    ub = schedule.getEarliestEnd(S) - 1;
                    std::cout << std::setw(10) << (ub + 1);
                    S.displayProgress(std::cout);
                }
            }
            S.restoreState(st);
        }

        S.set(schedule.end.before(ub));
    }
}

#endif //WARMSTART_HPP
