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
                S.post(schedule.end.before(schedule.getEarliestEnd(S)));
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




template<typename T>
void warmstartEstSlack(Solver<T> &S, const Interval<T> &schedule, std::vector<Interval<T>> intervals, T &ub) {
    
    T trivial_ub{0};
    
    SparseSet<> remaining_intervals;
    remaining_intervals.reserve(intervals.size());
    remaining_intervals.fill();
    
    for(auto& i : intervals) {
        trivial_ub += i.duration.min(S);
    }
    
    //    for(auto& i : intervals) {
    //        std::cout << i.id()<< ": " << i.start.min(S) << "..(" << i.duration.min(S) << ").." << i.end.max(S) << std::endl;
    //    }
    
    S.initializeSearch();
    S.post(schedule.end <= trivial_ub);
    
    //    for(auto& i : intervals) {
    //        std::cout << i.id()<< ": " << i.start.min(S) << ".." << i.start.max(S) << "..(" << i.duration.min(S) << ").." << i.end.min(S) << ".." << i.end.max(S) << std::endl;
    //    }
    
    
    auto st{S.saveState()};
    
    auto sat{true};
    while(sat and not remaining_intervals.empty()) {
        int best{-1};
        T min_est{Constant::Infinity<T>};
        T min_slack{Constant::Infinity<T>};
        for(auto i : remaining_intervals) {
            auto est{intervals[i].start.min(S)};
            auto slack{intervals[i].start.max(S)};
            if(est < min_est or (est==min_est and slack < min_slack)) {
                best = i;
                min_est = est;
                min_slack = slack;
            }
        }
        
        //        std::cout << "Schedule task " << intervals[best].id() << std::endl;
        
        try {
            S.post(intervals[best].start <= min_est);
        } catch(Failure<T>& f) {
            //            std::cout << "FAILED :(\n";
            sat = false;
        }
        
        remaining_intervals.remove_back(best);
        
        
        //        for(auto& i : intervals) {
        //            std::cout << i.id()<< ": " << i.start.min(S) << ".." << i.start.max(S) << "..(" << i.duration.min(S) << ").." << i.end.min(S) << ".." << i.end.max(S) << std::endl;
        //        }
        
    }
    
    //        auto sat{greedy_insertion.runLex()};
    if (sat) {
        //        std::cout << "SUCCESS!!\n";
        if (schedule.getEarliestEnd(S) <= ub) {
            S.post(schedule.end.before(schedule.getEarliestEnd(S)));
            S.boolean.saveSolution();
            S.numeric.saveSolution();
            ub = schedule.getEarliestEnd(S) - 1;
            std::cout << std::setw(10) << (ub + 1);
            S.displayProgress(std::cout);
        }
    }
    //    else {
    //        std::cout << "Failure" << std::endl;
    //    }
    S.restoreState(st);
}



template<typename T>
void warmstartSlackEst(Solver<T> &S, const Interval<T> &schedule, std::vector<Interval<T>> intervals, T &ub) {
    
    T trivial_ub{0};
    
    SparseSet<> remaining_intervals;
    remaining_intervals.reserve(intervals.size());
    remaining_intervals.fill();
    
    for(auto& i : intervals) {
        trivial_ub += i.duration.min(S);
    }
    
    //    for(auto& i : intervals) {
    //        std::cout << i.id()<< ": " << i.start.min(S) << "..(" << i.duration.min(S) << ").." << i.end.max(S) << std::endl;
    //    }
    
    S.initializeSearch();
    S.post(schedule.end <= trivial_ub);
    
    //    for(auto& i : intervals) {
    //        std::cout << i.id()<< ": " << i.start.min(S) << ".." << i.start.max(S) << "..(" << i.duration.min(S) << ").." << i.end.min(S) << ".." << i.end.max(S) << std::endl;
    //    }
    
    
    auto st{S.saveState()};
    
    auto sat{true};
    while(sat and not remaining_intervals.empty()) {
        int best{-1};
        T min_est{Constant::Infinity<T>};
        T min_slack{Constant::Infinity<T>};
        for(auto i : remaining_intervals) {
            auto est{intervals[i].start.min(S)};
            auto slack{intervals[i].start.max(S)};
            if(slack < min_slack or (slack==min_slack and est < min_est)) {
                best = i;
                min_est = est;
                min_slack = slack;
            }
        }
        
        //        std::cout << "Schedule task " << intervals[best].id() << std::endl;
        
        try {
            S.post(intervals[best].start <= min_est);
        } catch(Failure<T>& f) {
            //            std::cout << "FAILED :(\n";
            sat = false;
        }
        
        remaining_intervals.remove_back(best);
        
        
        //        for(auto& i : intervals) {
        //            std::cout << i.id()<< ": " << i.start.min(S) << ".." << i.start.max(S) << "..(" << i.duration.min(S) << ").." << i.end.min(S) << ".." << i.end.max(S) << std::endl;
        //        }
        
    }
    
    //        auto sat{greedy_insertion.runLex()};
    if (sat) {
        //        std::cout << "SUCCESS!!\n";
        if (schedule.getEarliestEnd(S) <= ub) {
            S.post(schedule.end.before(schedule.getEarliestEnd(S)));
            S.boolean.saveSolution();
            S.numeric.saveSolution();
            ub = schedule.getEarliestEnd(S) - 1;
            std::cout << std::setw(10) << (ub + 1);
            S.displayProgress(std::cout);
        }
    }
    //    else {
    //        std::cout << "Failure" << std::endl;
    //    }
    S.restoreState(st);
}


}

#endif //WARMSTART_HPP
