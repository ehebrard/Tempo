/************************************************
 * Tempodisjunctive scheduling solver
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

#include <iostream>
#include <vector>

#include "Solver.hpp"
#include "heuristics/Greedy.hpp"
#include "constraints/Cardinality.hpp"
#include "util/parsing/tsptw.hpp"

using namespace tempo;

template <typename T>
void warmstart(Solver<T> &S, Interval<T> &schedule,
               std::vector<Interval<T>> intervals,
               std::vector<NoOverlapExpression<>> &resources, T &ub) {
  // try to get a better ub with an initial upper bound insertion heuristic
  Greedy greedy_insertion(S);
  greedy_insertion.addIntervals(intervals);
  for (auto &R : resources) {
    greedy_insertion.addResource(R.begDisjunct(), R.endDisjunct());
  }

  // the insertion heuristic is randomized so multiple runs can be useful
  S.initializeSearch();
  S.propagate();
  for (auto i{0}; i < S.getOptions().greedy_runs; ++i) {
    auto st{S.saveState()};
    auto sat{greedy_insertion.runEarliestStart()};
    if (sat) {
      if (schedule.getEarliestEnd(S) <= ub) {
        S.boolean.saveSolution();
        ub = schedule.getEarliestEnd(S) - 1;
        std::cout << std::setw(10) << (ub + 1);
        S.displayProgress(std::cout);
      }
    }
    S.restoreState(st);
  }

  // set the ub (again)
  S.set(schedule.end.before(ub));
}

template<typename T>
std::string prettyJob(const Interval<T>& task, const Solver<T>& S, const bool dur_flag) {
  std::stringstream ss;

  auto est{S.numeric.lower(task.start)};
  auto lst{S.numeric.upper(task.start)};
  auto ect{S.numeric.lower(task.end)};
  auto lct{S.numeric.upper(task.end)};

    
    if(S.boolean.value(task.exist)) {
        ss << "[";
        
        if (est == lst)
            ss << est;
        else
            ss << est << "-" << lst;
        ss << "..";
        if (ect == lct)
            ss << ect;
        else
            ss << ect << "-" << lct;
        ss << "]";
        
        if (dur_flag) {
            auto pmin{S.numeric.lower(task.duration)};
            auto pmax{S.numeric.upper(task.duration)};
            ss << " (" << pmin << "-" << pmax << ")";
        }
    } else {
        ss << "removed";
    }

  return ss.str();
}

template<typename T>
void printJobs(const Solver<T>& S, const std::vector<Interval<T>>& intervals) {
  int i{0};
  for (auto task : intervals) {
    std::cout << "job" << ++i << " (" << task.id()
              << "): " << prettyJob(task, S, true) << std::endl;
  }
}

template<typename T>
void printResources(const Solver<T>& S, const std::vector<Interval<T>>& intervals, const std::vector<std::vector<size_t>>& resource_tasks, std::vector<std::vector<std::vector<T>>>& resource_transitions) {
    int i{0};
    
    for (auto &tasks : resource_tasks) {
        
        if (tasks.size() > 1) {
            
            std::vector<index_t> order;
            for(index_t j{0}; j<static_cast<index_t>(tasks.size()); ++j) {
                if(S.boolean.value(intervals[tasks[j]].exist))
                    order.push_back(j);
            }
            
            std::sort(order.begin(), order.end(),
                      [&](const index_t a, const index_t b) {
                return S.numeric.lower(intervals[tasks[a]].start) < S.numeric.lower(intervals[tasks[b]].start);
            });
            
            std::cout << "resource " << i << ":";
            index_t n{Constant::NoIndex};
            for (auto o : order) {
                if(n != Constant::NoIndex) {
                    std::cout << " -> " << resource_transitions[i][order[n]][o] << " -> " ;
                    ++n;
                } else {
                    n = 0;
                }
                
                std::cout << " job" << o << ":" << prettyJob(intervals[tasks[o]], S, false);
            }
            std::cout << std::endl;
        }
        
        ++i;
    }
}



// implementation of a scheduling solver
int main(int argc, char *argv[]) {

  Options opt = tempo::parse(argc, argv);
  Solver<> S(opt);

  // an interval standing for the makespan of schedule
  auto schedule{S.newInterval(0, Constant::Infinity<int>, 0, 0, 0,
                              Constant::Infinity<int>)};
  //    auto s{S.newNumeric()};
  //    auto e{S.newNumeric()};
  //    auto schedule{S.newInterval(s,e)};

  std::vector<NoOverlapExpression<>> resources;
  std::vector<std::vector<size_t>> resource_tasks;
  std::vector<Interval<>> intervals;
  std::vector<int> weights;
  std::vector<std::vector<std::vector<int>>> resource_transitions;

  tsptw::parse(opt.instance_file, S, schedule, intervals, weights,
               resource_tasks, resource_transitions, true);

  index_t i{0};
  std::vector<Interval<int>> scope;
  for (auto &tasks : resource_tasks) {
    for (auto j : tasks) {
      scope.push_back(intervals[j]);
    }
    auto no_overlap{NoOverlap(schedule, scope, resource_transitions[i++])};
    resources.push_back(no_overlap);
    S.post(no_overlap);
    scope.clear();
  }

    std::vector<BooleanVar<int>> selection;
    for(auto i : intervals) {
        selection.push_back(i.exist);
        S.addToSearch(i.exist);
    }

  if (opt.print_mod) {
    std::cout << S << std::endl;
  }

  // search
  //  S.maximize(schedule.duration);
  S.maximize(Sum(selection, weights));

  if (opt.print_sol) {
    printJobs(S, intervals);
    printResources(S, intervals, resource_tasks, resource_transitions);
  }
}
