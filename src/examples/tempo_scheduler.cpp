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
#include "util/parsing/jsp.hpp"
#include "util/parsing/jstl.hpp"
#include "util/parsing/osp.hpp"
#include "util/parsing/path.hpp"
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
string prettyJob(const Interval<T>& task, const Solver<T>& S, const bool dur_flag) {
  std::stringstream ss;

  auto est{S.numeric.lower(task.start)};
  auto lst{S.numeric.upper(task.start)};
  auto ect{S.numeric.lower(task.end)};
  auto lct{S.numeric.upper(task.end)};

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
void printResources(const Solver<T>& S, const std::vector<Interval<T>>& intervals, std::vector<std::vector<size_t>>& resource_tasks) {
  int i{0};
    
  for (auto &tasks_idx : resource_tasks) {
      ++i;
    if (tasks_idx.size() > 1) {
        
        std::vector<index_t> order;
        for(auto j : tasks_idx) {
            auto job{intervals[j]};
            if(S.boolean.value(job.exist))
                order.push_back(j);
        }
        
        std::sort(order.begin(), order.end(),
                  [&](const index_t a, const index_t b) {
                    return S.numeric.lower(intervals[a].start) < S.numeric.lower(intervals[b].start);
                  });
 
        std::cout << "resource " << i << ":";
        for (auto o : order) {
          std::cout << " job" << (o+1) << ":" << prettyJob(intervals[o], S, false);
        }
        std::cout << std::endl;
      }
  }
}


// implementation of a scheduling solver
int main(int argc, char *argv[]) {

  Options opt = tempo::parse(argc, argv);
  Solver<> S(opt);

  // an interval standing for the makespan of schedule
  auto schedule{S.newInterval(0, Constant::Infinity<int>, 0, 0, 0,
                              Constant::Infinity<int>)};

  // depending on the option "input-format", parse a disjunctive scheduling
  // instance, and collect resources and interval objects
  std::vector<NoOverlapExpression<>> resources;
    std::vector<std::vector<size_t>> resource_tasks;
  std::vector<Interval<>> intervals;
  std::vector<std::vector<std::vector<int>>> resource_transitions;

  if (opt.input_format == "osp") {
    osp::parse(opt.instance_file, S, schedule, intervals, resource_tasks);
  } else if (opt.input_format == "jsp") {
    jsp::parse(opt.instance_file, S, schedule, intervals, resource_tasks);
  } else if (opt.input_format == "path") {
    path::parse(opt.instance_file, S, schedule, intervals, resource_tasks);
  } else if (opt.input_format == "tsptw") {
    tsptw::parse(opt.instance_file, S, schedule, intervals, resource_tasks,
                 resource_transitions);
  } else if (opt.input_format == "jstl") {
    jstl::parse(opt.instance_file, S, schedule, intervals, resource_tasks);
  }

  resource_transitions.resize(resource_tasks.size());

    
        index_t i{0};
        std::vector<Interval<int>> scope;
        for (auto &tasks : resource_tasks) {
            for(auto j : tasks) {
                scope.push_back(intervals[j]);
            }
          auto no_overlap{NoOverlap(schedule, scope, resource_transitions[i++])};
          resources.push_back(no_overlap);
          S.post(no_overlap);
            scope.clear();
        }

  // set a trivial (and the user-defined) upper bound
  int total_duration{0};
  int max_start{0};
  for (auto &j : intervals) {
    if (j.maxDuration(S) == Constant::Infinity<int>) {
      total_duration = Constant::Infinity<int>;
      break;
    }
    total_duration += j.maxDuration(S);

    if (j.start.min(S) > max_start)
      max_start = j.start.min(S);
  }
  auto ub{std::min(opt.ub, max_start + total_duration)};

  //    S.post(schedule.end <= ub);
  S.post(schedule.end.before(ub));

  if (opt.print_mod) {
    std::cout << S << std::endl;
  }

    if(opt.greedy_runs > 0)
  warmstart(S, schedule, intervals, resources, ub);

  // search
  S.minimize(schedule.duration);

  if (opt.print_sol) {
    printJobs(S, intervals);
    printResources(S, intervals, resource_tasks);
  }
}
