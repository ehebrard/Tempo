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
    
    if(est==lst)
        ss << est;
    else
        ss << est << "-" << lst;
    ss << ".." ;
    if(ect == lct)
        ss << ect;
    else
        ss << ect << "-" << lct ;
    ss << "]";
    
    if(dur_flag) {
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
        std::cout << "job" << ++i << ": " << prettyJob(task, S, true) << std::endl;
    }
}

template<typename T>
void printResources(const Solver<T>& S, const std::vector<Interval<T>>& intervals, const std::vector<std::vector<Interval<T>>>& resource_tasks) {
    int i{0};
    std::vector<int> jobmap(S.numeric.size(), -1);
    for (auto task : intervals) {
        jobmap[task.id()] = ++i;
    }
    i = 0;
    for(auto &tasks : resource_tasks) {
        ++i;
        if(tasks.size() > 1) {
            if(tasks.size() > 0) {
                auto jobs{tasks};
                std::sort(jobs.begin(), jobs.end(), [&](const Interval<int>& a, const Interval<int>& b) {return S.numeric.lower(a.start) < S.numeric.lower(b.start); });
                
                std::cout << "resource " << i << ":";
                for (auto task : jobs) {
                    std::cout << " " << jobmap[task.id()] << ":" << prettyJob(task, S, false);
                }
                std::cout << std::endl;
            }
        }
    }
}


// implementation of a scheduling solver
int main(int argc, char *argv[]) {

  Options opt = tempo::parse(argc, argv);
  Solver<> S(opt);

  // an interval standing for the makespan of schedule
  auto schedule{S.newInterval(0,Constant::Infinity<int>,0,0,0,Constant::Infinity<int>)};

  // depending on the option "input-format", parse a disjunctive scheduling
  // instance, and collect resources and interval objects
  //  std::vector<DisjunctiveResource<>> resources;
  std::vector<NoOverlapExpression<>> resources;
  std::vector<std::vector<Interval<>>> resource_tasks;
  std::vector<Interval<>> intervals;

  //    SchedulingModel<T> model;
    
  if (opt.input_format == "osp") {
    osp::parse(opt.instance_file, S, schedule, intervals, resource_tasks);
  } else if (opt.input_format == "jsp") {
    jsp::parse(opt.instance_file, S, schedule, intervals, resource_tasks);
  } else if (opt.input_format == "path") {
    path::parse(opt.instance_file, S, schedule, intervals, resource_tasks);
  }
  //    else if (opt.input_format == "tsptw") {
  //        tsptw::parse(opt.instance_file, S, schedule, intervals,
  //        resources);
  //    }
  else if (opt.input_format == "jstl") {
    jstl::parse(opt.instance_file, S, schedule, intervals, resource_tasks);
  }

//  std::vector<NoOverlapExpression<>> res;
  for (auto &tasks : resource_tasks) {
    auto no_overlap{NoOverlap(schedule, tasks)};
    resources.push_back(no_overlap);
    S.post(no_overlap);
  }

  // set a trivial (and the user-defined) upper bound
  auto trivial_ub{0};
  for (auto &j : intervals) {
    if (j.maxDuration(S) == Constant::Infinity<int>) {
      trivial_ub = Constant::Infinity<int>;
      break;
    }
    trivial_ub += j.maxDuration(S);
  }
  auto ub{std::min(opt.ub, trivial_ub)};

  //    S.post(schedule.end <= ub);
  S.post(schedule.end.before(ub));

  if (opt.print_mod)
    std::cout << S << std::endl;

  //    for(auto i : intervals) {
  //
  //    }

  //  warmstart(S, schedule, intervals, resources, ub);

  // search
  S.minimize(schedule.duration);

    
  if (opt.print_sol) {
      printJobs(S, intervals);
      printResources(S, intervals, resource_tasks);
  }
}
