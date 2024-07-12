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
  }
  //    else if (opt.input_format == "tsptw") {
  //        tsptw::parse(opt.instance_file, S, schedule, Intervals,
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
  for (auto &j : intervals)
    trivial_ub += j.minDuration(S);
  auto ub{std::min(opt.ub, trivial_ub)};

  //    S.post(schedule.end <= ub);
  S.post(schedule.end.before(ub));

  std::cout << S << std::endl;
  S.propagate();
  std::cout << S << std::endl;

  //    for(auto i : intervals) {
  //
  //    }

  warmstart(S, schedule, intervals, resources, ub);

  // search
  S.minimize(schedule.duration);
}
