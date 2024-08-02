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
#include "util/parsing/rcpsp.hpp"


using namespace tempo;






// implementation of a scheduling solver
int main(int argc, char *argv[]) {

  Options opt = tempo::parse(argc, argv);
  Solver<int> S(opt);
  seed(opt.seed);

  // an interval standing for the makespan of schedule
  auto schedule{S.newInterval(0, Constant::Infinity<int>, 0, 0, 0,
                              Constant::Infinity<int>)};

  // depending on the option "input-format", parse a disjunctive scheduling
  // instance, and collect resources and interval objects
  std::vector<CumulativeExpression<>> resources;
  std::vector<std::vector<size_t>> tasks_requirements;
    std::vector<std::vector<int>> task_demands;
    std::vector<int> resource_capacities;
  std::vector<Interval<>> intervals;


    rcpsp::parse(opt.instance_file, S, schedule, intervals, tasks_requirements, task_demands, resource_capacities);

    //      for(auto i : intervals) {
    //          std::cout << i.id() << ": " << i << std::endl;
    //      }

    std::vector<std::vector<Interval<int>>> resource_tasks(resource_capacities.size());
    std::vector<std::vector<NumericVar<int>>> resource_demands(resource_capacities.size());
    for(size_t j{0}; j<tasks_requirements.size(); ++j) {
        for(size_t k{0}; k<tasks_requirements[j].size(); ++k) {
            auto m{tasks_requirements[j][k]};
            auto d{task_demands[j][k]};
            resource_tasks[m].push_back(intervals[j]);
            resource_demands[m].push_back(S.newConstant(d));
        }
    }
  
    for(size_t k{0}; k<resource_capacities.size(); ++k) {
        NumericVar<int> capacity{S.newConstant(resource_capacities[k])};
        resources.push_back(Cumulative<int>(capacity, resource_tasks[k], resource_demands[k]));
        S.post(resources.back());
    }
 
  if (opt.print_mod) {
    std::cout << S << std::endl;
  }
//
//  auto optimal{false};
//  if (opt.greedy_runs > 0)
//    try {
//      warmstart(S, schedule, intervals /*, resources*/, ub);
//    } catch (Failure<int> &f) {
//      //            std::cout << " optimal solution found in a greedy run\n";
//      optimal = true;
//    }
//
//  // search
//  if (not optimal)

  // set a trivial (and the user-defined) upper bound
  int total_duration{0};
  for (auto &j : intervals) {
    if (j.maxDuration(S) == Constant::Infinity<int>) {
      total_duration = Constant::Infinity<int>;
      break;
    }
    total_duration += j.maxDuration(S);
  }
  auto ub{std::min(opt.ub, total_duration)};

  S.post(schedule.end.before(opt.ub));

  S.minimize(schedule.duration);
  //
  //  if (opt.print_sol) {
  //    printJobs(S, intervals);
  //    printResources(S, intervals, resource_tasks);
  //  }
}
