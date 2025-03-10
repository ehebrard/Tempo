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
#include "heuristics/warmstart.hpp"
#include "util/printing.hpp"

using namespace tempo;

// implementation of a scheduling solver
int main(int argc, char *argv[]) {

  Options opt = tempo::parse(argc, argv);
  Solver<> S(opt);

  auto makespan{S.newNumeric(0,Constant::Infinity<int>)};
  auto origin{S.newConstant(0)};
  auto schedule{S.between(origin, makespan)};

  std::vector<NoOverlapExpression<>> resources;
  std::vector<std::vector<size_t>> resource_tasks;
  std::vector<Interval<int>> intervals;
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
  

  auto obj{S.newNumeric()};

  
  S.post(obj == Sum(selection, weights));
    
    
    
  if (opt.print_mod) {
    std::cout << S << std::endl;
  }
  
  bool optimal{false};
  try {
      auto ub = Constant::Infinity<int>;
    heuristics::warmstartDisjunctive(S, schedule, intervals, ub);
  } catch (Failure<int> &f) {
    std::cout << "Greedy found optimal solution" << std::endl;
    optimal = true;
  }  
  
  if(!optimal){
    // search
    //  S.maximize(schedule.duration);
    std::cout << "Start maximisation" << std::endl;
    S.maximize(obj);
  }

  if (opt.print_sol) {
    printJobs(S, intervals);
    printResources(S, intervals, resource_tasks, resource_transitions);
  }
}
