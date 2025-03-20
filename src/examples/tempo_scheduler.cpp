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
#include "heuristics/Static.hpp"
#include "heuristics/warmstart.hpp"
#include "util/parsing/jsp.hpp"
#include "util/parsing/jstl.hpp"
#include "util/parsing/osp.hpp"
#include "util/parsing/path.hpp"
#include "util/parsing/tsptw.hpp"
#include "util/parsing/jssdst.hpp"
#include "util/printing.hpp"
//#include "util/parsing/airbus.hpp"
#include "helpers/cli.hpp"
#include "util/Profiler.hpp"
#include "helpers/shell.hpp"
#include "helpers/git_sha.hpp"
#include "Solution.hpp"

using namespace tempo;


// implementation of a scheduling solver
int main(int argc, char *argv[]) {
        
  auto parser = tempo::getBaseParser();
  bool profileHeuristic;
  cli::detail::configureParser(parser, cli::SwitchSpec("heuristic-profiling", "activate heuristic profiling",
                                                       profileHeuristic, false));
    
    std::string save_solution{""};
    parser.getCmdLine().add<TCLAP::ValueArg<std::string>>(save_solution, "", "save-solution", "save solution to file", false, "", "string");
    
    std::string load_solution{""};
    parser.getCmdLine().add<TCLAP::ValueArg<std::string>>(load_solution, "", "load-solution", "load solution from file", false, "", "string");
    
  parser.parse(argc, argv);
  Options opt = parser.getOptions();
    
  if(opt.print_par)
    std::cout << opt << std::endl;
    
  seed(opt.seed);
    
  Solver<> S(opt);

  // an interval standing for the makespan of schedule
//  auto schedule{S.newInterval(0, Constant::Infinity<int>, 0, 0, 0,
//                              Constant::Infinity<int>)};
    auto makespan{S.newNumeric(0,Constant::Infinity<int>)};
    auto origin{S.newConstant(0)};
    auto schedule{S.between(origin, makespan)};

  // depending on the option "input-format", parse a disjunctive scheduling
  // instance, and collect resources and interval objects
  std::vector<NoOverlapExpression<>> resources;
  
    std::vector<std::vector<size_t>> resource_tasks;
    
    
  std::vector<Interval<>> intervals;
  std::vector<int> weights;
  std::vector<std::vector<std::vector<int>>> resource_transitions;

  if (opt.input_format == "osp") {
    osp::parse(opt.instance_file, S, schedule, intervals, resource_tasks);
  } else if (opt.input_format == "jsp") {
    jsp::parse(opt.instance_file, S, schedule, intervals, resource_tasks);
  } else if (opt.input_format == "path") {
    path::parse(opt.instance_file, S, schedule, intervals, resource_tasks);
  } else if (opt.input_format == "tsptw") {
    tsptw::parse(opt.instance_file, S, schedule, intervals, weights,
                 resource_tasks, resource_transitions);
  } else if (opt.input_format == "jstl") {
    jstl::parse(opt.instance_file, S, schedule, intervals, resource_tasks);
  } else if (opt.input_format == "jssdst") {
      jssdst::parse(opt.instance_file, S, schedule, intervals, resource_tasks, resource_transitions);
    }
//  else if (opt.input_format == "airbus") {
//      airbus::parse(opt.instance_file, S, schedule, intervals, resource_tasks);
//    }

  //    for(auto i : intervals) {
  //        std::cout << i.id() << ": " << i << std::endl;
  //    }
  //    exit(1);
    
    
//    std::vector<size_t> task2machine(intervals.size());
//    
//    size_t i{0};
//    for(size_t r{0}; r<resource_tasks.size(); ++r) {
//        for(size_t k{0}; k<resource_tasks[r].size(); ++k) {
//            task2machine[i++] = r;
//        }
//    }
    
    
//    std::cout << S << std::endl;

  resource_transitions.resize(resource_tasks.size());

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

//  // set a trivial (and the user-defined) upper bound
//  int total_duration{0};
//  int max_start{0};
//  for (auto &j : intervals) {
//    if (j.maxDuration(S) == Constant::Infinity<int>) {
//      total_duration = Constant::Infinity<int>;
//      break;
//    }
//    total_duration += j.maxDuration(S);
//
//    if (j.start.min(S) > max_start)
//      max_start = j.start.min(S);
//  }
//  auto ub{std::min(opt.ub, max_start + total_duration)};
//
//  //    S.post(schedule.end <= ub);

  if (opt.ub != Constant::Infinity<int>) {
    S.post(schedule.end.before(opt.ub));
    //        S.set(schedule.end.before(schedule.start, -opt.ub));
    //        S.propagate();

    //        std::cout << (schedule.end.before(schedule.start, -opt.ub)) <<
    //        std::endl;
  }

  // HACK!!: add an edge from origin to end to allow propagation of
  // shortest path
//  S.set({0, 1, S.numeric.upper(1)});

  int num_restart{0};
  SubscriberHandle handlerToken;
  FullTransitivity<int>* primal{NULL};
  if (opt.full_transitivity or opt.primal_boost) {
    primal = S.postFullTransitivity(resources.begin(), resources.end());
    if (opt.primal_boost)
      handlerToken = S.SearchRestarted.subscribe_handled([&](bool) {
            ++num_restart;

            if (S.num_solutions > 2 or S.num_fails >= 250) {
              S.relax(primal);
              handlerToken.unregister();
            }
          });
  }

  if (opt.print_mod) {
    std::cout << S << std::endl;
  }

  util::Profiler profiler;
  if (profileHeuristic) {
      using namespace tempo::heuristics;
      using T = decltype(S)::Time;
      util::ProfiledHeuristic<VariableHeuristic<T>> varBranching(profiler, make_variable_heuristic(S));
      util::ProfiledHeuristic<ValueHeuristic<T>> valBranching(profiler, make_value_heuristic(S));
      S.setBranchingHeuristic(make_compound_heuristic(std::move(varBranching), std::move(valBranching)));
  }

  auto optimal{false};
  if (opt.greedy_runs > 0) {

    //        std::vector<Interval<int>> by_resource;
    //        for(auto &tasks : resource_tasks) {
    //            for(auto j : tasks) {
    //                by_resource.push_back(intervals[j]);
    //            }
    //        }

    try {
      auto ub{Constant::Infinity<int>};
//      heuristics::warmstartDisjunctive(S, schedule, intervals, ub);
        
//        heuristics::warmstartEstSlack(S, schedule, intervals, ub);
        heuristics::warmstartSlackEst(S, schedule, intervals, ub);
//        exit(1);
      //            warmstart(S, schedule, by_resource, ub);
    } catch (Failure<int> &f) {
      //            std::cout << " optimal solution found in a greedy run\n";
      optimal = true;
    }
  }


  // search
  if (not optimal) {
    S.minimize(schedule.duration);
  }

  if (S.numeric.hasSolution()) {
      std::cout << "-- makespan " << S.numeric.solutionLower(schedule.duration) << std::endl;
  }

  if (opt.print_sol) {
    printJobs(S, intervals);
    printResources(S, intervals, resource_tasks);
  }

  profiler.printAll<std::chrono::microseconds>(std::cout);
  std::cout << "-- date: " << shell::getTimeStamp() << std::endl;
  std::cout << "-- commit: " << GitSha << std::endl;
    
    if(save_solution != "") {
        std::ofstream cl_file(save_solution, std::ofstream::out);
        Solution<int> best(S);
        
        std::cout << "save " << best << std::endl;
        cl_file << best << std::endl;
    }

}
