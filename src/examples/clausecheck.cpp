/*************************************************************************
minicsp

Copyright 2010--2011 George Katsirelos

Minicsp is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

Minicsp is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with minicsp.  If not, see <http://www.gnu.org/licenses/>.

*************************************************************************/


#include <chrono>
#include <iostream>
#include <vector>
#include <filesystem>
#include <optional>

#include "Solver.hpp"
#include "heuristics/Greedy.hpp"
#include "util/parsing/jsp.hpp"
#include "util/parsing/jstl.hpp"
#include "util/parsing/osp.hpp"
#include "util/parsing/path.hpp"
#include "util/parsing/tsptw.hpp"


using namespace tempo;


int main(int argc, char *argv[]) {
  Options opt = tempo::parse(argc, argv);
    
    size_t num_fails{0};
    size_t num_trivial{0};
    size_t num_search{0};

  seed(opt.seed);
    
    std::ifstream cl_file(opt.dbg_file, std::ifstream::in);

    opt.dbg_file = "";
    opt.learning = false;
    
    
    
    
    Solver<int> S(opt);
    
    // an interval standing for the makespan of schedule
    auto schedule{S.newInterval(0, Constant::Infinity<int>, 0, 0, 0,
                                Constant::Infinity<int>)};

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
    
    
//    
//    Solver<> S(opt);
//
//      auto schedule{S.newJob()};
//      std::vector<DisjunctiveResource<>> resources;
//      std::vector<Job<>> jobs;
//      
//      osp::parse(opt.instance_file, S, schedule, jobs, resources);
//      
//      std::vector<BooleanVar<>> vars;
//      for(auto &R : resources) {
//          R.createOrderVariables(S, vars);
//      }
//          
//      for (auto x : vars)
//        S.addToSearch(x);
//    
//    S.initializeSearch();
//    
//    
//    S.propagate();
//    
////    std::cout << S << std::endl;
 
    var_t x, y;
  int t, n, d;

  std::vector<var_t> X;
  std::vector<var_t> Y;
  std::vector<int> D;

  int line{0};
  do {

    X.clear();
    Y.clear();
    D.clear();

//      std::cout << "save:\n" << S << std::endl;
      
//      std::cout << "\n\nclause " << (line+1) << std::endl;
      
    S.saveState();

    cl_file >> t;
    cl_file >> n;
    if (not cl_file.good())
      break;

    bool trivially_unsat{false};

    for (auto i{0}; i < n; ++i) {
      cl_file >> x;
      cl_file >> y;
      cl_file >> d;

      X.push_back(x);
      Y.push_back(y);
      D.push_back(d);

      //            std::cout << "add " << prettyEvent(y) << " - " <<
      //            prettyEvent(x) << " <= " << d << std::endl;

      try {
          DistanceConstraint<int> c{x, y, d};
        S.post(c);
      } catch (Failure<int> &f) {
        trivially_unsat = true;
      }
    }
      
//      std::cout << S << std::endl;
//      exit(1);

    //        std::cout << "\nsolve:\n" << S << std::endl;

    if (not trivially_unsat) {
      //            std::cout << "ok (trivial)\n";
      //        } else {

      //            std::cout << S << std::endl;
      //            exit(1);

      bool need_search{true};
      try {
        S.propagate();
      } catch (Failure<int> &f) {
        need_search = false;
      }

      if (need_search) {

//                  std::cout << S << std::endl;

        auto nf{S.num_fails};
//        S.search();

        if (S.satisfiable()) {
          std::cout << "cl " << line << " ("
                    << (t == 0 ? "minimized" : (t == 1 ? "reason" : "uip"))
                    << "): ";
          std::cout << "bug!\n";

          for (size_t i{0}; i < X.size(); ++i) {
              std::cout << "> " << DistanceConstraint<int>(X[i], Y[i], D[i]) << std::endl;
          }

          if (opt.print_sol) {
            for (auto i : intervals) {
              std::cout << "x" << i.start.id() << ": "
                        << S.numeric.lower(i.start) << ".."
                        << S.numeric.lower(i.end) << " ("
                        << S.numeric.lower(i.duration) << ")\n";
            }
          }

          exit(1);
        }
        //                else {
        //                    std::cout << "UNSAT\n";
        //                }

        ++num_search;
        num_fails += (S.num_fails - nf);
      }
      //            else {
      //                std::cout << "ok\n";
      //            }
    } else {
      ++num_trivial;
    }
    ++line;

    //        std::cout << S.env.level() << " --> 0"

    S.restoreState(0);
//      S.undo();
      
//      std::cout << S.numLiteral() << std::endl;

//            std::cout << "\nrestored:\n" << S << std::endl;

    //    if ((line % 100) == 0)
    std::cout << line << ": " << num_trivial << " trivial, "
              << (line - num_trivial - num_search) << " easy, " << num_search
              << " hard (" << num_fails / num_search << ")" << std::endl;

    //        if(line > 100)
    //            exit(1);
    } while(true);
 
}
