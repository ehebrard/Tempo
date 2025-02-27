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
#include "util/parsing/psplib.hpp"


using namespace tempo;


int main(int argc, char *argv[]) {
  Options opt = tempo::parse(argc, argv);
    
    size_t num_fails{0};
    size_t num_trivial{0};
    size_t num_search{0};

  seed(opt.seed);
    
    std::ifstream cl_file(opt.dbg_file, std::ifstream::in);

    opt.dbg_file = "";
//    opt.learning = false;
    opt.overlap_finding = false;
    opt.tt_edge_finding = false;
    opt.edge_finding = false;
    
    
    
    Solver<int> S(opt);
    
    // an interval standing for the makespan of schedule
    auto makespan{S.newNumeric(0,Constant::Infinity<int>)};
    auto origin{S.newConstant(0)};
    auto schedule{S.between(origin, makespan)};
//    auto schedule{S.newInterval(0, Constant::Infinity<int>, 0, 0, 0,
//                                Constant::Infinity<int>)};

    // depending on the option "input-format", parse a disjunctive scheduling
    // instance, and collect resources and interval objects
    std::vector<CumulativeExpression<>> resources;
    std::vector<std::vector<size_t>> tasks_requirements;
      std::vector<std::vector<int>> task_demands;
      std::vector<int> resource_capacities;
    std::vector<Interval<>> intervals;
    std::vector<std::pair<int, int>> precedences;
    std::vector<std::vector<int>> graph;

    psplib::parse(opt.instance_file, S, schedule, intervals, tasks_requirements,
                 task_demands, resource_capacities, precedences, graph);

    std::vector<std::vector<Interval<int>>> resource_tasks(
        resource_capacities.size());
    std::vector<std::vector<NumericVar<int>>> resource_demands(
        resource_capacities.size());
    for (size_t j{0}; j < tasks_requirements.size(); ++j) {
      for (size_t k{0}; k < tasks_requirements[j].size(); ++k) {
        auto m{tasks_requirements[j][k]};
        auto d{task_demands[j][k]};
        resource_tasks[m].push_back(intervals[j]);
        resource_demands[m].push_back(S.newConstant(d));
      }
    }

      for(size_t k{0}; k<resource_capacities.size(); ++k) {
          NumericVar<int> capacity{S.newConstant(resource_capacities[k])};
          resources.push_back(Cumulative<int>(
              schedule, capacity, resource_tasks[k], resource_demands[k]));
          S.post(resources.back());
      }
    
  
 
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

   
      try {
          DistanceConstraint<int> c{x, y, d};
        S.post(c);
      } catch (Failure<int> &f) {
        trivially_unsat = true;
      }
    }
 
    if (not trivially_unsat) {
  
      bool need_search{true};
      try {
        S.propagate();
      } catch (Failure<int> &f) {
        need_search = false;
      }

      if (need_search) {
        auto nf{S.num_fails};
        if (S.satisfiable()) {
          std::cout << "cl " << line << " ("
            << (t == 1 ? "minimized" : (t == 2 ? "reason" : "uip"))
//            << (t ? "expl" : "cut")
            << "): ";
          std::cout << "bug!\n";

          for (size_t i{0}; i < X.size(); ++i) {
              std::cout << "> " << DistanceConstraint<int>(X[i], Y[i], D[i]) << std::endl;
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
