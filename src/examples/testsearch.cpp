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

//#include "Objective.hpp"
#include "Scheduler.hpp"
#include "util/parsing/format.hpp"
#include "util/parsing/jsp.hpp"
#include "util/parsing/jstl.hpp"
#include "util/parsing/osp.hpp"
#include "util/parsing/tsptw.hpp"

using namespace tempo;

void load(std::string &solfile_name, std::vector<bool> &solution) {
  //    std::cout << "load\n";

  std::ifstream solfile(solfile_name.c_str(), std::ifstream::in);
  size_t n;
  bool val;
  solfile >> n;

  //    std::cout << n << std::endl;

  solution.resize(n);

  for (size_t i{0}; i < n; ++i) {
    solfile >> val;

    //        std::cout << val << std::endl;

    solution[i] = val;
  }
}

int main(int argc, char *argv[]) {

  auto start = std::chrono::system_clock::now();
  Options opt = tempo::parse(argc, argv);

  seed(opt.seed);

  ProblemInstance data;

  int ub{opt.ub};

//    std::cout << ub << " <> " << INFTY << std::endl;
    
  if (opt.input_format == "osp") {
    data = osp::read_instance(opt.instance_file);
    if (ub == INFTY)
      ub = osp::getUb<int>(data);
    else {
        std::cout << ub << " <> " << INFTY << std::endl;
        exit(1);
    }
  }
  else if (opt.input_format == "jsp") {
    data = jsp::read_instance(opt.instance_file);
    if (ub == INFTY)
      ub = jsp::getUb<int>(data);
  } else if (opt.input_format == "tsptw") {
    data = tsptw::read_instance(opt.instance_file);
  } else if (opt.input_format == "jstl") {
    data = jstl::read_instance(opt.instance_file);
    if (ub == INFTY)
      ub = jstl::getUb<int>(data);
  }

  if (opt.print_ins) {
    std::cout << data.durations.size() << " tasks:";
    for (auto d : data.durations) {
      std::cout << " " << d;
    }
    std::cout << std::endl;

    std::cout << data.resources.size() << " resources:";
    for (auto R : data.resources) {
      std::cout << " (";
      for (auto x : R) {
        std::cout << " " << x;
      }
      std::cout << ")";
    }
    std::cout << std::endl << "ub = " << ub << std::endl;
  }

  Scheduler<int> S(opt);
  std::vector<task> t;
  for (auto d : data.durations) {
    t.push_back(S.newTask(d, d));
  }

  //  S.setUpperBound(ub);

  for (auto [x, y, k] : data.constraints) {

    S.newPrecedence(x, y, k);
  }

  //  std::cout << data.resources.size() << std::endl;
  //  std::cout << data.resources[0].size() << std::endl;
  //  //    exit(1);

  std::vector<var> scope;
  std::vector<var> tasks;
  for (auto &job : data.resources) {
    for (size_t i{0}; i < job.size(); ++i) {
      tasks.push_back(t[job[i]]);
      for (size_t j{i + 1}; j < job.size(); ++j) {
        scope.push_back(S.newVariable(
            {START(t[job[i]]), END(t[job[j]]), -job.getTransitionTime(j, i)},
            {START(t[job[j]]), END(t[job[i]]), -job.getTransitionTime(i, j)}));
      }
    }

    //           for (auto ti{job.begin()}; ti!=job.end(); ++ti) {
    //          for (auto tj{ti+1}; tj!=job.end(); ++tj) {
    //            scope.push_back(S.newVariable({START(*ti), END(*tj), 0},
    //                                          {START(*tj), END(*ti), 0}));
    //          }
    //      }
    if (opt.edge_finding) {
      S.postEdgeFinding(tasks.begin(), tasks.end(), scope.begin(), scope.end());
    }
      if (opt.transitivity) {
        S.postTransitivityReasoning(tasks.begin(), tasks.end(), scope.begin(),
                                    scope.end());
      }
      tasks.clear();
      scope.clear();
  }

#ifdef DBG_SOL
  if (opt.dbg_file != "") {
    std::vector<bool> dbg_sol;
    load(opt.dbg_file, dbg_sol);
    S.load(dbg_sol);
    //    std::cout << "load";
    //    for (size_t i{0}; i < dbg_sol.size(); ++i)
    //      std::cout << " " << dbg_sol[i];
    //    std::cout << std::endl;
  }
#endif

  //    std::cout << S << std::endl;

  //  S.search();



  int sol;
  //  if (opt.input_format == "tsptw") {
  //    //        PathLength<int> tour(S, data.resources[0]);
  //    //        S.minimize(tour);
  //    //        sol = tour.upperBound();
  //
  //      if (opt.print_mod) {
  //        S.display(std::cout, true, true, true, true, false, false, true);
  //        std::cout << S << std::endl;
  //      }
  //
  //    sol = S.satisfiable();
  //      std::cout << (sol ? "SAT" : "UNSAT") << std::endl;
  //  } else {
  Makespan<int> makespan(S, ub);

  if (opt.print_mod) {
    S.display(std::cout, true, true, true, true, false, false, true);
    std::cout << S << std::endl;
  }

    //      makespan.setUpperBound(ub);
    if (opt.dichotomy)
      S.optimize_dichotomy(makespan);
    else
      S.optimize(makespan);
    
    sol = makespan.primalBound();
    //  }

    if (opt.verbosity > tempo::Options::SILENT) {
      auto end = std::chrono::system_clock::now();
      auto elapsed =
          std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
      std::cout << "Objective value = " << sol
                << " Execution time = " << elapsed.count() << " ms"
                << std::endl;
    }

  //  auto sol = S.getMakespan();
  //  auto sol = makespan.upperBound();
  if (KillHandler::instance().signalReceived()) {
    std::cout << "execution aborted, best solution found " << sol << std::endl;
  } else if (data.optimalSolution != -1 and sol != data.optimalSolution) {
    std::cerr << "Suboptimal solution! " << sol << " vs optimal "
              << data.optimalSolution << std::endl;
    std::exit(1);
  }

  //    vector<task>  the_tasks;
  //
  //    vector<vector<int>> order;
  //
  //    //order[i][j] -> indice de newVariable()
  //
  //    sort(the_tasks.begin(), the_tasks.end(), [](const task a, const task b)
  //    {
  //        return before[oder[a][b]];
  //    };
  //         );

  if (opt.print_sol) {
    auto before{S.getSolution()};
    std::cout << before.size();
    for (size_t i{0}; i < before.size(); ++i)
      std::cout << " " << before[i];
    std::cout << std::endl;
  }
}
