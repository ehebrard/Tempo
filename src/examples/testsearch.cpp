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

#include "Scheduler.hpp"
#include "util/parsing/format.hpp"
#include "util/parsing/jsp.hpp"
#include "util/parsing/osp.hpp"


using namespace tempo;


int main(int argc, char *argv[]) {
  auto start = std::chrono::system_clock::now();
  Options opt = tempo::parse(argc, argv);

  seed(opt.seed);

  ProblemInstance data;

  int ub;

  if (opt.input_format == "osp") {
    data = osp::read_instance(opt.instance_file);
    ub = osp::getUb(data.durations, data.resources);
  }
  else if (opt.input_format == "jsp") {
    data = jsp::read_instance(opt.instance_file);
    ub = jsp::getUb(data.durations, data.resources);
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
  for (auto d : data.durations) {
    S.newTask(
              //"t" + std::to_string(S.numTask()),
              d, d);
  }


  S.setUpperBound(ub);



    
  for (auto [x, y, k] : data.constraints) {

      assert(k==0);

    S.newPrecedence(x, y, k);

  }

  std::vector<var> scope;
  for (auto &job : data.resources) {
      for (auto ti{job.begin()}; ti!=job.end(); ++ti) {
          for (auto tj{ti+1}; tj!=job.end(); ++tj) {
            scope.push_back(S.newVariable({START(*ti), END(*tj), 0},
                                          {START(*tj), END(*ti), 0}));
          }
      }
      if (opt.edge_finding) {
          for(auto j : job) {
              std::cout << " t" << j ;
          }
          std::cout << std::endl;
          for(auto v : scope) {
              std::cout << " x" << v ;
          }
          std::cout << std::endl;
          
        S.postEdgeFinding(job.begin(), job.end(), scope.begin(), scope.end());
      }
      scope.clear();
  }
    
//    std::cout << S << std::endl;
    
  S.search();
    
  auto end = std::chrono::system_clock::now();
  auto elapsed =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "Temps d'execution = " << elapsed.count() << " ms" << std::endl;

  auto sol = S.getMakespan();
  if (KillHandler::instance().signalReceived()) {
    std::cout << "execution aborted, best solution found " << sol << std::endl;
  } else if (data.optimalSolution != -1 and sol != data.optimalSolution) {
    std::cerr << "Suboptimal solution! " << sol << " vs optimal "
              << data.optimalSolution << std::endl;
    std::exit(1);
  }
}
