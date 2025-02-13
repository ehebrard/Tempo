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
#include "util/parsing/airbus.hpp"
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
    
    opt.edge_finding = false;
    opt.transitivity = false;
    
  seed(opt.seed);
    
    
//    SchedulingInstance<int> data;
    Solver<int> solver(opt);
    
    airbus::parse(opt.instance_file, solver);
    
    if (opt.print_mod) {
        std::cout << solver << std::endl;
    }
    
        
//    solver.minimize(model.getScheduleInterval().duration);
    
    auto sat{solver.satisfiable()};
    
    std::cout << sat << std::endl;

}
