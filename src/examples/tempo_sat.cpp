/************************************************
 * Tempo sat solver
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
#include "util/parsing/dimacs.hpp"

using namespace tempo;

// implementation of SAT solver
int main(int argc, char *argv[]) {

  Options opt = tempo::parse(argc, argv);
  Solver<> S(opt);

  // parse an input file in dimacs cnf format and collect the variables
  std::vector<BooleanVar<>> X;
  dimacs::parse(opt.instance_file, S, X);

  // notify the solver to assign a value to all variables
  for (auto x : X)
    S.addToSearch(x);

  // search
  auto sat{S.satisfiable()};

  // output
  std::cout << (sat ? "SAT" : "UNSAT") << " #fails = " << S.num_fails
            << std::endl;
}
