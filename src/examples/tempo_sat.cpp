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


#include <iostream>
#include <vector>


#include "Solver.hpp"
#include "util/parsing/dimacs.hpp"

using namespace tempo;



int main(int argc, char *argv[]) {

  Options opt = tempo::parse(argc, argv);

    Solver<> S(opt);

    std::vector<BooleanVar<>> X;

    dimacs::parse(opt.instance_file, S, X);

    for (auto x : X)
      S.addToSearch(x);

    auto sat{S.satisfiable()};

    std::cout << (sat ? "SAT" : "UNSAT") << " #fails = " << S.num_fails << std::endl;
}
