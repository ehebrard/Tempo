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
#include <filesystem>

#include "util/Options.hpp"
#include "util/parsing/format.hpp"
#include "util/parsing/jsp.hpp"
#include "util/parsing/osp.hpp"
#include "util/parsing/tsptw.hpp"

#include <ilcp/cp.h>

using namespace tempo;

int main(int argc, char *argv[]) {
  Options opt = tempo::parse(argc, argv);

  ProblemInstance data;

  IloEnv env;
  IloModel model(env);

  if (opt.input_format == "osp") {
    data = osp::read_instance(opt.instance_file);
  } else if (opt.input_format == "jsp") {
    data = jsp::read_instance(opt.instance_file);
  } else if (opt.input_format == "tsptw") {
    data = tsptw::read_instance(opt.instance_file);
  }

  IloIntExprArray ends(env);
  for (auto &resource : data.resources) {
    IloIntervalVarArray scope;
    for (auto ti : resource) {
      IloIntervalVar tx(env, data.duration[ti]);
      scope.add(tx);
      ends.add(IloEndOf(tx));
    }
    model.add(IloNoOverlap(env, scope));
  }

  IloObjective objective = IloMinimize(env, IloMax(ends));
  model.add(objective);

  IloCP cp(model);
  cp.setParameter(IloCP::Workers, 1);
  cp.setParameter(IloCP::FailLimit, failLimit);
  cp.setParameter(IloCP::LogVerbosity, IloCP::Terse);
  cp.out() << "Instance \t: " << filename << std::endl;
  if (cp.solve()) {
    cp.out() << "Makespan \t: " << cp.getObjValue() << std::endl;
  } else {
    cp.out() << "No solution found." << std::endl;
  }

  env.end();
  return 0;
}
