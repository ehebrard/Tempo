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
#include "util/parsing/jsp.hpp"
#include "util/parsing/jstl.hpp"
#include "util/parsing/osp.hpp"

using namespace tempo;

// implementation of a scheduling solver
int main(int argc, char *argv[]) {

  Options opt = tempo::parse(argc, argv);
  Solver<> S(opt);

  // an interval standing for the makespan of schedule
  auto schedule{S.newInterval()};

  // depending on the option "input-format", parse a disjunctive scheduling
  // instance, and collect resources and interval objects
  std::vector<DisjunctiveResource<>> resources;
  std::vector<Interval<>> Intervals;
  if (opt.input_format == "osp") {
    osp::parse(opt.instance_file, S, schedule, Intervals, resources);
  } else if (opt.input_format == "jsp") {
    jsp::parse(opt.instance_file, S, schedule, Intervals, resources);
  }
  //    else if (opt.input_format == "tsptw") {
  //        tsptw::parse(opt.instance_file, S, schedule, Intervals,
  //        resources);
  //    }
  else if (opt.input_format == "jstl") {
    jstl::parse(opt.instance_file, S, schedule, Intervals, resources);
  }

  // create disjunct variables for each resource and post the propagators
  std::vector<BooleanVar<>> X;
  for (auto &R : resources) {
    auto s{X.size()};
    R.createOrderVariables(S, X);

    // post the propagator for the edge-finding rule
    if (opt.edge_finding) {
      S.postEdgeFinding(schedule, R.begin(), R.end(), X.begin() + s);
    }

    // post the propagator for the precedence reasoning rule
    if (opt.transitivity) {
      S.postTransitivity(schedule, R.begin(), R.end(), X.begin() + s);
    }
  }

  // notify the solver to assign a value to all disjunct variables
  for (auto x : X)
    S.addToSearch(x);

  // the objective
  //  MakespanObjective<int> duration(schedule, S);
  MinimizationObjective<int, TemporalVar<int>> duration(schedule.end);

  // set a trivial (and the user-defined) upper bound
  auto trivial_ub{0};
  for (auto &j : Intervals)
    trivial_ub += j.minDuration();
  auto ub{std::min(opt.ub, trivial_ub)};
  S.set(schedule.end.before(ub));

  // try to get a better ub with an initial upper bound insertion heuristic
  Greedy greedy_insertion(S);
  greedy_insertion.addIntervals(Intervals);
  size_t k{0};
  for (auto &R : resources) {
    auto n = k + R.size() * (R.size() - 1) / 2;
    greedy_insertion.addResource(X.begin() + k, X.begin() + n);
    k = n;
  }

  // the insertion heuristic is randomized so multiple runs can be useful
  S.initializeSearch();
  S.propagate();
  for (auto i{0}; i < opt.greedy_runs; ++i) {
    auto st{S.saveState()};
    auto sat{greedy_insertion.runEarliestStart()};
    if (sat) {
      if (schedule.getEarliestEnd(S) <= ub) {
        ub = schedule.getEarliestEnd(S) - 1;
        std::cout << std::setw(10) << (ub + 1);
        S.displayProgress(std::cout);
      }
    }
    S.restoreState(st);
  }

  // set the ub (again)
  S.set(schedule.end.before(ub));

  // search
  S.optimize(duration);
}
