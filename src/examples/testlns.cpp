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
#include <utility>
#include <string>

#include "Solver.hpp"
#include "heuristics/Greedy.hpp"
#include "util/parsing/jsp.hpp"
#include "util/parsing/jstl.hpp"
#include "util/parsing/osp.hpp"
#include "util/parsing/path.hpp"
#include "util/parsing/tsptw.hpp"
#include "util/parsing/jssdst.hpp"
#include "util/printing.hpp"
#include "helpers/cli.hpp"
#include "util/Profiler.hpp"
#include "helpers/shell.hpp"
#include "helpers/git_sha.hpp"
#include "helpers/scheduling_helpers.hpp"

#include "heuristics/LNS/relaxation_policy_factories.hpp"
#include "heuristics/LNS/PerfectRelaxationOracle.hpp"
#include "heuristics/warmstart.hpp"

using namespace tempo;


// implementation of a scheduling solver
int main(int argc, char *argv[]) {
  namespace lns = tempo::lns;
  auto parser = tempo::getBaseParser();
  bool profileHeuristic;
  lns::RelaxationPolicyParams policyParams{
    .decayConfig = lns::PolicyDecayConfig(), .numScheduleSlices = 4, .allTaskEdges = false
  };
  lns::RelaxPolicy policyType;
  bool useOracle = false;
  double oracleEpsilon = 0;
  double sporadicIncrement = 0;
  StatsConfig statsConfig{
    .displayStats = false, .regionTimeout = 0, .nRegionThreads = 0, .nRegionSaverFails = 0,
    .solPath = {}, .regionRecordPath = {}, .regionExecutionPolicy = lns::ExecutionPolicy::Lazy
  };
  cli::detail::configureParser(parser, cli::SwitchSpec("heuristic-profiling", "activate heuristic profiling",
                                                       profileHeuristic, false),
                               cli::ArgSpec("fix-decay", "relaxation ratio decay",
                                            false, policyParams.decayConfig.decay, 0.5),
                               cli::ArgSpec("fix-ratio", "initial relaxation ratio",
                                            false, policyParams.decayConfig.fixRatio, 0.5),
                               cli::ArgSpec("decay-min-fail", "lower bound solver failure rate for ratio decay config",
                                            false, policyParams.decayConfig.minFailRatio),
                               cli::ArgSpec("decay-max-fail", "upper bound solver failure rate for ratio decay config",
                                            false, policyParams.decayConfig.maxFailRatio),
                               cli::SwitchSpec("decay-on-success", "whether to decrease fix rate even on success",
                                               policyParams.decayConfig.decreaseOnSuccess, false),
                               cli::ArgSpec("retry-limit", "number of fails before decreasing relaxation ratio",
                                            false, policyParams.decayConfig.retryLimit),
                               cli::ArgSpec("decay-mode", "relaxation ratio decay mode on failure", false,
                                            policyParams.decayConfig.decayMode),
                               cli::ArgSpec("relax-slices", "number of schedule slices",
                                            false, policyParams.numScheduleSlices, 4),
                               cli::SwitchSpec("fix-all-task-edges",
                                               "whether to fix all task edges or only those between fixed tasks",
                                               policyParams.allTaskEdges, false),
                               cli::ArgSpec("lns-policy", "lns relaxation policy", true, policyType),
                               cli::ArgSpec("optimal-solution", "location of optimal solution (e.g. for oracle)", false,
                                            statsConfig.solPath),
                                            cli::SwitchSpec("oracle", "use perfect relaxation oracle", useOracle, false),
                               cli::ArgSpec("oracle-epsilon", "LNS oracle policy epsilon", false, oracleEpsilon),
                               cli::ArgSpec("sporadic-increment", "sporadic root search probability increment", false,
                                            sporadicIncrement),
                               cli::SwitchSpec("stats", "enable policy statistics", statsConfig.displayStats, false),
                               cli::ArgSpec("local-optimum-timeout",
                                            "timeout for search for local optimum (stats only)", false,
                                            statsConfig.regionTimeout),
                               cli::ArgSpec("local-optimum-threads",
                                            "number of threads for local optimum search (stats only)", false,
                                            statsConfig.nRegionThreads),
                               cli::ArgSpec("local-optimum-exec",
                                            "execution policy for local optimum search (stats only)", false,
                                            statsConfig.regionExecutionPolicy),
                               cli::ArgSpec("save-regions", "destination file where lns regions are saved to", false,
                                            statsConfig.regionRecordPath),
                               cli::ArgSpec("region-save-fails",
                                            "minimum number of fails before regions are recorded", false,
                                            statsConfig.nRegionSaverFails)
  );

  parser.parse(argc, argv);
  Options opt = parser.getOptions();
  if (opt.print_par) {
    std::cout << opt << std::endl;
  }
  Solver<> S(opt);
  seed(opt.seed);

  // an interval standing for the makespan of schedule
    auto makespan{S.newNumeric(0,Constant::Infinity<int>)};
    auto origin{S.newConstant(0)};
    auto schedule{S.between(origin, makespan)};
//  auto schedule{S.newInterval(0, Constant::Infinity<int>, 0, 0, 0,
//                              Constant::Infinity<int>)};

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

  //    for(auto i : intervals) {
  //        std::cout << i.id() << ": " << i << std::endl;
  //    }
  //    exit(1);

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

  if (opt.ub == -1) {
    // TODO this might produce invalid bounds for e.g. jstl
    int trivialUb = 0;
    for (const auto &t : intervals) {
        auto maxDur = t.maxDuration(S);
        if (maxDur == Constant::Infinity<Time>) {
            trivialUb = maxDur;
            break;
        }

        trivialUb += maxDur;
    }

    S.post(schedule.end.before(trivialUb));
  } else if (opt.ub != Constant::Infinity<int>) {
    S.post(schedule.end.before(opt.ub));
    //        S.set(schedule.end.before(schedule.start, -opt.ub));
    //        S.propagate();

    //        std::cout << (schedule.end.before(schedule.start, -opt.ub)) <<
    //        std::endl;
  }

  // HACK!!: add an edge from origin to end to allow propagation of
  // shortest path
  S.set({0, 1, S.numeric.upper(1)});

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
      heuristics::warmstartDisjunctive(S, schedule, intervals, ub);
      //            warmstart(S, schedule, by_resource, ub);
    } catch (Failure<int> &f) {
      //            std::cout << " optimal solution found in a greedy run\n";
      optimal = true;
    }
  }

    if(not optimal) {
        MinimizationObjective<int> objective(schedule.duration);
        if (not useOracle) {
            std::cout << "-- using relaxation policy " << policyType << std::endl;
            auto policy = lns::make_relaxation_policy(policyType, intervals, resources, policyParams, opt.verbosity);
            if (sporadicIncrement != 0) {
              std::cout << "-- root search probability increment " << sporadicIncrement << std::endl;
              runLNS(lns::make_sporadic_root_search(sporadicIncrement, std::move(policy)), S,
                     objective, statsConfig);
            } else {
              runLNS(policy, S, objective, statsConfig);
            }
        } else {
            std::cout << "-- using perfect relaxation oracle" << std::endl;
            const auto sol = serialization::deserializeFromFile<serialization::Solution<int>>(statsConfig.solPath);
            lns::PerfectRelaxationOracle policy(toSolution(sol, opt), schedule.duration,
                                                booleanVarsFromResources(resources),
                                                policyParams.decayConfig.fixRatio, oracleEpsilon);
            runLNS(policy, S, objective, statsConfig);
        }
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

}
