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
#include "helpers/cli.hpp"
#include "heuristics/Greedy.hpp"
#include "heuristics/relaxation_policy_factories.hpp"
#include "util/parsing/psplib.hpp"
#include "util/parsing/rcpsp.hpp"

using namespace tempo;



template <typename T>
std::string prettyJob(const Interval<T> &task, const Solver<T> &S,
                      const bool dur_flag) {
  std::stringstream ss;

  auto est{S.numeric.lower(task.start)};
  auto lst{S.numeric.upper(task.start)};
  auto ect{S.numeric.lower(task.end)};
  auto lct{S.numeric.upper(task.end)};
    
  ss << "[";

  if (est == lst)
    ss << est;
  else
    ss << est << "-" << lst;
  ss << "..";
  if (ect == lct)
    ss << ect;
  else
    ss << ect << "-" << lct;
  ss << "]";

  if (dur_flag) {
    auto pmin{S.numeric.lower(task.duration)};
    auto pmax{S.numeric.upper(task.duration)};
      if(pmin == pmax)
          ss << "(" << pmin << ")";
      else
          ss << " (" << pmin << "-" << pmax << ")";
  }

  return ss.str();
}

template<typename T>
void printJobs(const Solver<T>& S, const std::vector<Interval<T>>& intervals) {
  int i{0};
for(auto task : intervals)
    std::cout << "job" << ++i << ": " << prettyJob(task, S, true) << std::endl;
}


template<typename T>
void checkResource(const Solver<T>& S, const std::vector<Interval<T>>& tasks, const std::vector<NumericVar<T>>& demands, const T capacity) {

    std::vector<std::pair<T,T>> events;
    int i{0};
    for (auto job : tasks) {
        
        std::cout << job << std::endl;
        
        auto s{S.numeric.lower(job.start)};
        auto d{S.numeric.lower(job.duration)};
        auto e{S.numeric.lower(job.end)};
        auto dem{S.numeric.lower(demands[i])};
        
        if(s+d!=e) {
            std::cout << "bug span task " << i << ": " << job << std::endl;
        }
        
        events.emplace_back(s,dem);
        events.emplace_back(e,-dem);
        
        ++i;
    }
    
    std::sort(events.begin(), events.end(),
              [&](const std::pair<T,T>& a, const std::pair<T,T>& b) {
        return a.first < b.first or (a.first == b.first and a.second < b.second);
              });
    
    T profile{0};
    for(auto e : events) {
        profile += e.second;
        std::cout << e.first << " " << profile << std::endl;
        if(profile > capacity) {
            std::cout << "bug at time " << e.first << ": profile=" << profile << std::endl;
            exit(1);
        }
        
    }
    
}

template<typename T>
void printResources(const Solver<T>& S, const std::vector<std::vector<Interval<T>>>& resource_tasks, const std::vector<std::vector<NumericVar<T>>>& resource_demands, std::vector<int>& resource_capacities) {
  int i{0};

    auto demands{resource_demands.begin()};
  for (auto &tasks : resource_tasks) {
    ++i;
    if (tasks.size() > 1) {

      std::vector<int> order;
        int j{0};
      for (auto job : tasks) {
          if (S.boolean.value(job.exist)) {
              order.push_back(j);
//              order.push_back(job);
          }
          ++j;
      }

      std::sort(order.begin(), order.end(),
                [&](const int a, const int b) {
                  return S.numeric.lower(tasks[a].start) <
                         S.numeric.lower(tasks[b].start);
                });

      std::cout << "resource " << i << " (" << resource_capacities[i] << "):";
      for (auto o : order) {
        std::cout << " " << S.numeric.lower((*demands)[o]) << "x"
                  << prettyJob(tasks[o], S, false);
      }
      std::cout << std::endl;
    }
    
      checkResource(S, tasks, *demands, resource_capacities[i]);
      
      ++demands;
  }
}


// implementation of a scheduling solver
int main(int argc, char *argv[]) {
    
    namespace h = tempo::heuristics;
    auto parser = tempo::getBaseParser();
    bool useLNS;
//    bool tt_reasoning;
    h::RelaxationPolicyParams policyParams;
    h::RelaxPolicy policyType;
    cli::detail::configureParser(parser,
//                                 cli::SwitchSpec("no-ttef", "switch tt reasoning off in edge-finding", false, tt_reasoning, true),
                                 cli::SwitchSpec("lns", "activate large neighborhood search",
                                                         useLNS, false),
                                 cli::ArgSpec("relax-decay", "relaxation ratio decay",
                                              false, policyParams.ratioDecay, 0.5),
                                 cli::ArgSpec("relax-ratio", "initial relaxation ratio",
                                              false, policyParams.relaxRatio, 0.5),
                                 cli::ArgSpec("relax-slices", "number of schedule slices",
                                              false, policyParams.numScheduleSlices, 4),
                                 cli::ArgSpec("lns-policy", "lns relaxation policy", false, policyType, h::RelaxPolicy::RandomTasks));
    

    parser.parse(argc, argv);
    Options opt = parser.getOptions();
    Solver<> S(opt);
    seed(opt.seed);

  // an interval standing for the makespan of schedule
  auto schedule{S.newInterval(0, Constant::Infinity<int>, 0, 0, 0,
                              Constant::Infinity<int>)};

  // depending on the option "input-format", parse a disjunctive scheduling
  // instance, and collect resources and interval objects
  std::vector<CumulativeExpression<>> resources;
  std::vector<std::vector<size_t>> tasks_requirements;
    std::vector<std::vector<int>> task_demands;
    std::vector<int> resource_capacities;
  std::vector<Interval<>> intervals;
  std::vector<std::pair<int, int>> precedences;
    std::vector<std::vector<int>> graph;

  if (opt.input_format == "rcp")
    rcpsp::parse(opt.instance_file, S, schedule, intervals, tasks_requirements,
                 task_demands, resource_capacities, precedences, graph);
  else
    psplib::parse(opt.instance_file, S, schedule, intervals, tasks_requirements,
                  task_demands, resource_capacities, precedences, graph);
    
    for(auto &neighbors : graph) {
        std::sort(neighbors.begin(), neighbors.end());
    }

  //      for(auto i : intervals) {
  //          std::cout << i.id() << ": " << i << std::endl;
  //      }

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
 
  if (opt.print_mod) {
    std::cout << S << std::endl;
  }
    
    
    int ub_makespan{0};
      for (auto &j : intervals) {
        if (j.maxDuration(S) == Constant::Infinity<int>) {
            ub_makespan = Constant::Infinity<int>;
          break;
        }
          ub_makespan += j.maxDuration(S);
      }
    
//
  auto optimal{false};
    if (opt.greedy_runs > 0) {
        
        S.initializeSearch();
        
        ScheduleGenerationScheme<int> sgs(S, tasks_requirements, task_demands, resource_capacities, intervals, precedences, graph);
        
        for(auto i{0}; i<opt.greedy_runs; ++i) {
            auto makespan{sgs.run()};
            sgs.clear();
            
            if(makespan < ub_makespan) {
                
                std::cout << "-- load improving sgs solution " << makespan << std::endl;
                
                sgs.load();
                ub_makespan = makespan;
                S.num_choicepoints += sgs.num_insertions;
                if (opt.verbosity >= Options::NORMAL) {
                    std::cout << std::setw(10) << ub_makespan;
                    S.displayProgress(std::cout);
                }
            }
        }
        
//        ub_makespan = sgs.best_makespan;
        
        try {
            
//            std::cout << "post ub = " << ub_makespan-1 << "\n";
            
            S.post(schedule.duration < ub_makespan);
        } catch(Failure<int>& f) {
            if (opt.verbosity >= Options::QUIET)
                S.displaySummary(std::cout, "optimal");
            optimal = true;
        }
        
        S.num_choicepoints = 0;
    }

  // search
    if (not optimal) {
        
        auto ub{std::min(opt.ub, ub_makespan)};
        S.post(schedule.end.before(ub));
        
        if(useLNS) {
            
            MinimizationObjective<int> objective(schedule.duration);
            auto policy = h::make_relaxation_policy(policyType, intervals, resources, policyParams);
            std::cout << "-- using relaxation policy " << policyType << std::endl;
            S.largeNeighborhoodSearch(objective, policy);
            
        } else {
                
            S.minimize(schedule.duration);
            
        }
    }
  
    if (opt.print_sol) {
      printJobs(S, intervals);
      printResources(S, resource_tasks, resource_demands, resource_capacities);
    }
}
