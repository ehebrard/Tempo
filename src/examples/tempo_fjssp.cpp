/************************************************
 * Tempo scheduling solver
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
#include "constraints/Cardinality.hpp"
#include "util/parsing/fjssp.hpp"

using namespace tempo;


template <typename T>
std::string prettyJob(const Interval<T> &task, const Solver<T> &S,
                      const bool dur_flag) {
  std::stringstream ss;
    
    if(S.boolean.value(task.exist)) {
        
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
            ss << " (" << pmin << "-" << pmax << ")";
        }
    } else {
        ss << "not present (" << task.exist << "=false)";
    }

  return ss.str();
}

template<typename T>
void printJobs(const Solver<T>& S, const std::vector<Interval<T>>& intervals) {
  int i{0};
  for (auto task : intervals) {
    std::cout << "job" << ++i << " (" << task
              << "): " << prettyJob(task, S, true) << std::endl;
  }
}


// implementation of a scheduling solver
int main(int argc, char *argv[]) {

  Options opt = tempo::parse(argc, argv);
  Solver<> S(opt);

  // an interval standing for the makespan of schedule
  auto schedule{S.newInterval(0, Constant::Infinity<int>, 0, 0, 0,
                              Constant::Infinity<int>)};

    
    
    std::vector<std::vector<int>> resource_alloc;
    std::vector<std::vector<int>> duration_mode;
    std::vector<int> job_size;
    

    int num_machine;
    
  fjssp::parse(opt.instance_file, num_machine, job_size, resource_alloc, duration_mode);
    
    std::vector<NoOverlapExpression<int>> resources; //(num_machine);
    for(auto i{0}; i<num_machine; ++i) {
        resources.push_back(NoOverlap(schedule));
    }
    
    std::vector<std::vector<BooleanVar<int>>> allocation_vars;
    
    std::vector<Interval<int>> intervals;
    std::vector<Interval<int>> previous_jobs;
    
    
    
    int i{0};
    for(auto nj : job_size) {
//        previous_jobs.push_back(schedule.start);
        for(auto j{0}; j<nj; ++j) {
            auto s{S.newNumeric()};
            
            if(j==0) {
                S.post(s >= schedule.start);
            } else if(previous_jobs.size() == 1) {
                S.post(s >= previous_jobs.begin()->end);
            } else {
                int mindur{Constant::Infinity<int>};
                for(auto t : previous_jobs) {
                    mindur = std::min(mindur, t.duration.min(S));
                }
                S.post(s.after(previous_jobs.begin()->start, mindur));
                for(auto t : previous_jobs) {
                    if(t.duration.min(S) > mindur) {
                        S.post(t.exist.implies(s >= t.end));
                    }
                }
            }
            
//            std::cout << S << std::endl;
            
            previous_jobs.clear();
            std::vector<BooleanVar<int>> allocation;
            
            if(resource_alloc[i].size() == 1) {
                // non-flexible task
                auto dur{duration_mode[i][0]};
                auto mach{resource_alloc[i][0]};
                
                auto t{S.between(s, s+dur)};
                intervals.push_back(t);
                resources[mach].push_back(t);
                previous_jobs.push_back(t);
                
            } else {
//                allocation.clear();
                for(unsigned k{0}; k<resource_alloc[i].size(); ++k) {
                    allocation.push_back(S.newBoolean());
                    S.addToSearch(allocation.back());
                    
                    auto dur{duration_mode[i][k]};
                    auto mach{resource_alloc[i][k]};
                    
                    auto t{S.between(s, s+dur, allocation.back())};
                    
                    
//                    std::cout << " ---> " << t.id() << ": " << allocation.back() << " / " << t.exist << std::endl;
                    
                    intervals.push_back(t);
                    resources[mach].push_back(t);
                    previous_jobs.push_back(t);
                }
                
                
//                S.post(BigOr(allocation));
                S.post(AtMost(static_cast<int>(allocation.size()-1),allocation));
                S.post(AtLeast(1,allocation));
                
                
            }
            
            allocation_vars.push_back(allocation);
            
            if(j == nj-1) {
                if(previous_jobs.size() == 1) {
                    S.post(schedule.end >= previous_jobs.begin()->end);
                } else {
                    int mindur{Constant::Infinity<int>};
                    for(auto t : previous_jobs) {
                        mindur = std::min(mindur, t.duration.min(S));
                    }
                    S.post(schedule.end.after(previous_jobs.begin()->start, mindur));
                    for(auto t : previous_jobs) {
                        if(t.duration.min(S) > mindur) {
                            S.post(t.exist.implies(schedule.end >= t.end));
                        }
                    }
                }
//
//                for(auto t : previous_jobs) {
//                    S.post(t.end <= schedule.end);
//                }
            }
            
            ++i;
        }
    }
    
    for(auto& resource_constraint : resources) {
        
//        std::cout << resource_constraint.size() << std::endl;
        
        S.post(resource_constraint);
    }
    
   

  if (opt.print_mod) {
    std::cout << S << std::endl;
  }

  // search
  S.minimize(schedule.duration);

  if (opt.print_sol) {
    printJobs(S, intervals);
//    printResources(S, intervals, resource_alloc);
  }
}
