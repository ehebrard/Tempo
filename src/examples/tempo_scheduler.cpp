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
#include "heuristics/Greedy.hpp"
#include "util/parsing/osp.hpp"
#include "util/parsing/jsp.hpp"
#include "util/parsing/jstl.hpp"
//#include "util/parsing/dimacs.hpp"

using namespace tempo;





int main(int argc, char *argv[]) {


    
    
  Options opt = tempo::parse(argc, argv);

    Solver<> S(opt);

    auto schedule{S.newJob()};
    std::vector<DisjunctiveResource<>> resources;
    std::vector<Job<>> jobs;

    if (opt.input_format == "osp") {
        osp::parse(opt.instance_file, S, schedule, jobs, resources);
    }
    else if (opt.input_format == "jsp") {
        jsp::parse(opt.instance_file, S, schedule, jobs, resources);
    }
//    else if (opt.input_format == "tsptw") {
//        tsptw::parse(opt.instance_file, S, schedule, jobs, resources);
//    } 
    else if (opt.input_format == "jstl") {
        jstl::parse(opt.instance_file, S, schedule, jobs, resources);
    }

    std::vector<BooleanVar<>> X;
    for (auto &R : resources) {
      //        auto disjuncts{X.end()};
      auto s{X.size()};
      R.createOrderVariables(S, X);
      if (opt.edge_finding) {
        S.postEdgeFinding(schedule, R.begin(), R.end(), X.begin() + s);
      }
      if (opt.transitivity) {
        S.postTransitivity(schedule, R.begin(), R.end(), X.begin() + s);
      }
    }

    for (auto x : X)
      S.addToSearch(x);

    MakespanObjective<int> duration(schedule, S);
    
    
//    std::cout << S << std::endl;
//    S.displayPrecedences(std::cout);
//
//    exit(1);
        
    Greedy greedy_insertion(S);
    
    greedy_insertion.addJobs(jobs);
    
    size_t k{0};
    for (auto &R : resources) {
        auto n = k + R.size() * (R.size()-1) / 2;
        greedy_insertion.addResource(X.begin()+k, X.begin()+n);
        k = n;
    }
    
    auto trivial_ub{0};
    for(auto& j : jobs)
        trivial_ub += j.minDuration();
    
    // set the user defined ub
    auto ub{std::min(opt.ub, trivial_ub)};
    
//    std::cout << ub << std::endl;
    
    
    S.set(schedule.end.before(ub));
    
    S.initializeSearch();
    
    S.propagate();
    
//    std::cout << S << std::endl;
    
//    exit(1);
    
        
    // find a ub
    for(auto i{0}; i<opt.greedy_runs; ++i) {
        auto st{S.saveState()};
//        auto sat{((tempo::random()%2) ? greedy_insertion.runEarliestStart() : greedy_insertion.runLatestEnd())};
//        auto sat{greedy_insertion.runLatestEnd()};
        auto sat{greedy_insertion.runEarliestStart()};
        if(sat) {
            if(schedule.getEarliestEnd(S) <= ub) {
                ub = schedule.getEarliestEnd(S) - 1;
                std::cout << std::setw(10) << (ub+1);
                S.displayProgress(std::cout);
            }
        }
        S.restoreState(st);
    }
    
    // set the ub (again)
    S.set(schedule.end.before(ub));

    S.optimize(duration);
}
