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

#include "Solver.hpp"
#include "heuristics/Greedy.hpp"
#include "util/parsing/jsp.hpp"
#include "util/parsing/jstl.hpp"
#include "util/parsing/osp.hpp"
#include "util/parsing/path.hpp"
#include "util/parsing/tsptw.hpp"


using namespace tempo;


enum outcome { TRIVIALLY_UNSAT = 0, PROPAG_UNSAT, SEARCH_UNSAT, SAT };


outcome check_clause(std::vector<DistanceConstraint<int>>& clause, Solver<int>& S) {
    
//    std::cout << std::endl;
    auto init = S.saveState();
    auto result{outcome::SAT};
    for(auto c : clause) {
        try {
            
//            std::cout << "post " << c << std::endl;
            S.post(c);
        } catch (Failure<int> &f) {
            result = outcome::TRIVIALLY_UNSAT;
        }
    }
    if (result != outcome::TRIVIALLY_UNSAT) {
        try {
            S.propagate();
        } catch (Failure<int> &f) {
            result = outcome::PROPAG_UNSAT;
        }
        if (result != outcome::PROPAG_UNSAT) {
            
////            std::cout << S << std::endl;
//            S.displayDomains(std::cout);
//            S.displayBranches(std::cout);
////            std::cout << std::endl;
            
            if (not S.satisfiable()) {
                result = outcome::SEARCH_UNSAT;
            }
        }
    }
    S.restoreState(init);
    S.clauses.clear();
    
    return result;
}



int main(int argc, char *argv[]) {
    Options opt = tempo::parse(argc, argv);
    
    size_t num_fails{0};
    size_t num_trivial{0};
    size_t num_search{0};
    
    seed(opt.seed);
    
    std::ifstream cl_file(opt.dbg_file, std::ifstream::in);
    
    opt.dbg_file = "";
//    opt.learning = false;
    opt.verbosity = 3;
//    opt.restart_policy = "no";
    
    
    
    
    Solver<int> S(opt);
    
    // an interval standing for the makespan of schedule
    auto makespan{S.newNumeric(0,Constant::Infinity<int>)};
    auto origin{S.newConstant(0)};
    auto schedule{S.between(origin, makespan)};
    //    auto schedule{S.newInterval(0, Constant::Infinity<int>, 0, 0, 0,
    //                                Constant::Infinity<int>)};
    
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
    }
    
    resource_transitions.resize(resource_tasks.size());
    
    
    index_t i{0};
    std::vector<Interval<int>> scope;
    for (auto &tasks : resource_tasks) {
        for(auto j : tasks) {
            scope.push_back(intervals[j]);
        }
        auto no_overlap{NoOverlap(schedule, scope, resource_transitions[i++])};
        resources.push_back(no_overlap);
        S.post(no_overlap);
        scope.clear();
    }
    
    var_t x, y;
    int t, n, d;
    
    std::vector<DistanceConstraint<int>> previous_clause;
    std::vector<DistanceConstraint<int>> current_clause;
    
    std::map<std::pair<int,int>,int> lit_map;
    
    int line{0};
    do {
        current_clause.clear();
        cl_file >> t;
        cl_file >> n;
        if (not cl_file.good()) {
            std::cout << "end of file\n";
            break;
        }
        
        for (auto i{0}; i < n; ++i) {
            cl_file >> x;
            cl_file >> y;
            cl_file >> d;
            current_clause.emplace_back(x,y,d);
        }
        
        auto nf{S.num_fails};
        
//        if(line == 47) {
//            S.debug_flag = 1;
//        }
        
        auto res = check_clause(current_clause, S);
        
        if(res == TRIVIALLY_UNSAT) {
            ++num_trivial;
        } else if(res == SEARCH_UNSAT) {
            num_fails += (S.num_fails - nf);
            ++num_search;
        } else if(res == SAT) {
            std::cout << "cl " << line << " ("
            << (t == 1 ? "minimized" : (t == 2 ? "reason" : "uip"))
            << "): ";
            std::cout << "bug!\n";
            for(auto c : current_clause) {
                std::cout << "> " << c << std::endl;
            }
            
            if (opt.print_sol) {
                for (auto i : intervals) {
                    std::cout << "t" << i.start.id() << ": "
                    << S.numeric.solutionLower(i.start) << ".."
                    << S.numeric.solutionLower(i.end) << " ("
                    << S.numeric.solutionLower(i.duration) << ")\n";
                }
            }
            
            exit(1);
            
            if(t == 1) {
                std::cout << "minimized clause, check literals:\n current:";
                for(auto c : current_clause) {
                    std::cout << " " << c;
                }
                std::cout << "\nprevious:";
                for(auto c : previous_clause) {
                    std::cout << " " << c;
                }
                std::cout << "\n";
                
                
                lit_map.clear();
                for(auto c : current_clause) {
                    lit_map[{c.from,c.to}] = c.distance;
                }
                std::vector<DistanceConstraint<int>> clause;
                for(auto c : previous_clause) {
                    auto minimized{false};
                    if(not lit_map.contains({c.from,c.to})) {
                        std::cout << c << " has been removed\n";
                        minimized = true;
                    } else if(lit_map[{c.from,c.to}] > c.distance) {
                        auto d{DistanceConstraint(c.from, c.to, lit_map[{c.from,c.to}])};
                        std::cout << c << " has been changed to " << d << "\n";
                        clause.push_back(d);
                        minimized = true;
                    }
                    if(minimized) {
                        for(auto d : previous_clause) {
                            if(d != c)
                                clause.push_back(d);
                        }
                        std::cout << "check:";
                        for(auto d : clause) {
                            std::cout << " " << d;
                        }
                        std::cout << "\n";
                        
                        auto res = check_clause(clause, S);
                        
                        if(res == SAT) {
                            std::cout << "bug here!\n";
                            break;
                        }
                        
                        std::cout << "\n";
                        
                    }
                    clause.clear();
                }
                
            }
            
            std::cout << "error?\n";
            exit(1);
        }
    
        
        std::swap(current_clause, previous_clause);
    
        ++line;
        
        std::cout << line << ": " << num_trivial << " trivial, "
        << (line - num_trivial - num_search) << " easy, " << num_search
        << " hard (" << num_fails / num_search << ")" << std::endl;
    } while(true);
    
}
