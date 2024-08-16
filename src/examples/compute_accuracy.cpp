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
#include "util/parsing/path.hpp"
#include "util/parsing/tsptw.hpp"
#include "helpers/cli.hpp"
#include "util/Profiler.hpp"
#include "helpers/shell.hpp"
#include "helpers/git_sha.hpp"

using namespace tempo;

void build_model(Solver<>& S, Interval<>& schedule) {
    
    auto& opt(S.getOptions());
    
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
      for (auto j : tasks) {
        scope.push_back(intervals[j]);
      }
      auto no_overlap{NoOverlap(schedule, scope, resource_transitions[i++])};
      resources.push_back(no_overlap);
      S.post(no_overlap);
      scope.clear();
    }
    
}


boolean_state test_branch(Options& opt, const int makespan, std::vector<bool>& signs, std::vector<var_t>& vars) {
    
    Solver<> S(opt);
    
    auto schedule{S.newInterval(0, makespan, 0, 0, 0, makespan)};

    build_model(S, schedule);
    
    
    auto x{vars.begin()};
    auto s{signs.begin()};
    while(x + 1 < vars.end()) {
        auto constraint{S.boolean.getLiteral(*s, *x)};
        S.set(constraint);
        ++x;
        ++s;
    }
    
    try {
        auto constraint{S.boolean.getLiteral(not *s, *x)};
        S.post(constraint);
    } catch(Failure<int> &f) {
        return FalseState;
    }
    
    return S.satisfiable();
}

int solve(Options& gopt) {
    Options opt{gopt};
    opt.restart_policy = "no";
    opt.primal_boost = false;
    opt.greedy_runs = 0;
//    opt.instance_file = ifilename;
    
    long unsigned int num_right_decisions{0};
    long unsigned int num_wrong_decisions{0};
    long unsigned int cp_in_wasted_restarts{0};
    long unsigned int previous_num_cp{0};
    
    std::vector<std::vector<Literal<int>>> right_branches;
    std::ofstream outfile(opt.dbg_file);
    
    
    Solver<> S(opt);
    
    auto schedule{S.newInterval(0, Constant::Infinity<int>, 0, 0, 0, Constant::Infinity<int>)};
    
    
    SubscriberHandle solutionHandler(S.SolutionFound.subscribe_handled([&](const auto &) {
        num_right_decisions += S.numDecision();
          unsigned long num_wrong{0};
          for (auto &branches : right_branches)
          num_wrong += branches.size();
          num_wrong_decisions += num_wrong;
          outfile << S.numeric.lower(schedule.duration.id()) << " " << S.num_choicepoints - cp_in_wasted_restarts << " " << S.numDecision() << " " << num_wrong;
          for(auto decision : S.getDecisions())
              outfile << " " << decision.sign() << " " << decision.variable();
          outfile << std::endl;
    }));
    
    SubscriberHandle failHandler(S.ClauseAdded.subscribe_handled([&](const auto &learnt_clause) {
        right_branches.resize(S.numDecision() + 1);
        right_branches.back().push_back(learnt_clause[0]);
    }));
    
    SubscriberHandle restartHandler(S.SearchRestarted.subscribe_handled([&](const auto on_solution) {
        if(on_solution) {
            previous_num_cp = S.num_choicepoints;
        } else {
            cp_in_wasted_restarts += (S.num_choicepoints - previous_num_cp);
        }
    }));
    
        
        

        build_model(S, schedule);
    
    
    S.minimize(schedule.duration);
    
    return S.numeric.lower(schedule.duration);
}

 
void crunch_numbers(Options& opt) {
    opt.verbosity = 0;
    
    
    std::ifstream infile(opt.dbg_file);
    
        
    int makespan;
    unsigned long num_cp;
    unsigned branch_length;
    unsigned num_wrong;
    bool sign;
    var_t var;
    
    unsigned total_right{0};
    unsigned total_wrong{0};
    unsigned total_irrelevant_sat{0};
    
    std::vector<var_t> vars;
    std::vector<bool> signs;
    
    while(true) {
        
        unsigned irrelevant{0};
        
        infile >> makespan;
        infile >> num_cp;
        
        if(not infile.good())
            break;
        
        infile >> branch_length;
        infile >> num_wrong;
        for(unsigned i{0}; i<branch_length; ++i) {
            infile >> sign;
            infile >> var;
            
            signs.push_back(sign);
            vars.push_back(var);
            
            
            
            auto sat{test_branch(opt, makespan, signs, vars)};
            
            if(sat == TrueState) {
                ++irrelevant;
            }
        }
        
        
        total_wrong += num_wrong;
        total_right += (branch_length - irrelevant);
        total_irrelevant_sat += irrelevant;
        
        
        std::cout << "obj = " << makespan << " relevant = " << total_right + total_wrong << "/" << num_cp
        << " " << total_wrong << " unsat trees of avg size "
        << (num_cp - total_right - total_irrelevant_sat)/total_wrong
        << " accuracy = " ;
        if(total_right + total_wrong == 0)
            std::cout << "n/a";
        else
            std::cout << static_cast<double>(total_right) / static_cast<double>(total_right + total_wrong) ;
        std::cout << std::endl;
        
        
        vars.clear();
        signs.clear();
    }
    
}



// implementation of a scheduling solver
int main(int argc, char *argv[]) {
    auto parser = tempo::getBaseParser();
    
    parser.parse(argc, argv);
    Options opt = parser.getOptions();
    
    seed(opt.seed);
    
    solve(opt);
    
    crunch_numbers(opt);
    
}
