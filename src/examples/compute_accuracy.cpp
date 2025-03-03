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
#include "util/IntFinity.hpp"

//#define VERBOSE true

using namespace tempo;

void build_model(Solver<> &S, Interval<> &schedule) {

  auto &opt(S.getOptions());

  // depending on the option "input-format", parse a disjunctive scheduling
  // instance, and collect resources and interval objects
  std::vector<NoOverlapExpression<>> resources;
  std::vector<std::vector<size_t>> resource_tasks;
  std::vector<Interval<>> intervals;
  std::vector<int> weights;
  std::vector<std::vector<std::vector<int>>> resource_transitions;
    
//    std::cout << "parsing\n";

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
    
//    std::cout << "parse ok\n";

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
    
//    std::cout << "resources ok\n";
}

void read_branch(std::istream &in, 
                 unsigned& branch_length,
                 unsigned& num_wrong,
                 std::vector<bool> &dsigns,
                 std::vector<var_t> &dvars,
                 std::vector<std::vector<bool>> &rsigns,
                 std::vector<std::vector<var_t>> &rvars) {

  dsigns.clear();
  dvars.clear();
  rsigns.clear();
  rvars.clear();
    bool end_on_righ_flag;

    in >> end_on_righ_flag;
  in >> branch_length;
  unsigned wrongcount;
  in >> wrongcount;
    
    num_wrong = wrongcount;

//    std::cout << "end_on_righ_flag = " << end_on_righ_flag << std::endl;
//  std::cout << "branch_length = " << branch_length << std::endl;
//    std::cout << "wrongcount = " << wrongcount << std::endl;

  int s;
  int v;

  for (unsigned i{0}; i < branch_length + end_on_righ_flag; ++i) {
    rsigns.resize(rsigns.size() + 1);
    rvars.resize(rvars.size() + 1);
      
//      std::cout << rsigns.size() << " " << rsigns.back().size() << std::endl;

//      std::cout << "level " << i << std::endl;
      
    int nwrong;
    in >> nwrong;
      
//      std::cout << "w:" << nwrong << " ";
      
      wrongcount -= nwrong;
    for (auto j{0}; j < nwrong; ++j) {
      in >> s;
      in >> v;
      rsigns.back().push_back(s);
      rvars.back().push_back(v);
        
//        std::cout << s << ":" << v << " ";
    }

      if(i<branch_length) {
          in >> s;
          in >> v;
          
//          std::cout << s << ":" << v << " ";
          dsigns.push_back(s);
          dvars.push_back(v);
      }
//      std::cout << std::endl;
  }

  unsigned total{0};
  for (auto &r : rvars) {
    total += r.size();
  }

  if (wrongcount != 0 or total != num_wrong) {
    std::cout << "bug read branch (" << total << "/" << num_wrong << ")\n";
    exit(1);
  }
}

//#define VERBOSE true

bool satisfiable(Solver<>& S, Literal<int> constraint) {
    try {
        
#ifdef VERBOSE
        std::cout << " + " << constraint << std::endl;
#endif
        
      S.post(constraint);
    } catch (Failure<int> &f) {
        
#ifdef VERBOSE
        std::cout << "direct failure\n" ;
#endif
        
      return false;
    }
    
    return (S.satisfiable() == TrueState);
}

int test_branches(Options &opt, const int makespan, std::vector<bool> &dsigns,
                  std::vector<var_t> &dvars,
                  std::vector<std::vector<bool>> &rsigns,
                  std::vector<std::vector<var_t>> &rvars,
                  std::vector<int> &dec_level, std::vector<int> &num_mistakes,
                  const bool side) {

  int num_irrelevant{0};

  if (rsigns.size() < dsigns.size()) {
    std::cout << "BUG!!\n";
    exit(1);
    }

#ifdef VERBOSE
  std::cout << "\ntest branch\n";
#endif

  for (size_t i{0}; i < rsigns.size(); ++i) {
    Solver<> S(opt);

    auto end_sched{S.newNumeric(0, makespan)};
    auto schedule{S.between(S.zero(), end_sched)};

    build_model(S, schedule);

    if (num_mistakes.size() < S.boolean.size())
      num_mistakes.resize(S.boolean.size(), 0);

#ifdef VERBOSE
    std::cout << "up to level " << (i + 1) << std::endl;
#endif
      
    for (size_t j{0}; j < i; ++j) {
      // add the right branches
        
      for (size_t k{0}; k < rsigns[j].size(); ++k) {
        auto constraint{S.boolean.getLiteral(rsigns[j][k], rvars[j][k])};

#ifdef VERBOSE
        std::cout << " - deduction " << constraint << std::endl;
#endif
          
        S.set(constraint);
      }

      // add the left branch
      auto constraint{S.boolean.getLiteral(dsigns[j], dvars[j])};

#ifdef VERBOSE
      std::cout << " - branch " << constraint << std::endl;
#endif
        
      S.set(constraint);
    }

      auto n{rsigns[i].size()};
      if(i == dsigns.size())
          --n;
      
      // add the right branches
      for (size_t k{0}; k < n; ++k) {
        auto constraint{S.boolean.getLiteral(rsigns[i][k], rvars[i][k])};

#ifdef VERBOSE
        std::cout << " - deduction " << constraint << std::endl;
#endif
          
        S.set(constraint);
      }
      
      Literal<int> constraint;
      
      if(i == dsigns.size())
          constraint = S.boolean.getLiteral(rsigns[i].back(), rvars[i].back());
      else
          constraint = S.boolean.getLiteral(dsigns[i], dvars[i]);

      if(side == true) {
          if(not satisfiable(S, constraint)) {
              std::cout << "bug branch!!\n";
              exit(1);
          } 
#ifdef VERBOSE
          else {
              std::cout << "OK\n";
          }
#endif
      } else {
        if (satisfiable(S, ~constraint)) {
          ++num_irrelevant;
        } else {
          ++dec_level[i];
        }
      }
      
  }

  return num_irrelevant;
}

boolean_state test_branch(Options &opt, const int makespan,
                          std::vector<bool> &signs, std::vector<var_t> &vars) {

  Solver<> S(opt);

//    auto makespan{S.newNumeric(0,Constant::Infinity<int>)};
    auto origin{S.newConstant(0)};
    auto schedule{S.between(origin, makespan)};
//  auto schedule{S.newInterval(0, makespan, 0, 0, 0, makespan)};

  build_model(S, schedule);

  auto x{vars.begin()};
  auto s{signs.begin()};
  while (x + 1 < vars.end()) {
    auto constraint{S.boolean.getLiteral(*s, *x)};
    try {
      S.set(constraint);
    } catch (Failure<int> &f) {
      std::cout << "BUG: Failure when branchong on " << constraint << " b/c "
                << f.reason << std::endl;
      exit(1);
    }
    ++x;
    ++s;
  }

  try {
    auto constraint{S.boolean.getLiteral(not *s, *x)};
    S.post(constraint);
  } catch (Failure<int> &f) {
    return FalseState;
  }

  auto r{S.satisfiable()};

  return r;
}

int solve(Options& gopt, std::string& record_file) {
  Options opt{gopt};
  opt.restart_policy = "no";
  opt.primal_boost = false;
  opt.greedy_runs = 0;
  opt.learning = 0;
  //    opt.instance_file = ifilename;

  long unsigned int num_correct_decisions{0};
  long unsigned int num_wrong_decisions{0};
  long unsigned int cp_in_wasted_restarts{0};
  long unsigned int previous_total_cp{0};

  std::vector<std::vector<Literal<int>>> right_branches;

  std::stringstream buffer;

  std::ofstream outfile(record_file);

  Solver<> S(opt);

    auto makespan{S.newNumeric(0,Constant::Infinity<int>)};
    auto origin{S.newConstant(0)};
    auto schedule{S.between(origin, makespan)};
//  auto schedule{S.newInterval(0, Constant::Infinity<int>, 0, 0, 0,
//                              Constant::Infinity<int>)};

  SubscriberHandle solutionHandler(
      S.SolutionFound.subscribe_handled([&](const auto &) {
        num_correct_decisions += S.numDecision();
        unsigned long num_wrong{0};
        for (auto &branches : right_branches)
          num_wrong += branches.size();
        num_wrong_decisions += num_wrong;
        buffer << S.numeric.solutionLower(schedule.duration) << " "
               << S.num_choicepoints - cp_in_wasted_restarts << " "
               << (right_branches.size() > S.numDecision()) << " "
          << S.numDecision() << " "
          << num_wrong;
        for (unsigned i{0}; i < S.numDecision(); ++i) {
          if (right_branches.size() > i) {
            buffer << " " << right_branches[i].size();
            for (auto l : right_branches[i]) {
              buffer << " " << l.sign() << " " << l.variable();
            }
          } else {
            buffer << " 0";
          }
          buffer << " " << S.getDecisions()[i].sign() << " "
                 << S.getDecisions()[i].variable();
        }
          if(right_branches.size() > S.numDecision()) {
              buffer << " " << right_branches.back().size();
              for (auto l : right_branches.back()) {
                buffer << " " << l.sign() << " " << l.variable();
              }
          }
          
        buffer << std::endl;
        right_branches.clear();
      }));

  SubscriberHandle failCLHandler(
      S.ClauseAdded.subscribe_handled([&](const auto &learnt_clause) {
        right_branches.resize(S.numDecision() + 1);
        right_branches.back().push_back(learnt_clause[0]);
      }));

  SubscriberHandle failNOCLHandler(
      S.DeductionMade.subscribe_handled([&](const auto &lit) {
        right_branches.resize(S.numDecision());
        right_branches.back().push_back(lit);
      }));

  SubscriberHandle restartHandler(
      S.SearchRestarted.subscribe_handled([&](const auto on_solution) {
        if (on_solution) {
          previous_total_cp = S.num_choicepoints;
        } else {
          cp_in_wasted_restarts += (S.num_choicepoints - previous_total_cp);
        }
      }));

  build_model(S, schedule);

    if (opt.lb > -Constant::Infinity<int>) {
        S.post(schedule.duration >= opt.lb);
    }
    

  S.minimize(schedule.duration);

  auto obj{S.numeric.solutionLower(schedule.duration)};

  outfile << obj << std::endl << buffer.str();

  return obj;
}

void crunch_numbers(Options& opt, std::string& analyse_file) {
  auto v{opt.verbosity};
  opt.verbosity = 0;

  std::ifstream infile(analyse_file);

  int obj;
  infile >> obj;

    
    int bias{2};
//      std::cout << "obj = " << obj << std::endl;

  unsigned long prev_cp;
  intfinity<int> prev_makespan;
  int makespan{Constant::Infinity<int>};

  unsigned branch_length;

  unsigned num_wrong;
  unsigned num_correct;
  unsigned long num_cp;
  unsigned irrelevant_correct;

  unsigned total_correct{0};
  unsigned total_wrong{0};
  unsigned long total_cp{0};
  //    unsigned total_irrelevant_correct{0};

  unsigned total_branches{0};
    
    double biased_avg_correct{0};
    double biased_avg_wrong{0};
    

  std::vector<int> dec_level;
  std::vector<int> num_mistakes;

  //  bool sign;
  //  var_t var;

  std::vector<var_t> vars;
  std::vector<bool> signs;

  std::vector<bool> dsigns;
  std::vector<var_t> dvars;
  std::vector<std::vector<bool>> rsigns;
  std::vector<std::vector<var_t>> rvars;

  std::vector<double> var_ratio;
  int num_steps{5};
  var_ratio.resize(num_steps);
  std::iota(var_ratio.begin(), var_ratio.end(), 0);

  std::cout << "  obj."
            << "      gap |"
            << " branches"
            << "  rel."
            << "   #U"
            << "  size"
            << "   acc. |"
            << " branches(c)"
            << "  rel.(c)"
            << "   #U(c)"
            << "  size(c)"
            << "   acc.(c)\n";

  //    int branch_i{0};
  while (true) {
    if (KillHandler::instance().signalReceived()) {
        std::cout << "-- killed" << std::endl;
        break;
    }

    prev_cp = total_cp;
    prev_makespan = makespan;

    infile >> makespan;
    infile >> total_cp;

    //      std::cout << "makespan " << makespan << " total_cp " << total_cp <<
    //      std::endl;

    num_cp = total_cp - prev_cp;

    if (not infile.good())
      break;

    read_branch(infile, branch_length, num_wrong, dsigns, dvars, rsigns, rvars);

    if (dec_level.size() < rsigns.size())
      dec_level.resize(rsigns.size(), 0);

    irrelevant_correct = test_branches(
        opt, (prev_makespan - 1).get(), dsigns, dvars, rsigns, rvars, dec_level, num_mistakes, false);

    //      unsigned ttotal{0};
    for (size_t l{0}; l < rvars.size(); ++l) {
      if (not rvars[l].empty()) {
        for (auto x : rvars[l]) {

          //                  std::cout << "hello " << x << "/" <<
          //                  num_mistakes.size() << "\n";

          ++num_mistakes[x];
          ++dec_level[l];
        }
      }
      //          ttotal += rvars[l].size();
    }
      
      
      auto prev_total{static_cast<double>(total_correct + total_wrong)};
      

    num_correct = (branch_length - irrelevant_correct);
    total_correct += num_correct;
    total_wrong += num_wrong;
      
      
      if(biased_avg_correct == 0)
          biased_avg_correct = static_cast<double>(num_correct);
      else {
          biased_avg_correct = prev_total * biased_avg_correct + static_cast<double>((num_wrong + num_correct) * bias * num_correct);
          biased_avg_correct /= static_cast<double>(total_wrong + total_correct + (num_wrong + num_correct) * (bias-1));
      }
      
      if(biased_avg_wrong == 0)
          biased_avg_wrong = static_cast<double>(num_wrong);
      else {
          biased_avg_wrong = prev_total * biased_avg_wrong + static_cast<double>((num_wrong + num_correct) * bias * num_wrong);
          biased_avg_wrong /= static_cast<double>(total_wrong + total_correct + (num_wrong + num_correct) * (bias-1));
      }
      


    total_branches += branch_length;

    for (auto t{0}; t < num_steps; ++t) {
      unsigned nvars{0};
      for (auto n : num_mistakes) {
        if (n > t)
          ++nvars;
      }

      var_ratio[t] =
          static_cast<double>(nvars) / static_cast<double>(num_mistakes.size());
    }

    intfinity<std::size_t> avg_dec_level{0};
    for (size_t i{0}; i < dec_level.size(); ++i) {
      avg_dec_level += i * dec_level[i];
    }
    avg_dec_level /= (total_correct + total_wrong);

    //      if(ttotal != num_wrong) {
    //          std::cout << "bug (1)!\n" << ttotal << " " << num_wrong <<
    //          std::endl; exit(1);
    //      }
    //
    //      unsigned total{0};
    //      for(auto n : num_mistakes) {
    //          total += n;
    //      }
    //      if(total != total_wrong) {
    //          std::cout << "bug (2)!\n" << total << " " << total_wrong <<
    //          std::endl; exit(1);
    //      }

    std::cout << std::setw(6) << makespan << std::setw(9)
              << std::setprecision(3)
              << (static_cast<double>(makespan - obj) /
                  static_cast<double>(obj))
              << "  " << std::setw(9) << num_cp << std::setw(6)
              << (num_correct + num_wrong) << std::setw(5) << num_wrong;
    if (num_wrong == 0)
      std::cout << "   n/a";
    else
      std::cout << std::setw(6) << ((num_cp - branch_length) / num_wrong);

    if ((num_correct + num_wrong) == 0)
      std::cout << "    n/a";
    else
      std::cout << std::setw(7) << std::setprecision(4)
                << static_cast<double>(num_correct) /
                       static_cast<double>(num_correct + num_wrong);

    std::cout << "  " << std::setw(12) << total_cp << std::setw(9)
              << total_correct + total_wrong << std::setw(8) << total_wrong;

    if (total_wrong == 0)
      std::cout << "      n/a";
    else
      std::cout << std::setw(9) << ((total_cp - total_branches) / total_wrong);
    if ((total_correct + total_wrong) == 0)
      std::cout << "       n/a";
    else
      std::cout << std::setw(10) << std::setprecision(7)
                << static_cast<double>(total_correct) /
                       static_cast<double>(total_correct + total_wrong);
      
      

      if ((biased_avg_correct + biased_avg_wrong) == 0)
        std::cout << "       n/a";
      else
        std::cout << std::setw(10) << std::setprecision(7)
                  << static_cast<double>(biased_avg_correct) /
                         static_cast<double>(biased_avg_correct + biased_avg_wrong);

      std::cout << " " << std::setw(5) << avg_dec_level;
      for(int t{0}; t<num_steps; ++t)
          std::cout << " " << std::setw(7) << std::setprecision(4) << var_ratio[t];
    

    std::cout << std::endl;
  }

  //    std::cout << "\nlevels:\n";
  //    for(size_t i{0}; i<dec_level.size(); ++i) {
  //        if(dec_level[i] != 0)
  //            std::cout << "level_" << std::setw(3) << std::left << i << " "
  //            << dec_level[i] << std::endl;
  //    }
  //    std::cout << "\nmistakes:\n";
  //    for(size_t i{0}; i<num_mistakes.size(); ++i) {
  //        if(num_mistakes[i] != 0)
  //            std::cout << "b" << std::setw(3) << std::left << i << " " <<
  //            num_mistakes[i] << std::endl;
  //    }

  opt.verbosity = v;
}



// implementation of a scheduling solver
int main(int argc, char *argv[]) {
  auto parser = tempo::getBaseParser();
    
    std::string record_file{""};
    std::string analyse_file{""};
    
    parser.getCmdLine().add<TCLAP::ValueArg<std::string>>(record_file, "", "record",
                                   "record search tree to file", false, "",
                                   "string");
    
    parser.getCmdLine().add<TCLAP::ValueArg<std::string>>(analyse_file, "", "analyse",
                                   "analyse search tree in file", false, "",
                                   "string");

  parser.parse(argc, argv);
  Options opt = parser.getOptions();

  seed(opt.seed);

    if(record_file != "")
        solve(opt, record_file);

    if(analyse_file != "")
        crunch_numbers(opt, analyse_file);
}
