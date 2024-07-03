
/************************************************
 * Tempo Options.hpp
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

#ifndef __CMDLINE_HPP
#define __CMDLINE_HPP

#include <tclap/CmdLine.h>

namespace tempo {

class Options {

public:
    Options() {}
    
  // the actual options
  std::string cmdline{""}; // for reference
  std::string instance_file{"../data/osp/hurley/j7per0-0.txt"};

  enum verbosity { SILENT = 0, QUIET, NORMAL, YACKING, SOLVERINFO };
  int verbosity{verbosity::NORMAL};

  int seed{1};

  int ub{std::numeric_limits<int>::max()};

  bool print_sol{false};
  bool print_par{false};
  bool print_mod{false};
  bool print_ins{false};
  bool print_sta{false};
  bool print_cmd{false};
  std::string dbg_file{""};

  bool learning{true};
  bool edge_finding{true};
  bool transitivity{true};

  bool dichotomy{false};

  bool full_up{false};
  bool order_bound_watch{false};

  enum class ChoicePointHeuristics { Tightest = 0, WeightedDegree, VSIDS };
  ChoicePointHeuristics choice_point_heuristics{ChoicePointHeuristics::VSIDS};

  enum class PolarityHeuristic { Tightest, SolutionGuided, Random };

  double polarity_epsilon{0};

  PolarityHeuristic polarity_heuristic{PolarityHeuristic::Tightest};

  double vsids_decay{0.999};
  double vsids_epsilon{0.05};

  double forgetfulness{0.3};

  //  enum class Minimization { None = 0, Greedy, QuickXplain };

  int minimization{1};

  int greedy_runs{1};

  enum class LiteralScore {
    Size = 0,
    Looseness,
    Activity,
    LoosenessOverActivity
  };

  LiteralScore forget_strategy{3};

  double restart_factor{1.2};
  int restart_base{128};
  std::string restart_policy{"geom"};
  std::string input_format{"osp"};
};

Options parse(int argc, char *argv[]);

static Options no_option;
} // namespace tempo

#endif // __CMDLINE_HPP
