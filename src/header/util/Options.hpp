

#ifndef __CMDLINE_HPP
#define __CMDLINE_HPP

#include <tclap/CmdLine.h>

namespace tempo {

class Options {

public:
  // the actual options
  std::string cmdline; // for reference
  std::string instance_file;

  enum verbosity { SILENT = 0, QUIET, NORMAL, YACKING, SOLVERINFO };
  int verbosity;

  int seed;

  bool print_sol;
  bool print_par;
  bool print_mod;
  bool print_ins;
  bool print_sta;
  bool print_cmd;
  std::string dbg_file;

  bool learning;
  bool edge_finding;
  bool transitivity;

  enum class ChoicePointHeuristics {
    Tightest = 0,
    WeightedDegree,
    WeightedCriticalPath,
    VSIDS
  };
  ChoicePointHeuristics choice_point_heuristics;

  enum class PolarityHeuristic { Identity, LocalExploration, Tightest };

  PolarityHeuristic polarity_heuristic;

  double vsids_decay;
  double vsids_epsilon;

  double restart_factor;
  int restart_base;
  std::string restart_policy;
  std::string input_format;
};

Options parse(int argc, char *argv[]);

static Options no_option;
} // namespace tempo

#endif // __CMDLINE_HPP
