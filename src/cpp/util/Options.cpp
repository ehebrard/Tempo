
#include <memory>
#include <numeric>

#include <Global.hpp>
#include <util/Options.hpp>

// using namespace schedcl;


tempo::Options tempo::parse(int argc, char *argv[]) {
  using namespace TCLAP;
  using namespace std::string_literals;
  Parser p = getBaseParser();
  p.getOptions().cmdline =
          accumulate(argv, argv + argc, ""s, [&](const std::string& acc, const char *arg) {
              return acc + " " + arg;
          });


  p.parse(argc, argv);
  return p.getOptions();
}

auto tempo::getBaseParser() -> Parser {
    using namespace TCLAP;
    Parser p;
    auto &cmd = p.getCmdLine();
    auto &opt = p.getOptions();
    cmd.add<UnlabeledValueArg<std::string>>(opt.instance_file, "file",
                                            "instance file", true, "", "string");

    cmd.add<ValueArg<int>>(
            opt.verbosity, "", "verbosity",
            "verbosity level (0:silent,1:quiet,2:improvements only,3:verbose", false,
            2, "int");

    cmd.add<ValueArg<int>>(opt.seed, "", "seed", "random seed", false, 1, "int");

    cmd.add<ValueArg<int>>(opt.ub, "", "ub", "initial ub", false,
                           std::numeric_limits<int>::max(), "int");

    cmd.add<ValueArg<int>>(opt.lb, "", "lb", "initial lb", false,
                           std::numeric_limits<int>::min(), "int");

    cmd.add<SwitchArg>(opt.print_sol, "", "print-sol",
                       "print the best found schedule", false);

    cmd.add<SwitchArg>(opt.print_par, "", "print-par", "print the paramters",
                       false);

    cmd.add<SwitchArg>(opt.print_mod, "", "print-mod", "print the model", false);

    cmd.add<SwitchArg>(opt.print_ins, "", "print-ins", "print the instance",
                       false);

    cmd.add<SwitchArg>(opt.print_sta, "", "print-sta", "print the statistics",
                       false);

    cmd.add<SwitchArg>(opt.print_cmd, "", "print-cmd", "print the command-line",
                       false);

    cmd.add<ValueArg<std::string>>(opt.dbg_file, "", "dbg",
                                   "Clause-learning dbg file []", false, "",
                                   "string");

    cmd.add<SwitchArg>(opt.learning, "", "learning", "learn clauses", false);
    cmd.add<SwitchArg>(opt.learning, "", "no-learning",
                       "do not use clause learning", true);

    cmd.add<SwitchArg>(opt.edge_finding, "", "edge-finding", "use edge-finding",
                       false);
    cmd.add<SwitchArg>(opt.edge_finding, "", "no-edge-finding",
                       "do not use edge-finding", true);
    
    cmd.add<SwitchArg>(opt.overlap_finding, "", "overlap-finding", "use overlap-finding",
                       false);
    cmd.add<SwitchArg>(opt.overlap_finding, "", "no-overlap-finding",
                       "do not use overlap-finding", true);
    
    cmd.add<SwitchArg>(opt.tt_edge_finding, "", "tt-ef", "use timetabling reasoning within edge-finding",
                       false);
    cmd.add<SwitchArg>(opt.tt_edge_finding, "", "no-tt-ef",
                       "do not use timetabling reasoning within edge-finding", true);
    
    cmd.add<ValueArg<int>>(opt.incomplete_edge_finding, "", "incomplete-edge-finding", "stop edge-finding at level median - value", false, -std::numeric_limits<int>::max(), "int");
    
    cmd.add<SwitchArg>(opt.time_tabling, "", "time-tabling", "use time-tabling",
                       false);
    cmd.add<SwitchArg>(opt.time_tabling, "", "no-time-tabling",
                       "do not use time-tabling", true);

    cmd.add<SwitchArg>(opt.transitivity, "", "transitivity",
                       "use transitivity reasoning", false);
    cmd.add<SwitchArg>(opt.transitivity, "", "no-transitivity",
                       "do not use transitivity reasoning", true);

    cmd.add<SwitchArg>(opt.dichotomy, "", "dichotomy", "use dichotomic search",
                       false);

    cmd.add<SwitchArg>(opt.full_up, "", "full-up",
                       "unit-propagate bound literals", false);

    cmd.add<SwitchArg>(opt.order_bound_watch, "", "order-watched",
                       "order bound watched lists", false);
    
    cmd.add<SwitchArg>(opt.full_transitivity, "", "full-transitivity",
                       "use full-transitivity", false);

    cmd.add<SwitchArg>(opt.primal_boost, "", "no-primal-boost",
                       "do not use primal boost", true);

    cmd.add<ValueArg<int>>(opt.choice_point_heuristics, "", "cp-heuristic",
                           "type of heuristic used for choice point selection "
                           "(0: Tightest, 1: WDEG, 2: VSIDS (default), 3: Random)",
                           false, 2, "int");
    cmd.add<ValueArg<int>>(
            opt.polarity_heuristic, "", "polarity-heuristic",
            "type "
            "of heuristic used for choice point polarity selection "
            "(0: tightest (default), 1: solution guided, 2: random)",
            false, 0, "int");

    cmd.add<ValueArg<double>>(
            opt.polarity_epsilon, "", "polarity-epsilon",
            "epsilon greedy value for value selection. probability in [0, 1]."
            "0 means that the value selected by the heuristic is always"
            "chosen, 1 means always random. default: 0.01",
            false, 0.01, "double");

    cmd.add<ValueArg<double>>(
            opt.vsids_decay, "", "vsids-decay",
            "decay value for the vsids heuristic, only effective if VSIDS is used "
            "as choice point heuristic",
            false, 0.999, "double");

    cmd.add<ValueArg<double>>(opt.forgetfulness, "", "forgetfulness",
                              "clause base reduction factor (0.3)", false, 0.3,
                              "double");

    cmd.add<ValueArg<int>>(
            opt.minimization, "", "clause-minimization",
            "depth for clause minimization (default 1)",
            false, 1, "int");

    cmd.add<ValueArg<int>>(
            opt.greedy_runs, "", "greedy-runs",
            "number of randomized greedy runs (default 1000)",
            false, 1000, "int");

    cmd.add<ValueArg<unsigned long>>(
        opt.search_limit, "", "search-limit",
        "limit in number of fails (default = a lot)", false, std::numeric_limits<unsigned long>::max(),
        "unsigned long");

    cmd.add<ValueArg<int>>(opt.forget_strategy, "", "forget-strategy",
                           "strategy for clause forgetting "
                           "(0: size (default), 1: literal looseness,"
                           "2: literal activity 3: looseness / activity",
                           false, 3, "int");

    cmd.add<ValueArg<std::string>>(opt.restart_policy, "", "restart",
                                   "choice of restart policy (no, luby, geom)",
                                   false, "geom", "string");

    cmd.add<ValueArg<int>>(
            opt.restart_base, "", "restart-base",
            "base of the sequence for geometric/luby restarts (default=100)", false,
            128, "int");

    cmd.add<ValueArg<double>>(opt.restart_factor, "", "restart-factor",
                              "factor of the geometric sequence (default=1.2)",
                              false, 1.2, "double");

    cmd.add<ValueArg<double>>(opt.vsids_epsilon, "", "vsids-epsilon",
                              "epsilon value for the epsilon greedy vsids "
                              "heuristic, only effective if epsilon greedy "
                              "VSIDS is used as choice point heuristic",
                              false, 0.05, "double");

    cmd.add<ValueArg<std::string>>(
            opt.input_format, "", "input-format",
            "format of input file "
            "(osp: Open Shop (default), jsp: Job Shop "
            "(Lawrence), tsptw: TSP with Time Windows, jstl: JSP with Time Lags)",
            false, "osp", "string");
    return p;
}

cmdline::cmdline(const std::string &message, const char delimiter, const std::string &version, bool helpAndVersion)
        : cmd(message, delimiter, version, helpAndVersion) {}

void cmdline::parse(int argc, char **argv) {
    cmd.parse(argc, argv);
    for (auto &arg: args) {
        arg->assign();
    }
}

tempo::Parser::Parser(const std::string &name) : options(std::make_unique<Options>()),
                                                 cmdLine(std::make_unique<cmdline>(name, ' ')) {}

auto tempo::Parser::getOptions() noexcept -> tempo::Options & {
    return *options;
}

auto tempo::Parser::getOptions() const noexcept -> const tempo::Options & {
    return *options;
}

auto tempo::Parser::getCmdLine() noexcept -> cmdline & {
    return *cmdLine;
}

auto tempo::Parser::getCmdLine() const noexcept -> const cmdline & {
    return *cmdLine;
}

void tempo::Parser::parse(int argc, char **argv) {
    cmdLine->parse(argc, argv);
}
