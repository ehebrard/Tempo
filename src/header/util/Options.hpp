
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
#include <limits>
#include <memory>

#include "enum.hpp"


struct argbase {
    virtual ~argbase() {}
    virtual void assign() = 0;
};

template <typename Opt, typename ClapArg, typename E = void>
struct arg : public argbase {
    ClapArg carg;
    Opt &opt;

    template <typename... T>
    arg(TCLAP::CmdLine &cmd, Opt &opt, T &&...args)
            : carg(std::forward<T>(args)...), opt(opt) {
        cmd.add(carg);
    }

    void assign() override { opt = carg.getValue(); }
};

template<typename Opt, typename ClapArg>
struct arg<Opt, ClapArg, typename std::enable_if<std::is_enum<Opt>{}>::type> : public argbase {
    ClapArg carg;
    Opt &opt;

    template<typename... T>
    arg(TCLAP::CmdLine &cmd, Opt &opt, T &&...args): carg(std::forward<T>(args)...), opt(opt) {
        cmd.add(carg);
    }

    void assign() override {
        opt = static_cast<typename std::remove_reference<Opt>::type>(carg.getValue());
    }
};


struct cmdline {
    TCLAP::CmdLine cmd;
    std::vector<std::unique_ptr<argbase>> args;

    explicit cmdline(const std::string &message, char delimiter = ' ',
                     const std::string &version = "none", bool helpAndVersion = true);

    template <typename ClapArg, typename Opt, typename... T>
    void add(Opt &opt, T &&...clapargs) {
        args.emplace_back(std::move(std::make_unique<arg<Opt, ClapArg>>(
        cmd, opt, std::forward<T>(clapargs)...)));
    }

    void parse(int argc, char *argv[]);
};


namespace tempo {
    namespace detail {
        PENUM(CPH, Tightest, WeightedDegree, VSIDS, Random, LRB, VSIDSHeap)

        PENUM(PH, Tightest, TSG, Random, RSG)
    }

class Options {

public:
    Options() = default;
    
  // the actual options
  std::string cmdline{}; // for reference
  std::string instance_file{"../data/osp/hurley/j7per0-0.txt"};

  enum verbosity { SILENT = 0, QUIET, NORMAL, YACKING, SOLVERINFO };
  int verbosity{verbosity::NORMAL};

  int seed{1};

  int ub{std::numeric_limits<int>::max()};
  int lb{std::numeric_limits<int>::min()};

  bool print_sol{false};
  bool print_par{false};
  bool print_mod{false};
  bool print_ins{false};
  bool print_sta{false};
  bool print_cmd{false};
  std::string dbg_file{""};
//  std::string order_file{""};

  bool learning{true};
  bool edge_finding{true};
    bool overlap_finding{true};
    bool time_tabling{true};
    bool tt_edge_finding{true};
    int incomplete_edge_finding{-std::numeric_limits<int>::max()};
  bool transitivity{true};

  bool dichotomy{false};

  bool full_up{true};
  bool order_bound_watch{false};
    
    bool full_transitivity{false};
    bool primal_boost{true};
    bool ground_update{false};

    using ChoicePointHeuristics = detail::CPH;
    ChoicePointHeuristics choice_point_heuristics{ChoicePointHeuristics::VSIDS};

    using PolarityHeuristic = detail::PH;

    double polarity_epsilon{0};

    PolarityHeuristic polarity_heuristic{PolarityHeuristic::TSG};
    
    double sgd_ratio{.01};

    double vsids_decay{0.999};
    double vsids_epsilon{0.05};
    bool vsids_reasons{false};

    double forgetfulness{0.3};

      enum class Cut { UIP = 0, Booleans, Decisions };
    Cut cut_type{Cut::UIP};
    int minimization{3};
    bool shrinking{false};

    int greedy_runs{10};

    unsigned long search_limit{std::numeric_limits<unsigned long>::max()};
    //    double time_limit;

    enum class LiteralScore {
      Size = 0,
      Looseness,
      Activity,
      LoosenessOverActivity,
        Glue,
        GlueTimesActivity
    };

    LiteralScore forget_strategy{3};

    double restart_factor{1.2};
    int restart_base{128};
    std::string restart_policy{"luby"};
    std::string input_format{"osp"};
};

Options parse(int argc, char *argv[]);

class Parser {
    std::unique_ptr<Options> options;
    std::unique_ptr<cmdline> cmdLine;
public:
    explicit Parser(const std::string &name = "tempo");
    auto getOptions() noexcept -> Options &;
    [[nodiscard]] auto getOptions() const noexcept -> const Options &;
    auto getCmdLine() noexcept -> cmdline &;
    [[nodiscard]] auto getCmdLine() const noexcept -> const cmdline &;
    void parse(int argc, char **argv);
};

auto getBaseParser() -> Parser;


static Options no_option;
} // namespace tempo

#endif // __CMDLINE_HPP
