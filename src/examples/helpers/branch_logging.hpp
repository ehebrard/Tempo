/**
 * @author Emmanuel Hebrard
 * @date 18.03.25
 * @file branch_logging.hpp
 * @brief Runs the solver and records branching decision for further analysis
 */

#ifndef BRANCH_LOGGING_HPP
#define BRANCH_LOGGING_HPP

#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Solver.hpp"
#include "util/traits.hpp"

namespace tempo {
/**
 * Runs the solver and records branching decision for further analysis
 * @tparam T Timing type
 * @param S ready-to-run solver
 * @param schedule schedule interval
 * @param record_file destination record file
 * @return makespan
 */
template <concepts::scalar T>
int solve(Solver<T> &S, const Interval<T> &schedule,
          const std::string &record_file) {
  const auto &opt = S.getOptions();
  if (opt.restart_policy != "no" or opt.primal_boost) {
    throw std::runtime_error(
        "please set restart policy to 'no' and disable primal boost");
  }

  long unsigned int num_correct_decisions{0};
  long unsigned int num_wrong_decisions{0};
  long unsigned int cp_in_wasted_restarts{0};
  long unsigned int previous_total_cp{0};

  std::vector<std::vector<Literal<int>>> right_branches;
  std::stringstream buffer;
  std::ofstream outfile(record_file);
  auto solution_count{0};

#ifdef OLD
  SubscriberHandle solutionHandler(
      S.SolutionFound.subscribe_handled([&](const auto &) {
        ++solution_count;

        num_correct_decisions += S.numDecision();
        unsigned long num_wrong{0};
        for (auto &branches : right_branches)
          num_wrong += branches.size();
        num_wrong_decisions += num_wrong;
        buffer << S.numeric.solutionLower(schedule.duration) << " "
               << S.num_choicepoints - cp_in_wasted_restarts << " "
               << (right_branches.size() > S.numDecision()) << " "
               << S.numDecision() << " " << num_wrong;
        for (unsigned i{0}; i < S.numDecision(); ++i) {
          if (right_branches.size() > i) {
            buffer << " " << right_branches[i].size();
            for (auto l : right_branches[i]) {
              buffer << " " << l.sign() << " " << l.variable();
            }
          } else {
            buffer << " 0";
          }

          assert(not S.getDecisions()[i].isNumeric());

          buffer << " " << S.getDecisions()[i].sign() << " "
                 << S.getDecisions()[i].variable();
        }
        if (right_branches.size() > S.numDecision()) {
          buffer << " " << right_branches.back().size();
          for (auto l : right_branches.back()) {
            buffer << " " << l.sign() << " " << l.variable();
          }
        }
#else
  SubscriberHandle solutionHandler(
      S.SolutionFound.subscribe_handled([&](const auto &) {
        ++solution_count;

        num_correct_decisions += S.numDecision();
        unsigned long num_wrong{0};
        for (auto &branches : right_branches)
          num_wrong += branches.size();
        num_wrong_decisions += num_wrong;
        buffer << S.numeric.solutionLower(schedule.duration) << " "
               << S.num_choicepoints - cp_in_wasted_restarts << " "
               << (right_branches.size() > S.numDecision()) << " "
               << S.numDecision() << " " << num_wrong;
        for (unsigned i{0}; i < S.numDecision(); ++i) {
          if (right_branches.size() > i) {
            buffer << " " << right_branches[i].size();
            for (auto l : right_branches[i]) {
              buffer << " " << l.isNumeric() << " " << l.sign() << " "
                     << l.variable();
              if (l.isNumeric())
                buffer << " " << l.value();
              else
                buffer << " " << l.semantic();
            }
          } else {
            buffer << " 0";
          }

          assert(not S.getDecisions()[i].isNumeric());

          buffer << " 0 " << S.getDecisions()[i].sign() << " "
                 << S.getDecisions()[i].variable() << " "
                 << S.getDecisions()[i].semantic();
        }
        if (right_branches.size() > S.numDecision()) {
          buffer << " " << right_branches.back().size();
          for (auto l : right_branches.back()) {
            buffer << " " << l.isNumeric() << " " << l.sign() << " "
                   << l.variable();
            if (l.isNumeric())
              buffer << " " << l.value();
            else
              buffer << " " << l.semantic();
          }
        }
#endif

        buffer << std::endl;
        right_branches.clear();

#ifdef VERBOSE
        std::cout << "solution #" << solution_count << std::endl;

#endif
      }));

  SubscriberHandle failCLHandler(
      S.ClauseAdded.subscribe_handled([&](const auto &solver) {
        right_branches.resize(S.numDecision() + 1);
        right_branches.back().push_back(solver.lastLearnt()[0]);

#ifdef VERBOSE
        std::cout << "\nfail\n";
        for (size_t i{0}; i < S.numDecision(); ++i) {
          std::cout << S.getDecisions()[i];
          for (size_t j{0}; j < right_branches[i].size(); ++j) {
            std::cout << " " << right_branches[i][j];
          }
          std::cout << std::endl;
        }
        std::cout << "F";
        for (size_t j{0}; j < right_branches.back().size(); ++j) {
          std::cout << " " << right_branches.back()[j];
        }
        std::cout << std::endl;
#endif
      }));

  SubscriberHandle failNOCLHandler(
      S.DeductionMade.subscribe_handled([&](const auto &lit) {
        right_branches.resize(S.numDecision());
        right_branches.back().push_back(lit);

#ifdef VERBOSE
        std::cout << "\nfail\n";
        for (size_t i{0}; i < S.numDecision() - 1; ++i) {
          std::cout << S.getDecisions()[i];
          for (size_t j{0}; j < right_branches[i].size(); ++j) {
            std::cout << " " << right_branches[i][j];
          }
          std::cout << std::endl;
        }
        std::cout << "F";
        for (size_t j{0}; j < right_branches.back().size(); ++j) {
          std::cout << " " << right_branches.back()[j];
        }
        std::cout << std::endl;
#endif
      }));

  SubscriberHandle restartHandler(
      S.SearchRestarted.subscribe_handled([&](const auto on_solution) {
        if (on_solution) {
          previous_total_cp = S.num_choicepoints;
        } else {
          cp_in_wasted_restarts += (S.num_choicepoints - previous_total_cp);
        }
      }));

  if (opt.lb > -Constant::Infinity<int>) {
    S.post(schedule.duration >= opt.lb);
  }

  S.minimize(schedule.duration);
  if (not S.numeric.hasSolution()) {
    throw std::runtime_error("no solution found yet");
  }

  auto obj{S.numeric.solutionLower(schedule.duration)};
  outfile << obj << std::endl << buffer.str();
  return obj;
}
} // namespace tempo

#endif // BRANCH_LOGGING_HPP
