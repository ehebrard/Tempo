//
// Created by tim on 15.11.22.
//

#ifndef TEMPO_VSIDS_HPP
#define TEMPO_VSIDS_HPP

#include "GapOverActivity.hpp"

namespace tempo::heuristics {
    /**
     * @brief VSIDS (Variable State Independent Decaying Sum) heuristic
     * @details @copybrief
     * Selects variables which occur in the largest number of clauses (proportionally to highest activity). Activity for
     * a literal l is calculated the following way:
     * \f$ A_l = \sum_{i=1}^n \gamma^{n - i} w_i \f$ where \f$ \gamma \f$ is a constant decay factor and \f$ n \f$ is
     * the the number of activity updates
     * The final score of a literal is inversely proportional to its activity
     */
    class VSIDS : public GapOverActivity {
    public:


        /**
         * CTor. Infers parameters from given scheduler. Subscribes to events of interest
         * @tparam T timing type
         * @param solver target solver
         */
      template <concepts::scalar T>
      explicit VSIDS(Solver<T> &solver)
          : GapOverActivity(solver, solver.ConflictExtracted.subscribe_handled(
                                        [this](const auto &arg) {
                                          this->updateActivity(arg);
                                        })) {}

      template <concepts::scalar T>
      void updateActivity(const Solver<T> &solver) {
#ifdef OLDVSIDS
        GapOverActivity::activity.update(solver);
#else
        bool_buffer.clear();
        num_buffer.clear();
          
//          std::vector<Literal<T>> lit_buffer;
          
        for (auto l : solver.lastLearnt()) {
          if (l.isNumeric()) {
            num_buffer.push_back(l.variable());
          } else {
            bool_buffer.push_back(l.variable());
          }
//            auto exp{solver.getReason(l)};
//            exp.explain(l, lit_buffer);
        }
          
//          if(solver.getOptions())
//          for (auto l : lit_buffer) {
//              if (l.isNumeric()) {
//                num_buffer.push_back(l.variable());
//              } else {
//                bool_buffer.push_back(l.variable());
//              }
//          }

        for (auto i : solver.cut.cached_) {
          auto l{solver.getLiteral(i)};
          if (l.isNumeric()) {
            num_buffer.push_back(l.variable());
          } else {
            bool_buffer.push_back(l.variable());
          }
        }

        GapOverActivity::numeric_activity.update(num_buffer);
        GapOverActivity::boolean_activity.update(bool_buffer);
#endif
      }

      std::vector<var_t> bool_buffer;
      std::vector<var_t> num_buffer;
    };
}

#endif //TEMPO_VSIDS_HPP
