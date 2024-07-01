/************************************************
 * Tempo Restart.hpp
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

#ifndef __TEMPO_GREEDY_HPP
#define __TEMPO_GREEDY_HPP


#include "Solver.hpp"
#include "Model.hpp"

namespace tempo {

/**********************************************
 * Greedy Primal Heuristic
 **********************************************/

template <typename T> class Greedy  {
public:
    
    Greedy(Solver<T>& s) : solver(s) {
      Interval_map.resize(solver.numeric.size(), -1);
    }

    void addIntervals(std::vector<Interval<T>> &J) {
      Intervals = J;
      unscheduled_Intervals.reserve(Intervals.size());
      unscheduled_Intervals.fill();
      precedences.resize(Intervals.size());
      int i{0};
      for (auto &j : Intervals) {
        Interval_map[j.start.id()] = i;
        Interval_map[j.end.id()] = i;
        ++i;
      }
    }

    void addResource(const std::vector<DisjunctVar<T>>::iterator bx,
                     const std::vector<DisjunctVar<T>>::iterator ex) {

      for (auto xi{bx}; xi != ex; ++xi) {

        auto l{solver.boolean.getLiteral(true, *xi)};
        auto c{solver.boolean.getEdge(l)};

        precedences[Interval_map[c.to]].push_back(l);
        precedences[Interval_map[c.from]].push_back(~l);
      }
    }

    bool runEarliestStart();
    bool runLatestEnd();
    
private:
    Solver<T>& solver;
    std::vector<Interval<T>> Intervals;
    SparseSet<> unscheduled_Intervals;
    std::vector<int> Interval_map;
    std::vector<std::vector<Literal<T>>> precedences;

    //    std::vector<SparseSet<>> unscheduled_Intervals_of;
};


template <typename T>
bool Greedy<T>::runEarliestStart() {
 
    solver.propagate();
    
//    std::cout << solver << std::endl;

    while (not unscheduled_Intervals.empty()) {
      ++solver.num_choicepoints;

      int next{-1};
      for (auto j : unscheduled_Intervals) {
        if (next == -1) {
          next = j;
        } else if (Intervals[next].getEarliestStart(solver) >
                   Intervals[j].getEarliestStart(solver)) {
          next = j;
        } else if (Intervals[next].getEarliestStart(solver) ==
                       Intervals[j].getEarliestStart(solver) and
                   Intervals[next].getLatestEnd(solver) >
                       Intervals[j].getLatestEnd(solver)) {
          next = j;
        } else if (Intervals[next].getEarliestStart(solver) ==
                       Intervals[j].getEarliestStart(solver) and
                   Intervals[next].getLatestEnd(solver) ==
                       Intervals[j].getLatestEnd(solver) and
                   (random() % 2) == 1) {
          next = j;
        }
      }

      try {

        //            std::cout << std::endl << "next="<< Intervals[next] <<
        //            std::endl;

        unscheduled_Intervals.remove_back(next);
        for (auto p : precedences[next]) {
          if (solver.boolean.isUndefined(p.variable())) {
            solver.set(p);
          }
        }
        if (unscheduled_Intervals.backsize() == 1)
          solver.set(Intervals[next].end.before(
              Intervals[next].getEarliestEnd(solver)));

        solver.propagate();

        //            std::cout << solver << std::endl;

      } catch (Failure<T> &f) {
        //            std::cout << "FAILED!\n";
        break;
      }
    }

    bool r{unscheduled_Intervals.empty()};
    unscheduled_Intervals.fill();
    return r;
}



template <typename T>
bool Greedy<T>::runLatestEnd() {
 
    solver.propagate();
    
//    std::cout << solver << std::endl;

    while (not unscheduled_Intervals.empty()) {
      ++solver.num_choicepoints;

      int next{-1};
      for (auto j : unscheduled_Intervals) {
        if (next == -1) {
          next = j;
        } else if (Intervals[next].getLatestEnd(solver) >
                   Intervals[j].getLatestEnd(solver)) {
          next = j;
        } else if (Intervals[next].getLatestEnd(solver) ==
                       Intervals[j].getLatestEnd(solver) and
                   Intervals[next].getEarliestStart(solver) >
                       Intervals[j].getEarliestStart(solver)) {
          next = j;
        } else if (Intervals[next].getEarliestStart(solver) ==
                       Intervals[j].getEarliestStart(solver) and
                   Intervals[next].getLatestEnd(solver) ==
                       Intervals[j].getLatestEnd(solver) and
                   (random() % 2) == 1) {
          next = j;
        }
      }

      try {

        //            std::cout << std::endl << "next="<< Intervals[next] <<
        //            std::endl;

        unscheduled_Intervals.remove_back(next);
        for (auto p : precedences[next]) {
          if (solver.boolean.isUndefined(p.variable())) {
            solver.set(p);
          }
        }
        if (unscheduled_Intervals.backsize() == 1)
          solver.set(Intervals[next].end.before(
              Intervals[next].getEarliestEnd(solver)));

        solver.propagate();

        //            std::cout << solver << std::endl;

      } catch (Failure<T> &f) {
        std::cout << "FAILED!\n";
        break;
      }
    }

    bool r{unscheduled_Intervals.empty()};
    unscheduled_Intervals.fill();
    return r;
}




}

#endif // __GREEDY_HPP
