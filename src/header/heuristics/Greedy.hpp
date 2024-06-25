
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

    void addResource(/*DisjunctiveResource<T>& R,*/ std::vector<BooleanVar<>>::iterator bx, std::vector<BooleanVar<>>::iterator ex) {
      //        unscheduled_Intervals_of.resize(unscheduled_Intervals_of.size()+1);
      //        unscheduled_Intervals_of.
      for (auto xi{bx}; xi != ex; ++xi) {
        auto l{solver.boolean.getLiteral(true, *xi)};
        auto c{solver.boolean.getEdge(l)};
        precedences[Interval_map[c.to]].push_back(l);
        precedences[Interval_map[c.from]].push_back(~l);
      }
      //        for(auto& Interval : R) {
      //            unscheduled_Intervals_of.back().add()
      //        }
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
