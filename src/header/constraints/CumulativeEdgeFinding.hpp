/************************************************
 * Tempo CumulativeEdgeFinding.hpp
 * Implementation of the "strong edge-finding" algorithm as described in
 * Vincent Gingras and Claude-Guy Quimper. 2016. Generalizing the edge-finder
 *rule for the cumulative constraint. In Proceedings of the Twenty-Fifth
 *International Joint Conference on Artificial Intelligence (IJCAI'16). AAAI
 *Press, 3103–3109.
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

#ifndef TEMPO_CUMULATIVEEDGEFINDING_HPP
#define TEMPO_CUMULATIVEEDGEFINDING_HPP

#include <cassert>
#include <map>
#include <vector>
#include <sstream>
#include <numeric>


#include "Explanation.hpp"
#include "Global.hpp"
#include "Model.hpp"
#include "constraints/Constraint.hpp"
#include "util/List.hpp"
namespace tempo {

template <typename T = int> struct Timepoint {

  Timepoint() {}
  Timepoint(T time, T increment, T incrementMax)
      : time(time), increment(increment),
        incrementMax(incrementMax) {}

  T time{0};
  T increment{0};
  T incrementMax{0};
    
    void merge(const Timepoint<T>& t) {
        assert(t.time == time);
        
        increment += t.increment;
        incrementMax += t.incrementMax;
    }

  std::ostream &display(std::ostream &os) const {
    os << "(t=" << time << "|∂=" << increment << "|∂max=" << incrementMax << ")";
    return os;
  }
};

template <typename T>
std::ostream &operator<<(std::ostream &os, const Timepoint<T> &x) {
  return x.display(os);
}

template<typename T>
class Solver;

template <typename T> class CumulativeEdgeFinding : public Constraint<T> {
private:
  Solver<T> &solver;
  NumericVar<T> capacity;
  Interval<T> schedule;
  std::vector<Interval<T>> task;
  std::vector<NumericVar<T>> demand;

  // mapping from the tasks that need to be adjusted to the corresponding
  // left-cut interval
  std::vector<int> prec;

  // mapping the lower bound of the set explaining the adjustment
  std::vector<int> contact;

  // tasks that have been found to be in conflict with a left-cut interval
  SparseSet<> in_conflict;

  List<Timepoint<T>> profile;
  std::vector<int>
      est_; // pointer to the profile element used for est of the task
  std::vector<int>
      ect_; // pointer to the profile element used for ect of the task
  std::vector<int>
      lct_; // pointer to the profile element used for lct of the task

  std::vector<int> est_shared; // points to the possibly shared element whose
                               // t is the est of the task
  std::vector<int> ect_shared; // points to the possibly shared element whose
                               // t is the ect of the task
  std::vector<int> lct_shared; // points to the possibly shared element whose
                               // t is the lct of the task

  size_t leftcut_pointer;

  int &get_shared_pointer(const int evt) {
    auto evti{evt - 1};
    auto i{(evti) / 3};
    switch ((evti % 3)) {
    case 0:
      return est_shared[i];
    case 1:
      return ect_shared[i];
    case 2:
      return lct_shared[i];
    }
  }

    int get_unique_pointer(const int evt) const {
    auto evti{evt - 1};
    auto i{(evti) / 3};
    switch ((evti % 3)) {
    case 0:
      return est_[i];
    case 1:
      return ect_[i];
    case 2:
      return lct_[i];
    }
  }

  std::vector<int> event_ordering;
  std::vector<int> lct_order;
  int alpha;
  int beta;

  std::vector<std::vector<Literal<T>>> explanation;
  Reversible<size_t> num_explanations;

public:
  template <typename ItTask, typename ItNVar>
  CumulativeEdgeFinding(Solver<T> &solver, const Interval<T> sched,
                        const NumericVar<T> capacity, const ItTask beg_task,
                        const ItTask end_task, const ItNVar beg_dem);
  virtual ~CumulativeEdgeFinding();

  // helpers
  T est(const unsigned i) const;
  T lst(const unsigned i) const;
  T ect(const unsigned i) const;
  T lct(const unsigned i) const;
  T minduration(const unsigned i) const;
  T maxduration(const unsigned i) const;
  T mindemand(const unsigned i) const;
  T maxdemand(const unsigned i) const;
  //    T energy(const unsigned i) const;
  //    bool hasFixedPart(const unsigned i) const;
  //    T overlapedFixedPartEnergy(const unsigned i, const unsigned j) const;
  //    std::vector<T> overlapedFixedPartEnergy();

      bool notify(const Literal<T>, const int rank) override;
      void post(const int idx) override;
      void propagate() override;
    
    void forwardDetection();
    
  //    void computeBound(const int i);

  //    /* Tools for the adjustments */
  //    T scheduleOmegaDetection(
  //                             const T C, const int i
  //    //      ,const T max_lct
  //    ); // scheduling omega used for the detection
  //    T scheduleOmegaAdjustment(
  //                              const T C, const int i
  //    //      ,const T max_lct
  //    ); // scheduling omega used for the adjustment
  //    int computeEstPrime(
  //                        const int i); // end of the scheduling in the upper
  //                        part of the profile
  //    void computeMaximumOverflow(
  //                                const int i); // energy to spend at the
  //                                upper part of the profile
  //    std::vector<T>
  //    consOverSlackAvail(const int time); // retrive the attributes of the
  //    time
  //    // points not belong to the profile

  // initialise the profile with every tasks
  void initialiseProfile();

  // set the current profile to contain the tasks lct_order[0],..,lct_order[i]
    void setLeftCutToTask(const int i);
    void shrinkLeftCutToTime(const T t);
    void rmTask(const int i);
    void addTask(const int i);

      // function used in explanation
      void xplain(const Literal<T> l, const hint h,
                  std::vector<Literal<T>> &Cl) override;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

  std::string prettyTask(const int i) const;
  std::string asciiArt(const int i) const;
};

template <typename T>
std::string CumulativeEdgeFinding<T>::prettyTask(const int i) const {
  std::stringstream ss;
  ss << "t" << task[i].id() << ": [" << est(i) << ".." << lct(i) << "] (" << mindemand(i)
     << "x" << minduration(i) << ")";
  return ss.str();
}

template <typename T>
std::string CumulativeEdgeFinding<T>::asciiArt(const int i) const {
  std::stringstream ss;
  ss << std::setw(3) << std::right << mindemand(i) << "x" << std::setw(3)
     << std::left << minduration(i) << " " << std::right;
  for (auto k{0}; k < est(i); ++k) {
    ss << " ";
  }
  ss << "[";
  for (auto k{est(i) + 1}; k < ect(i); ++k) {
    ss << "=";
  }
  if (lct(i) == Constant::Infinity<T>) {
    ss << "=... " << est(i) << "...";
  } else {
    for (auto k{ect(i)}; k < lct(i); ++k) {
      ss << ".";
    }
    ss << "] " << est(i) << ".." << lct(i);
  }
  return ss.str();
}

template <typename T> T CumulativeEdgeFinding<T>::est(const unsigned i) const {
  return task[i].getEarliestStart(solver);
}

template <typename T> T CumulativeEdgeFinding<T>::lst(const unsigned i) const {
  return task[i].getLatestStart(solver);
}

template <typename T> T CumulativeEdgeFinding<T>::ect(const unsigned i) const {
  return task[i].getEarliestEnd(solver);
}

template <typename T> T CumulativeEdgeFinding<T>::lct(const unsigned i) const {
  return task[i].getLatestEnd(solver);
}

template <typename T>
T CumulativeEdgeFinding<T>::minduration(const unsigned i) const {
  return task[i].minDuration(solver);
}

template <typename T>
T CumulativeEdgeFinding<T>::maxduration(const unsigned i) const {
  return task[i].maxDuration(solver);
}

template <typename T>
T CumulativeEdgeFinding<T>::mindemand(const unsigned i) const {
  return demand[i].min(solver);
}

template <typename T>
T CumulativeEdgeFinding<T>::maxdemand(const unsigned i) const {
  return demand[i].max(solver);
}

//template <typename T>
//T CumulativeEdgeFinding<T>::energy(const unsigned i) const {
//  return mindemand(i) * minduration(i);
//}
//
//template <typename T>
//bool CumulativeEdgeFinding<T>::hasFixedPart(const unsigned i) const {
//  return lst(i) < ect(i);
//}

// template <typename T>
// T CumulativeEdgeFinding<T>::overlapedFixedPartEnergy(const unsigned i, const
// unsigned j) const {
//     auto energy{0};
//     if (lct(i) > lct(j) and hasFixedPart(i) and lst(i) < lct(j)) {
//         energy = (std::max(0, std::min(ect(i), lct(j)) - lst(i))) *
//         mindemand(i);
//     }
//     return energy;
// }

// template <typename T>
// std::vector<T> CumulativeEdgeFinding<T>::overlapedFixedPartEnergy() {
//     std::vector<T> energy;
//     energy.resize(task.size());
//     for (auto j : lct_order) {
//         auto E{0};
//         for (auto i = j + 1; i < task.size(); ++i) {
//             E += overlapedFixedPartEnergy(i, j);
//         }
//         energy[j] = E;
//         while (j + 1 < task.size() and lct(j + 1) == lct(j)) {
//             energy[j + 1] = energy[j];
//             ++j;
//         }
//     }
//     return energy;
// }

template <typename T>
template <typename ItTask, typename ItNVar>
CumulativeEdgeFinding<T>::CumulativeEdgeFinding(
    Solver<T> &solver, const Interval<T> sched, const NumericVar<T> cap,
    const ItTask beg_task, const ItTask end_task, const ItNVar beg_dem)
    : solver(solver), num_explanations(0, &(solver.getEnv())) {
  schedule = sched, capacity = cap;

  Constraint<T>::priority = Priority::Low;

  auto dp{beg_dem};
  for (auto jp{beg_task}; jp != end_task; ++jp) {
    task.push_back(*jp);
    demand.push_back(*dp);
    ++dp;
  }

  auto ip{task.size()};

  in_conflict.reserve(ip);
  prec.resize(ip);
  est_.resize(ip);
  ect_.resize(ip);
  lct_.resize(ip);

  est_shared.resize(ip);
  ect_shared.resize(ip);
  lct_shared.resize(ip);

  contact.resize(ip);
  //    contact_time.resize(ip);
  //    minEct.resize(ip);
  //    isFeasible.resize(ip);
  //    maxOverflow.resize(ip);
  //    tp_attributes_est_i.resize(3);
  //    tp_attributes_ect_i.resize(3);

  auto maxcap{capacity.max(solver)};

  for (unsigned i = 0; i < ip; ++i) {
    est_[i] = profile.create_element();
    ect_[i] = profile.create_element();
    lct_[i] = profile.create_element();
  }

  lct_order.resize(task.size());
  std::iota(lct_order.begin(), lct_order.end(), 0);
  event_ordering.resize(profile.size());
  std::iota(event_ordering.begin(), event_ordering.end(), 1);
}

template <typename T> CumulativeEdgeFinding<T>::~CumulativeEdgeFinding() {}


template <typename T> void CumulativeEdgeFinding<T>::post(const int idx) {

  Constraint<T>::cons_id = idx;
  Constraint<T>::idempotent = true;

#ifdef DBG_CEDGEFINDING
  if (DBG_CEDGEFINDING) {
    std::cout << "post " << *this << std::endl;
  }
#endif

  for (size_t i{0}; i < task.size(); ++i) {
    solver.wake_me_on(lb<T>(task[i].end.id()), this->id());
    solver.wake_me_on(ub<T>(task[i].end.id()), this->id());
    solver.wake_me_on(lb<T>(demand[i].id()), this->id());
  }
}

template <typename T>
bool CumulativeEdgeFinding<T>::notify(const Literal<T>, const int) {
  return true;
}

template <typename T> void CumulativeEdgeFinding<T>::initialiseProfile() {

  profile.clear();

  auto C{capacity.max(solver)};

  for (auto i : lct_order) {
      
#ifdef DBG_SEF
      if (DBG_SEF) {
          std::cout << " init timepoints for " << prettyTask(i) << std::endl;
      }
#endif
      
      est_shared[i] = est_[i];
      ect_shared[i] = ect_[i];
      lct_shared[i] = lct_[i];
      
    profile[est_[i]].time = est(i);
    profile[ect_[i]].time = ect(i);
    profile[lct_[i]].time = lct(i);
    profile[est_[i]].increment = mindemand(i);
    profile[est_[i]].incrementMax = mindemand(i);
    profile[ect_[i]].increment = -mindemand(i);
    profile[lct_[i]].incrementMax = -mindemand(i);
  }
    
//    for(auto evt : event_ordering) {
//        std::cout << evt << ": " << profile[evt] << " " << get_unique_pointer
//    }

  std::sort(event_ordering.begin(), event_ordering.end(),
            [this](const int i, const int j) {
              return this->profile[i].time < this->profile[j].time;
            });

    
#ifdef DBG_SEF
      if (DBG_SEF) {
          std::cout << " remove repetitions " << std::endl;
      }
#endif
    
  profile.add_front(event_ordering[0]);
    auto previous{event_ordering[0]};
    
  for (unsigned i{1}; i < event_ordering.size(); ++i) {
    auto current{event_ordering[i]};
      
//      std::cout << "compare " << profile[previous] << " <> " << profile[current] << std::endl;
      
    if (profile[previous].time == profile[current].time) {
      profile[previous].merge(profile[current]);
        get_shared_pointer(current) = previous; //get_unique_pointer(previous);
        
//        std::cout << "merge\n";
    } else {
        get_shared_pointer(current) = current; //get_unique_pointer(current);
      profile.add_after(previous, current);
      previous = current;
        
//        std::cout << "new\n";
    }
  }
    
    
//    std::cout << profile << std::endl;
    
    for(auto i : lct_order) {
        
//        std::cout << prettyTask(i) << ": " << est(i) << " / " << est_shared[i] << " / " << profile[est_shared[i]].time << std::endl;
//        
        assert(profile[est_shared[i]].time == est(i));
        assert(profile[ect_shared[i]].time == ect(i));
        assert(profile[lct_shared[i]].time == lct(i));
    }
    
    

  leftcut_pointer = static_cast<int>(lct_order.size());
}


template <typename T> void CumulativeEdgeFinding<T>::rmTask(const int i) {
    profile[est_shared[i]].increment -= mindemand(i);
    profile[est_shared[i]].incrementMax -= mindemand(i);
    profile[ect_shared[i]].increment += mindemand(i);
    profile[lct_shared[i]].incrementMax += mindemand(i);
}

template <typename T> void CumulativeEdgeFinding<T>::addTask(const int i) {
    profile[est_shared[i]].increment += mindemand(i);
    profile[est_shared[i]].incrementMax += mindemand(i);
    profile[ect_shared[i]].increment -= mindemand(i);
    profile[lct_shared[i]].incrementMax -= mindemand(i);
}

template <typename T> void CumulativeEdgeFinding<T>::setLeftCutToTask(const int i) {
    if(i < leftcut_pointer) {
        // we need to remove the tasks in the interval [i, leftcut_pointer)
        while(leftcut_pointer --> i) {
            rmTask(leftcut_pointer);
        }
    } else {
        // we need to add the tasks in the interval [leftcut_pointer, i)
        while(leftcut_pointer > i) {
            addTask(leftcut_pointer++);
        }
    }
}

template <typename T> void CumulativeEdgeFinding<T>::shrinkLeftCutToTime(const T t) {
        // we need to remove the tasks that have a lct equal to or larger than t
        while(t <= lct(lct_order[leftcut_pointer])) {
            rmTask(leftcut_pointer--);
        }
}

template <typename T> void CumulativeEdgeFinding<T>::propagate() {

  in_conflict.clear();

  std::sort(lct_order.begin(), lct_order.end(),
            [this](const int i, const int j) { return lct(i) < lct(j); });

#ifdef DBG_SEF
  if (DBG_SEF) {
    std::cout << "\n\nstart propagation (" << capacity.max(solver) << ")\n";
    for (auto j : lct_order) {
      std::cout << "task " << j << ": " << asciiArt(j) << std::endl;
    }
  }
#endif

  initialiseProfile();
    
    
    forwardDetection();

}



template <typename T> void CumulativeEdgeFinding<T>::forwardDetection() {
    
#ifdef DBG_SEF
    if (DBG_SEF) {
        std::cout << "\n\nstart forward detection\n";
        for (auto j : lct_order) {
            std::cout << "task " << j << ": " << asciiArt(j) << std::endl;
        }
    }
#endif
    
    //    std::cout << " first = " << *(lct_order.rend()) << std::endl;
    
    auto cap{capacity.max(solver)};
    auto stop{lct_order.rend()};
    --stop;
    
    // explore the tasks by decreasing lct
    for (auto ii{lct_order.rbegin()}; ii != stop;) {
        
#ifdef DBG_SEF
        if (DBG_SEF) {
            std::cout << " - analyse tasks whose lct is " << lct(*ii) << std::endl;
        }
#endif
        
        // remove tasks whose lct is larger than or equal to lct(*ii)
        shrinkLeftCutToTime(lct(*ii));
        
        
        
        // if there are no more tasks, all those in the current level have the same
        // lct and we can stop
        if (ii == lct_order.rend()) {
            break;
        }
        
        // lct of the task that precedes all the task whose lct is lct(i)
        auto j{*ii};
        auto lct_j{lct(j)};
        
        // tasks in the range [leftcut_pointer, *ii) all have a lct equal to lct(*ii)
        // - add their "prime" versions one by one and run scheduleOmega
        do {
            
#ifdef DBG_SEF
            if (DBG_SEF) {
                std::cout << " - analyse task " << *ii << ": schedule tasks on [0," << lct_j << ")\n";
            }
#endif
            
        } while(*(ii--) != leftcut_pointer);
        
    }
}
      
//      
//    while (*ii !=) {
//
//      auto i{*is};
//
//      if (ect(i) < lct(i)) {
//        // std::cout<< profile[est_[task.size()]] << std::endl;
//        for (auto p{profile.begin()}; p != profile.end(); ++p) {
//          p->capacity = cap;
//        }
//        auto omega_ect{0};
//        if (addPrime(i)) {
//          omega_ect = scheduleOmegaDetection(
//              cap, i
//              //                                             , lct_j
//          );
//          rmPrime();
//#ifdef DBG_SEF
//          if (DBG_SEF) {
//            std::cout << "  * rm " << i << "'\n";
//          }
//#endif
//
//#ifdef DBG_SEF
//          if (DBG_SEF) {
//            std::cout << " ect^H = " << omega_ect << " / lct(S) = " << lct_j
//                      << std::endl;
//          }
//#endif
//
//          if (omega_ect > lct_j) {
//
//#ifdef DBG_SEF
//            if (DBG_SEF) {
//              std::cout << " task " << i
//                        << " is in conflict with task interval "
//                        << *(lct_order.begin()) << ".." << j << std::endl;
//            }
//#endif
//            prec[i] = j;
//            in_conflict.add(i);
//            // rmPrime();
//          }
//
//          else {
//
//#ifdef DBG_SEF
//            if (DBG_SEF) {
//              std::cout << " compute bound\n";
//            }
//#endif
//
//            auto ti{ii};
//            computeBound(i);
//
//            if (beta != -1) {
//#ifdef DBG_SEF
//              if (DBG_SEF) {
//                std::cout << "  - beta = " << beta << std::endl;
//              }
//#endif
//              for (; *ti != beta; ++ti) {
//                rmTask(*ti);
//                ++lc_ptr;
//              }
//              for (auto p{profile.begin()}; p != profile.end(); ++p) {
//                p->capacity = cap;
//              }
//              auto ect_i_H{0};
//              if (addPrime(i)) {
//                ect_i_H = scheduleOmegaDetection(cap, i
//                                                 //, lct(beta)
//                );
//                rmPrime();
//#ifdef DBG_SEF
//                if (DBG_SEF) {
//                  std::cout << "  * rm " << i << "' (beta)\n";
//                }
//#endif
//              }
//
//              if (ect_i_H > lct(beta)) {
//#ifdef DBG_SEF
//                if (DBG_SEF) {
//                  std::cout
//                      << " task " << i << " is in conflict with task interval "
//                      << *(lct_order.begin()) << ".." << beta << std::endl;
//                }
//#endif
//                prec[i] = beta;
//                in_conflict.add(i);
//              }
//            }
//            if (prec[i] == -1 and alpha != -1) {
//
//#ifdef DBG_SEF
//              if (DBG_SEF) {
//                std::cout << "  - alpha = " << alpha << std::endl;
//              }
//#endif
//              for (; *ti != alpha; ++ti) {
//                rmTask(*ti);
//                ++lc_ptr;
//              }
//              for (auto p{profile.begin()}; p != profile.end(); ++p) {
//                p->capacity = cap;
//              }
//              auto ect_i_H{0};
//              if (addPrime(i)) {
//                ect_i_H = scheduleOmegaDetection(
//                    cap, i
//                    //                                                 ,
//                    //                                                 lct(alpha)
//                );
//                rmPrime();
//#ifdef DBG_SEF
//                if (DBG_SEF) {
//                  std::cout << "  * rm " << i << "' (alpha)\n";
//                }
//#endif
//              }
//
//              if (ect_i_H > lct(alpha)) {
//
//#ifdef DBG_SEF
//                if (DBG_SEF) {
//                  std::cout
//                      << " task " << i << " is in conflict with task interval "
//                      << *(lct_order.begin()) << ".." << alpha << std::endl;
//                }
//#endif
//
//                prec[i] = alpha;
//                in_conflict.add(i);
//              }
//            }
//            for (; ti != ii;) {
//              --ti;
//              --lc_ptr;
//              addTask(*ti);
//            }
//          }
//        }
//      }
//      ++is;
//    }
//  }
//
//#ifdef DBG_SEF
//  if (DBG_SEF) {
//    std::cout << "end forward detection\n";
//  }
//#endif
//}
//
//template <typename T>
//void CumulativeEdgeFinding<T>::computeMaximumOverflow(const int i) {
//
//#ifdef DBG_SEF
//  if (DBG_SEF) {
//    std::cout << "compute max overflow of " << i << " with prec[" << i
//              << "]=" << prec[i] << std::endl;
//    std::cout << " contact = " << contact[i] << " (" << contact_time[i]
//              << ")\n";
//  }
//#endif
//
//  auto e_idx{num_explanations};
//  if (explanation.size() <= e_idx) {
//    explanation.resize(e_idx + 1);
//  } else {
//    explanation[e_idx].omega.clear();
//  }
//  auto ECT{Constant::Infinity<T>};
//  auto minEnergy{Constant::Infinity<T>};
//  auto lb{Constant::Infinity<T>};
//  auto k{prec[i]};
//
//  for (auto j : lct_order) {
//    if (lct(j) > lct(k))
//      break;
//    if (ect(j) >= contact_time[i] and lct(j) <= lct(k)) {
//      ECT = std::min(ect(j), ECT);
//      minEnergy = std::min(energy(j), minEnergy);
//      lb = std::min(lb, est(j));
//
//      explanation[e_idx].omega.push_back(j);
//    }
//  }
//
//  explanation[e_idx].lb_omega = lb;
//  explanation[e_idx].ub_omega = lct(k);
//  explanation[e_idx].i = i;
//  explanation[e_idx].lb_i = std::min(est(i), contact_time[i]);
//
//  minEct[i] = ECT;
//
//  auto cont{contact_time[i]};
//  auto ip{static_cast<int>(task.size())};
//
//  if (ect(i) < lct(k)) {
//    // std::cout << "i:" << i << " // " << cont << " <> " << est(i) <<
//    // std::endl;
//
//    /*if (cont < est(i)) {
//    cont = est(i);
//    }*/
//    //      std::cout << "consOverSlackAvail(" << cont << ")\n";
//    auto overlap1{0};
//    auto slackUnder1{0};
//    auto available1{0};
//
//    if (cont > est(i)) {
//      auto t{profile[contact[i]]};
//      overlap1 = t.overlap;
//      slackUnder1 = t.slackUnder;
//      available1 = t.available;
//    } else {
//      auto t{profile[est_[ip]]};
//      overlap1 = t.overlap;
//      slackUnder1 = t.slackUnder;
//      available1 = t.available;
//    }
//
//    //      std::cout << "consOverSlackAvail(" << ect(i) << ")\n";
//    // auto attributes2{consOverSlackAvail(ect(i))};
//    // auto overlap2{tp_attributes_ect_i[0]};
//    auto t{profile[ect_[ip]]};
//    auto overlap2 = t.overlap;
//
//    //      std::cout << "consOverSlackAvail(" << lct(k) << ")\n";
//    // auto attributes3{consOverSlackAvail(lct(k))};
//    t = profile[lct_[k]];
//    auto slackUnder2{t.slackUnder};
//    auto overlap3{t.overlap};
//    auto available2{t.available};
//    if (available2 - available1 < energy(i)) {
//      if (isFeasible[i] and minEnergy > slackUnder2 - slackUnder1) {
//        maxOverflow[i] = overlap3 - overlap1;
//      } else {
//        maxOverflow[i] = overlap3 - overlap1 - (slackUnder2 - slackUnder1);
//      }
//    } else {
//      if (isFeasible[i] and minEnergy > slackUnder2 - slackUnder1) {
//        maxOverflow[i] = overlap2 - overlap1;
//      } else {
//        maxOverflow[i] = overlap2 - overlap1 - (slackUnder2 - slackUnder1);
//      }
//    }
//  } else {
//    auto overlap1{0};
//    auto slackUnder1{0};
//    if (cont > est(i)) {
//      auto t{profile[contact[i]]};
//      overlap1 = t.overlap;
//      slackUnder1 = t.slackUnder;
//    } else {
//      auto t{profile[est_[ip]]};
//      overlap1 = t.overlap;
//      slackUnder1 = t.slackUnder;
//    }
//    auto t = profile[lct_[k]];
//    auto slackUnder2{t.slackUnder};
//    auto overlap2{t.overlap};
//    if (isFeasible[i] and minEnergy > slackUnder2 - slackUnder1) {
//      maxOverflow[i] = overlap2 - overlap1;
//    } else {
//      maxOverflow[i] = overlap2 - overlap1 - (slackUnder2 - slackUnder1);
//    }
//  }
//
//  if (maxOverflow[i] < 0 /*or solver.num_cons_propagations == 70*/) {
//    std::cout << "bug (" << solver.num_cons_propagations << " est_i " << est(i)
//              << " ect_i " << ect(i) << ")\n";
//
//    std::cout << "max_Overflow_is_equal = " << maxOverflow[i]
//              << " contact time " << contact_time[i] << profile[contact[i]]
//              << " --- " << profile[est_[ip]] << ")\n";
//
//    for (auto j : lct_order) {
//      if (lct(j) > lct(k))
//        break;
//      // if (ect(i) >= contact_time[i])
//      std::cout << std::setw(3) << j << ": " << asciiArt(j) << std::endl;
//    }
//    std::cout << std::endl
//              << std::setw(3) << i << ": " << asciiArt(i) << std::endl;
//
//    std::cout << "capacity = " << capacity.max(solver) << " id = " << this->id()
//              << std::endl;
//
//    std::cout << profile << std::endl;
//
//    exit(1);
//  }
//}
//
//template <typename T>
//void CumulativeEdgeFinding<T>::computeBound(const int i) {
//  T E{0};
//  alpha = -1;
//  beta = -1;
//  T minSlack[2] = {Constant::Infinity<T>, Constant::Infinity<T>};
//  for (auto j : lct_order) {
//    if (lct(j) == lct(i))
//      break;
//    E += energy(j);
//    if (lct(j) <= ect(i) and est(i) < lct(j)) {
//      auto slack{(capacity.max(solver) - mindemand(i)) * lct(j) - E};
//      if (slack < minSlack[0] and profile[lct_[j]].overflow > 0) {
//        minSlack[0] = slack;
//        alpha = j;
//      }
//    }
//    // else
//    if (lct(j) > ect(i)) {
//      auto slack{capacity.max(solver) * lct(j) - E};
//      if (slack < minSlack[1] and profile[lct_[j]].overflow > 0) {
//        minSlack[1] = slack;
//        beta = j;
//      }
//    }
//  }
//}
//
//// template <typename T> T CumulativeEdgeFinding<T>::scheduleOmega(const T C,
//// const T max_lct) {
////
////
////#ifdef DBG_SEF
////   if (DBG_SEF) {
////     std::cout << "[schedule tasks until " << *lc_ptr
////               << "] profile=" << std::endl
////               << profile;
////   }
////   assert(verify());
////#endif
////
//////
////
////  auto saved_size{profile.size()};
//////  auto sentinel{profile.end()};
//////  --sentinel;
////
////  auto next{profile.begin()};
////  overflow = 0;
////  T omega_ect{-Constant::Infinity<T>};
////  T S{0};
////  T h_req{0};
////
////  while (next->time < max_lct) {
////    auto t{next};
////    ++next;
////
////#ifdef DBG_SEF
////    if (DBG_SEF) {
////      std::cout << "jump to e_" << t.index << " = t_" << t->time << ":";
////    }
////#endif
////
////    t->overflow = overflow;
////    auto l = next->time - t->time;
////
////    // S is the sum of demands of the tasks that could be processed at time
////    S += t->incrementMax;
////    // h_max is the min between the resource's capacity and the total of the
////    // demands
////    auto h_max{std::min(S, C)};
////    // h_req is the total demand counting tasks processed at their earliest
////    h_req += t->increment;
////    // h_cons is the amount of resource actually used in the optimistic
////    scenario
////    // (min between what is required + due from earlier, and what is
////    available) auto h_cons{std::min(h_req + overflow, h_max)};
////
////#ifdef DBG_SEF
////    if (DBG_SEF) {
////      std::cout << " h_max=" << h_max << ", h_req=" << h_req
////                << ", h_cons=" << h_cons << ", ov=" << overflow;
////      if (overflow > 0) {
////        std::cout << "->" << overflow - ((h_cons - h_req) * l)
////                  << " @t=" << next->time;
////      }
////    }
////#endif
////
////    // there is some overflow, and it will be resorbed by the next time point
////    if (overflow > 0 and overflow < ((h_cons - h_req) * l)) {
////      // then we create a new time point for the moment it will be resorbed
////      // l = std::max(Gap<T>::epsilon(), ceil_division<T>(overflow, (h_cons -
////      // h_req)));
////      l = std::max(Gap<T>::epsilon(), overflow / (h_cons - h_req));
////      auto new_event{
////          profile.create_element(t->time + l, t->capacity, 0, 0, 0, 0, 0, 0,
////          0)};
////      profile.add_after(t.index, new_event);
////      next = profile.at(new_event);
////    }
////    // overflow is the deficit on resource for that period (because tasks are
////    // set to their earliest)profile[est_[ip]].time = _est;
////    overflow += (h_req - h_cons) * l;
////    // once there
////    t->capacity = C - h_cons;
////    if (overflow > 0)
////      omega_ect = Constant::Infinity<T>;
////    else if (t->capacity < C)
////      omega_ect = profile[profile.next(t.index)].time;
////
////#ifdef DBG_SEF
////    if (DBG_SEF) {
////      if (omega_ect != -Constant::Infinity<T>)
////        std::cout << ", ect=" << omega_ect;
////      std::cout << std::endl;
////    }
////#endif
////  }
////
////  while (profile.size() > saved_size) {
////    profile.pop_back();
////  }
////
////  return omega_ect;
////}
////
////
////
//
//template <typename T>
//T CumulativeEdgeFinding<T>::scheduleOmegaDetection(const T C, const int i
////                                                   ,const T max_lct
//) {
//
//  //
//
//  auto saved_size{profile.size()};
//
//  auto stop{profile.end()};
//  for (auto p{profile.begin()}; p != stop; ++p) {
//    p->reset(capacity.max(solver));
//  }
//
//#ifdef DBG_SEF
//  if (DBG_SEF) {
//    std::cout << "[schedule tasks until " << *lc_ptr
//              << "] profile=" << std::endl
//              << profile;
//  }
//  assert(verify());
//#endif
//
//  auto next{profile.begin()};
//  overflow = 0;
//  T omega_ect{-Constant::Infinity<T>};
//  T S{0};
//  T h_req{0};
//
//  // while (next->time <= max_lct and next != ) {
//  while (next != stop) {
//    auto t{next};
//    ++next;
//
//#ifdef DBG_SEF
//    if (DBG_SEF) {
//      std::cout << "jump to e_" << t.index << " = t_" << t->time << ":";
//    }
//#endif
//
//    t->overflow = overflow;
//    auto l = (next == stop ? Constant::Infinity<T> : next->time - t->time);
//    //      if(next->time == Constant::Infinity<T>) {
//    //
//    //      }
//
//    // S is the sum of demands of the tasks that could be processed at time
//    S += t->incrementMax;
//    // h_max is the min between the resource's capacity and the total of the
//    // demands
//    auto h_max{std::min(S, C)};
//    // h_req is the total demand counting tasks processed at their earliest
//    h_req += t->increment;
//    // h_cons is the amount of resource actually used in the optimistic scenario
//    // (min between what is required + due from earlier, and what is available)
//    auto h_cons{std::min(h_req + overflow, h_max)};
//
//    //      auto absorbtion{(h_cons > h_req) and (l > 0)};
//
//    //      if(absorbtion)
//    if (overflow > 0) {
//      //          if((h_cons - h_req) == 0) {
//      //
//      //#ifdef DBG_SEF
//      //    if (DBG_SEF) {
//      //      std::cout << " h_max=" << h_max << ", h_req=" << h_req
//      //                << ", h_cons=" << h_cons
//      //                << ", ov=" << (overflow - (h_req - h_cons) * l);
//      //      if (overflow > 0) {
//      //        std::cout << "->" << overflow << " @t=" << next->time;
//      //      }
//      //
//      //      if (omega_ect != -Constant::Infinity<T>)
//      //        std::cout << ", ect=" << omega_ect;
//      //      std::cout << std::endl;
//      //    }
//      //#endif
//      //
//      //              break;
//      //          }
//
//      if ((h_cons - h_req) != 0) {
//
//        auto al{std::max(Gap<T>::epsilon(), overflow / (h_cons - h_req))};
//
//        //      std::cout << ", l=" << l << ", al="  ;
//
//        // there is some overflow, and it will be resorbed by the next time
//        // point
//        if (
//            //              overflow > 0 and
//            al < l) {
//          //        overflow < ((h_cons - h_req) * l)) {
//          // then we create a new time point for the moment it will be resorbed
//          // l = std::max(Gap<T>::epsilon(), ceil_division<T>(overflow, (h_cons
//          // - h_req)));
//          l = al;
//
//          //#ifdef DBG_SEF
//          //    if (DBG_SEF) {
//          //      std::cout << " overflow absorbtion event @t" << t->time + l <<
//          //      std::endl;
//          //    }
//          //#endif
//
//          auto new_event{profile.create_element(t->time + l, t->capacity, 0, 0,
//                                                0, 0, 0, 0, 0)};
//          profile.add_after(t.index, new_event);
//          next = profile.at(new_event);
//        }
//      }
//    }
//
//    //      std::cout << l;
//
//    // overflow is the deficit on resource for that period (because tasks are
//    // set to their earliest)profile[est_[ip]].time = _est;
//    overflow += (h_req - h_cons) * l;
//    // once there
//    t->capacity = C - h_cons;
//
//    if (overflow > 0)
//      omega_ect = Constant::Infinity<T>;
//    else if (t->capacity < C)
//      omega_ect = profile[profile.next(t.index)].time;
//
//#ifdef DBG_SEF
//    if (DBG_SEF) {
//      std::cout << " h_max=" << h_max << ", h_req=" << h_req
//                << ", h_cons=" << h_cons
//                << ", ov=" << (overflow - (h_req - h_cons) * l);
//      if (overflow > 0) {
//        std::cout << "->" << overflow << " @t=" << next->time;
//      }
//
//      if (omega_ect != -Constant::Infinity<T>)
//        std::cout << ", ect=" << omega_ect;
//      std::cout << std::endl;
//    }
//#endif
//  }
//
//  if (overflow > 0) {
//    //      contact[i] = -1;
//    //      auto t{profile.reverse_at(next.index)};
//    //      while(t != )
//    //
//    //
//    //      auto t{next};
//    //      --t;
//    //      while(t )
//
//    contact[i] = -1;
//    auto previous{next};
//    --previous;
//    while (previous != profile.begin()) {
//      auto t{previous};
//      --previous;
//      if (t->overflow < overflow) {
//
//        contact_time[i] = t->time;
//
//        while (t.index > 3 * static_cast<int>(task.size()))
//          --t;
//
//        contact[i] = t.index;
//
//        break;
//      }
//    }
//    if (contact[i] == -1 and previous == profile.begin()) {
//      contact[i] = previous.index;
//      contact_time[i] = previous->time;
//    }
//
////
//#ifdef DBG_SEF
//    if (DBG_SEF) {
//      std::cout << "contact[" << i << "] = " << profile[contact[i]] << " ("
//                << contact[i] << ") contact time = " << contact_time[i]
//                << std::endl;
//    }
//#endif
//  }
//
//  while (profile.size() > saved_size) {
//    profile.pop_back();
//  }
//
//  return omega_ect;
//}
//
//template <typename T>
//T CumulativeEdgeFinding<T>::scheduleOmegaAdjustment(const T C, const int i
////                                                    ,const T max_lct
//) {
//
//#ifdef DBG_SEF
//  if (DBG_SEF) {
//    std::cout << "[schedule tasks until " << *lc_ptr
//              << "] profile=" << std::endl
//              << profile;
//  }
//  assert(verify());
//#endif
//
//  //
//
//  auto saved_size{profile.size()};
//  //  auto sentinel{profile.end()};
//  //  --sentinel;
//
//  auto stop{profile.end()};
//  for (auto p{profile.begin()}; p != stop; ++p) {
//    p->reset(capacity.max(solver));
//  }
//
//  auto next{profile.begin()};
//  overflow = 0;
//  T omega_ect{-Constant::Infinity<T>};
//  T S{0};
//  T h_req{0};
//  T overlap{0};
//  T slackUnder{0};
//  T available{0};
//  T h{mindemand(i)};
//  isFeasible[i] = true;
//  tp_attributes_est_i.resize(3);
//  tp_attributes_ect_i.resize(3);
//
//  T prev_cons{0};
//
//  // while (next->time <= max_lct) {
//  while (next != stop) {
//    auto t{next};
//    ++next;
//
//#ifdef DBG_SEF
//    if (DBG_SEF) {
//      std::cout << "jump to e_" << t.index << " = t_" << t->time << ":";
//    }
//#endif
//
//    t->overflow = overflow;
//    auto l = (next == stop ? Constant::Infinity<T> : next->time - t->time);
//    //    auto l = next->time - t->time;
//    next->overlap = overlap;
//    next->slackUnder = slackUnder;
//    next->available = available;
//
//    // S is the sum of demands of the tasks that could be processed at time
//    S += t->incrementMax;
//    // h_max is the min between the resource's capacity and the total of the
//    // demands
//    auto h_max{std::min(S, C)};
//    // h_req is the total demand counting tasks processed at their earliest
//    h_req += t->increment;
//    // h_cons is the amount of resource actually used in the optimistic scenario
//    // (min between what is required + due from earlier, and what is available)
//    auto h_cons{std::min(h_req + overflow, h_max)};
//    t->consumption = h_cons;
//
//    //      std::cout << ", l=" << l << ", al="  ;
//
//    // there is some overflow, and it will be resorbed by the next time point
//    if (overflow > 0) {
//
//      if ((h_cons - h_req) != 0) {
//
//        auto al{std::max(Gap<T>::epsilon(), overflow / (h_cons - h_req))};
//
//        if (al < l) {
//          //        overflow < ((h_cons - h_req) * l)) {
//          // then we create a new time point for the moment it will be resorbed
//          // l = std::max(Gap<T>::epsilon(), ceil_division<T>(overflow, (h_cons
//          // - h_req)));
//          //      l = std::max(Gap<T>::epsilon(), overflow / (h_cons - h_req));
//
//          l = al;
//          auto new_event{profile.create_element(t->time + l, t->capacity, 0, 0,
//                                                0, 0, 0, 0, 0)};
//          profile.add_after(t.index, new_event);
//          next = profile.at(new_event);
//        }
//      }
//    }
//    // overflow is the deficit on resource for that period (because tasks are
//    // set to their earliest)profile[est_[ip]].time = _est;
//    overflow += (h_req - h_cons) * l;
//    // once there
//    t->capacity = C - h_cons;
//    // check the feasibility of the scheduling
//    if (overflow > 0)
//      isFeasible[i] = false;
//
//    //      if(l>0) {
//    // update the values of overlap, slack under and available
//    overlap += (std::max(h_cons - (C - h), 0) * l);
//    slackUnder += (std::max(std::min(C - h, h_max) - h_cons, 0) * l);
//    available += (std::min(C - h_cons, h) * l);
//    if ((t->time) < est(i)) {
//      if (next->time == est(i)) {
//        tp_attributes_est_i[0] = overlap;
//        tp_attributes_est_i[1] = slackUnder;
//        tp_attributes_est_i[2] = available;
//      } else if (next->time > est(i)) {
//        tp_attributes_est_i[0] =
//            t->overlap + (std::max(h_cons - (C - h), 0) * (est(i) - t->time));
//        tp_attributes_est_i[1] =
//            t->slackUnder +
//            (std::max(std::min(C - h, h_max) - h_cons, 0) * (est(i) - t->time));
//        tp_attributes_est_i[2] =
//            t->available + (std::min(C - h_cons, h) * (est(i) - t->time));
//      }
//    }
//
//    if ((t->time) < ect(i)) {
//      if (next->time == ect(i)) {
//        tp_attributes_ect_i[0] = overlap;
//        tp_attributes_ect_i[1] = slackUnder;
//        tp_attributes_ect_i[2] = available;
//      } else if (next->time > ect(i)) {
//        tp_attributes_ect_i[0] =
//            t->overlap + (std::max(h_cons - (C - h), 0) * (ect(i) - t->time));
//        tp_attributes_ect_i[1] =
//            t->slackUnder +
//            (std::max(std::min(C - h, h_max) - h_cons, 0) * (ect(i) - t->time));
//        tp_attributes_ect_i[2] =
//            t->available + (std::min(C - h_cons, h) * (ect(i) - t->time));
//      }
//    }
//
//    //      }
//
//    if (overflow > 0)
//      omega_ect = Constant::Infinity<T>;
//    else if (t->capacity < C)
//      omega_ect = profile[profile.next(t.index)].time;
//
//#ifdef DBG_SEF
//    if (DBG_SEF) {
//
//      std::cout << "l=" << l << " h_max=" << h_max << ", h_req=" << h_req
//                << ", h_cons=" << h_cons << "|" << prev_cons
//                << ", ov=" << (overflow - (h_req - h_cons) * l);
//      if (overflow > 0) {
//        std::cout << "->" << overflow << " @t=" << next->time;
//      }
//      std::cout << ", overlap=" << overlap << ", slackUnder=" << slackUnder
//                << ", available=" << available;
//
//      if (omega_ect != -Constant::Infinity<T>)
//        std::cout << ", ect=" << omega_ect;
//      std::cout << std::endl;
//    }
//#endif
//
//    if (l > 0) {
//      prev_cons = h_cons;
//    }
//  }
//  // std::cout << " is Feasible " << " is Feasible "<< isFeasible[i] << "\n";
//
//  while (profile.size() > saved_size) {
//    profile.pop_back();
//  }
//
//  return omega_ect;
//}



template <typename T>
void CumulativeEdgeFinding<T>::xplain(const Literal<T> l, const hint h,
                                      std::vector<Literal<T>> &Cl) {

//
//#ifdef DBG_EXPLCE
//  std::cout << "explain (" << explanation[h].size() << ") " << solver.pretty(l)
//            << ":\n";
//#endif
//
//  if (l == Solver<T>::Contradiction) {
//    //        std::cout << "xplain contradiction: TODO\n";
//    //        exit(1);
//
//    for (auto i : explanation[h].omega) {
//
//#ifdef DBG_EXPLCE
//      std::cout << " * "
//                << solver.pretty(task[i].start.after(explanation_lb[h]))
//                << " and "
//                << solver.pretty(task[i].end.before(explanation_ub[h]))
//                << std::endl;
//#endif
//
//      Cl.push_back(task[i].end.after(explanation[h].lb_omega));
//      Cl.push_back(task[i].end.before(explanation[h].ub_omega));
//    }
//
//  } else if (l.variable() == schedule.end.id()) {
//    std::cout << "xplain global bound: TODO\n";
//    exit(1);
//  } else {
//
//    for (auto i : explanation[h].omega) {
//
//#ifdef DBG_EXPLCE
//      std::cout << " * "
//                << solver.pretty(task[i].start.after(explanation_lb[h]))
//                << " and "
//                << solver.pretty(task[i].end.before(explanation_ub[h]))
//                << std::endl;
//#endif
//
//      Cl.push_back(task[i].start.after(explanation[h].lb_omega));
//      Cl.push_back(task[i].end.before(explanation[h].ub_omega));
//    }
//
//#ifdef DBG_EXPLCE
//    std::cout << " AND "
//              << solver.pretty(
//                     task[explanation[h].i].end.after(explanation[h].lb_i))
//              << std::endl;
//#endif
//
//    Cl.push_back(task[explanation[h].i].start.after(explanation[h].lb_i));
//  }
//
//  //      Cl.push_back(geq<T>(l.variable(), explanation_lb[h]));
}

template <typename T>
std::ostream &CumulativeEdgeFinding<T>::display(std::ostream &os) const {
  os << "Cumulative Edge-Finding"; // "data/sample/j309_5.sm"

#ifdef DBG_CEDGEFINDING
  os << "[" << this->id() << "]";
#endif

  os << "(";
  for (auto &t : task) {
    std::cout << " t" << t.id();
  }
  std::cout << " )";
  return os;
}

template <typename T>
std::ostream &CumulativeEdgeFinding<T>::print_reason(std::ostream &os, const hint) const {
  os << "cumulative-edge-finding";
  return os;
}

} // namespace tempo

#endif
