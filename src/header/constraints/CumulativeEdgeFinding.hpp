/************************************************
 * Tempo CumulativeEdgeFinding.hpp
 * Implementation of the "strong edge-finding" algorithm as described in
 * Vincent Gingras and Claude-Guy Quimper. 2016. Generalizing the edge-finder rule for the cumulative constraint. In Proceedings of the Twenty-Fifth International Joint Conference on Artificial Intelligence (IJCAI'16). AAAI Press, 3103–3109.
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
#include <numeric>

#include "Explanation.hpp"
#include "Global.hpp"
#include "constraints/Constraint.hpp"
#include "util/Profile.hpp"
#include "Model.hpp"


namespace tempo {

template<typename T>
class Solver;

template <typename T> class CumulativeEdgeFinding : public Constraint<T> {
private:
  Solver<T> &m_solver;
    NumericVar<T> capacity;
  Interval<T> schedule;
  std::vector<Interval<T>> the_tasks;
    std::vector<NumericVar<T>> demand;
  std::vector<std::vector<Literal<T>>> precedence;

  // helpers
    List<Timepoint<T>> profile;
    std::vector<int> est_;
    std::vector<int> ect_;
    std::vector<int> lct_;
    int sentinel;
    
    std::vector<std::vector<int>> triggers;

 public:
    
    static int est_flag;
    static int ect_flag;
    static int lct_flag;
    static int dem_flag;
    
    
  template <typename ItTask, typename ItNVar, typename ItBVar>
  CumulativeEdgeFinding(Solver<T> &solver, const Interval<T> sched, const NumericVar<T> capacity, const ItTask beg_task,
                         const ItTask end_task, const ItNVar beg_dem, const ItBVar beg_disj);
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

  bool notify(const Literal<T>, const int rank) override;
  void post(const int idx) override;
  void propagate() override;

  void xplain(const Literal<T> l, const hint h,
              std::vector<Literal<T>> &Cl) override;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

  std::string prettyTask(const int i) const;
};

template <typename T>
int CumulativeEdgeFinding<T>::est_flag = 0;

template <typename T>
int CumulativeEdgeFinding<T>::ect_flag = 1;

template <typename T>
int CumulativeEdgeFinding<T>::lct_flag = 2;

template <typename T>
int CumulativeEdgeFinding<T>::dem_flag = 3;

template <typename T>
std::string CumulativeEdgeFinding<T>::prettyTask(const int i) const {
  std::stringstream ss;
  ss << "t" << the_tasks[i].id() << ": [" << est(i) << ".." << lct(i) << "] (" << mindemand(i) << "x" << minduration(i) << ")";
  return ss.str();
}



template <typename T> T CumulativeEdgeFinding<T>::est(const unsigned i) const {
  return the_tasks[i].getEarliestStart(m_solver);
}

template <typename T> T CumulativeEdgeFinding<T>::lst(const unsigned i) const {
  return the_tasks[i].getLatestStart(m_solver);
}

template <typename T> T CumulativeEdgeFinding<T>::ect(const unsigned i) const {
  return the_tasks[i].getEarliestEnd(m_solver);
}

template <typename T> T CumulativeEdgeFinding<T>::lct(const unsigned i) const {
  return the_tasks[i].getLatestEnd(m_solver);
}

template <typename T>
T CumulativeEdgeFinding<T>::minduration(const unsigned i) const {
  return the_tasks[i].minDuration(m_solver);
}

template <typename T>
T CumulativeEdgeFinding<T>::maxduration(const unsigned i) const {
  return the_tasks[i].maxDuration(m_solver);
}

template <typename T>
T CumulativeEdgeFinding<T>::mindemand(const unsigned i) const {
  return demand[i].min(m_solver);
}

template <typename T>
T CumulativeEdgeFinding<T>::maxdemand(const unsigned i) const {
  return demand[i].max(m_solver);
}

template <typename T>
template <typename ItTask, typename ItNVar, typename ItBVar>
CumulativeEdgeFinding<T>::CumulativeEdgeFinding(Solver<T> &solver, const Interval<T> sched, const NumericVar<T> cap, const ItTask beg_task,
                       const ItTask end_task, const ItNVar beg_dem, const ItBVar beg_disj)
    : m_solver(solver) {
        
    schedule = sched,
    capacity = cap;
    
    Constraint<T>::priority = Priority::Low;
    
    auto dp{beg_dem};
    for (auto jp{beg_task}; jp != end_task; ++jp) {
        the_tasks.push_back(*jp);
        demand.push_back(*dp);
        ++dp;
    }
        
    precedence.resize(the_tasks.size());
        est_.resize(the_tasks.size());
        ect_.resize(the_tasks.size());
        lct_.resize(the_tasks.size());
        
    
        auto maxcap{capacity.max(m_solver)};
    for (unsigned i = 0; i < the_tasks.size(); ++i) {
        precedence[i].resize(the_tasks.size());
        est_[i] = profile.create_element(est(i), maxcap, mindemand(i), mindemand(i), 0, 0, 0, 0);
        ect_[i] = profile.create_element(ect(i), maxcap, -mindemand(i), 0, 0, 0, 0, 0);
        lct_[i] = profile.create_element(lct(i), maxcap, 0, -mindemand(i), 0, 0, 0, 0);
    }
    
    auto ep{beg_disj};
    for (auto ip{beg_task}; ip != end_task; ++ip) {
        for (auto jp{ip + 1}; jp != end_task; ++jp) {
            auto x{*ep};
            
            auto i{std::distance(beg_task, ip)};
            auto j{std::distance(beg_task, jp)};
            precedence[i][j] = m_solver.boolean.getLiteral(false, x);
            
            x = *(++ep);
            
            precedence[j][i] = m_solver.boolean.getLiteral(false, x);
            ++ep;
        }
    }

        
        std::vector<int> event_ordering(profile.size());
        std::iota(event_ordering.begin(), event_ordering.end(), 1);
        
        std::sort(event_ordering.begin(), event_ordering.end(), [this](const int i, const int j) {return this->profile[i].time < this->profile[j].time;});
                  
        auto elt{event_ordering.begin()};
        profile.add_front(*elt);
        
//        std::cout << profile << std::endl;
        
        while(++elt != event_ordering.end()) {
            profile.add_after(*(elt-1), *elt);
            
//            std::cout << "add " <<  *elt << " after " << *(elt-1) << std::endl;
            
//            std::cout << profile << std::endl;
        }
        
        sentinel = profile.create_element(profile[*event_ordering.rbegin()].time, 0, 0, 0, 0, 0, 0, 0);
        
        profile.add_after(*event_ordering.rbegin(), sentinel);
        
        std::cout << profile << std::endl;
                  
}

template <typename T> CumulativeEdgeFinding<T>::~CumulativeEdgeFinding() {}

template <typename T> void CumulativeEdgeFinding<T>::post(const int idx) {

  Constraint<T>::cons_id = idx;
  Constraint<T>::idempotent = true;

#ifdef DBG_EDGEFINDING
  if (DBG_EDGEFINDING) {
    std::cout << "post " << *this << std::endl;
  }
#endif

    
    

    for (size_t i{0}; i < the_tasks.size(); ++i) {
        auto k{m_solver.wake_me_on(lb<T>(the_tasks[i].start.id()), this->id())};
        if(triggers.size() <= k) triggers.resize(k+1);
        triggers[k].push_back(est_flag + 4*i);
    }
    for (size_t i{0}; i < the_tasks.size(); ++i) {
        auto k{m_solver.wake_me_on(lb<T>(the_tasks[i].end.id()), this->id())};
        if(triggers.size() <= k) triggers.resize(k+1);
        triggers[k].push_back(ect_flag + 4*i);
    }
    for (size_t i{0}; i < the_tasks.size(); ++i) {
        auto k{m_solver.wake_me_on(ub<T>(the_tasks[i].end.id()), this->id())};
        if(triggers.size() <= k) triggers.resize(k+1);
        triggers[k].push_back(lct_flag + 4*i);
    }
    for (size_t i{0}; i < the_tasks.size(); ++i) {
        auto k{m_solver.wake_me_on(lb<T>(demand[i].id()), this->id())};
        if(triggers.size() <= k) triggers.resize(k+1);
        triggers[k].push_back(dem_flag + 4*i);
    }
    
}

template <typename T>
bool CumulativeEdgeFinding<T>::notify(const Literal<T> l, const int r) {
//    int n{static_cast<int>(the_tasks.size())};
    
#ifdef DBG_SEF
        std::cout << "notify " << m_solver.pretty(l) << std::endl;
#endif
    
    for(auto t : triggers[r]) {
        auto flag{t%4};
        auto i{t/4};
        
        if(flag == est_flag) {
            std::cout << "est of task " << the_tasks[i] << std::endl;
            
            // est of the_tasks[r] has changed (increased)
            auto t = profile[est_[i]].time = est(i);
            auto s{profile.at(est_[i])};
            auto j{s};
            do
                ++j;
             while((*j).time < t);
            --j;
            if(j!=s) {
    #ifdef DBG_SEF
            std::cout << " -insert " << profile[est_[i]] << " after " << *j << std::endl;
    //            std::cout << profile << std::endl;
    #endif
                
                profile.remove(est_[i]);
                profile.add_after(j.index, est_[i]);
            }
    #ifdef DBG_SEF
            else std::cout << " - " << profile[est_[i]] << "'s rank has not changed" << std::endl;
    #endif
            
        } else if(flag == ect_flag) {
            std::cout << "ect of task " << the_tasks[i] << std::endl;
            
            // ect of the_tasks[r-n] has changed (increased)
            auto t = profile[ect_[i]].time = ect(i);
            auto s{profile.at(ect_[i])};
            auto j{s};
            do
                ++j;
             while((*j).time < t);
            --j;
            if(j!=s) {
    #ifdef DBG_SEF
            std::cout << " -insert " << profile[est_[i]] << " after " << *j << std::endl;
    //            std::cout << profile << std::endl;
    #endif
                
                profile.remove(ect_[i]);
                profile.add_after(j.index, ect_[i]);
            }
    #ifdef DBG_SEF
            else std::cout << " - " << profile[ect_[i]] << "'s rank has not changed" << std::endl;
    #endif
            
        } else if(flag == lct_flag) {
            std::cout << "lct of task " << the_tasks[i] << std::endl;
            
            // lct of the_tasks[r-2*n] has changed (decreased)
            auto t = profile[lct_[i]].time = lct(i);
            auto j{profile.at(lct_[i])};
            
//            std::cout << *j << "\nin profile:\n" << profile << std::endl;
            
            
            do
                --j;
             while((*j).time > t);
            if(profile.next(j.index) != lct_[i]) {
    #ifdef DBG_SEF
            std::cout << " -insert " << profile[lct_[i]] << " after " << *j << std::endl;
    //            std::cout << profile << std::endl;
    #endif
                
                profile.remove(lct_[i]);
                profile.add_after(j.index, lct_[i]);
            }
    #ifdef DBG_SEF
            else std::cout << " - " << profile[lct_[i]] << "'s rank has not changed" << std::endl;
    #endif
            
        } else {
            std::cout << "demand of task " << the_tasks[i] << std::endl;
        }
    }
    
    
//
//    
//    if(r < n) {
//        auto i{r};
//        
//        std::cout << profile[est_[i]].time << " < " << est(i) << std::endl;
//        if(profile[est_[i]].time >= est(i)) {
//            std::cout << profile << std::endl;
//        }
//        
//        assert(profile[est_[i]].time < est(i));
//        
//        // est of the_tasks[r] has changed (increased)
//        auto t = profile[est_[r]].time = est(i);
//        auto s{profile.at(est_[i])};
//        auto j{s};
//        do
//            ++j;
//         while((*j).time < t);
//        --j;
//        if(j!=s) {
//#ifdef DBG_SEF
//        std::cout << " -insert " << profile[est_[i]] << " after " << *j << std::endl;
////            std::cout << profile << std::endl;
//#endif
//            
//            profile.remove(est_[i]);
//            profile.add_after(j.index, est_[i]);
//        }
//#ifdef DBG_SEF
//        else std::cout << " - " << profile[est_[i]] << "'s rank has not changed" << std::endl;
//#endif
//        
//    } else if(r < 2*n) {
//        
//        auto i{r-n};
//        
//        // ect of the_tasks[r-n] has changed (increased)
//        auto t = profile[ect_[r]].time = ect(i);
//        auto s{profile.at(ect_[i])};
//        auto j{s};
//        do
//            ++j;
//         while((*j).time < t);
//        --j;
//        if(j!=s) {
//#ifdef DBG_SEF
//        std::cout << " -insert " << profile[est_[i]] << " after " << *j << std::endl;
////            std::cout << profile << std::endl;
//#endif
//            
//            profile.remove(ect_[i]);
//            profile.add_after(j.index, ect_[i]);
//        }
//#ifdef DBG_SEF
//        else std::cout << " - " << profile[ect_[i]] << "'s rank has not changed" << std::endl;
//#endif
//        
//    } else if(r < 3*n) {
//        
//        auto i{r%n};
//        
//        // lct of the_tasks[r-2*n] has changed (decreased)
//        auto t = profile[lct_[r]].time = lct(i);
//        auto j{profile.at(lct_[i])};
//        do
//            --j;
//         while((*j).time > t);
//        if(profile.next(j.index) != lct_[i]) {
//#ifdef DBG_SEF
//        std::cout << " -insert " << profile[lct_[i]] << " after " << *j << std::endl;
////            std::cout << profile << std::endl;
//#endif
//            
//            profile.remove(lct_[i]);
//            profile.add_after(j.index, lct_[i]);
//        }
//#ifdef DBG_SEF
//        else std::cout << " - " << profile[lct_[i]] << "'s rank has not changed" << std::endl;
//#endif
//        
//    } else if(r < 4*n) {
//        // demand of the_tasks[r-3*n] has changed (increased)
//        
//    }
    
  return true;
}


template <typename T> void CumulativeEdgeFinding<T>::propagate() {
        
    T previous{-Constant::Infinity<T>};
    for(auto e{profile.begin()}; e!=profile.end(); ++e) {
        if((*e).time < previous) {
            std::cout << "not ordered (event " << *e << ") in\n" << profile << std::endl;
            exit(1);
        }
        previous = (*e).time;
    }
    
    std::cout << "ok\n";
}



template <typename T>
void CumulativeEdgeFinding<T>::xplain(const Literal<T>, const hint,
                                       std::vector<Literal<T>> &) {

 
}

template <typename T>
std::ostream &CumulativeEdgeFinding<T>::display(std::ostream &os) const {
  os << "Cumulative Edge-Finding";

#ifdef DBG_EDGEFINDING
  os << "[" << this->id() << "]";
#endif

  os << "(";
  for (auto &t : the_tasks) {
    std::cout << " t" << t.id();
  }
  std::cout << " )";
  return os;
}

template <typename T>
std::ostream &CumulativeEdgeFinding<T>::print_reason(std::ostream &os,
                                                      const hint) const {
  os << "cumulative-edge-finding";
  return os;
}

} // namespace tempo

#endif
