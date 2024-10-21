/************************************************
 * Tempo CumulativeTimetabling.hpp
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

#ifndef TEMPO_CUMULATIVETIMETABLING_HPP
#define TEMPO_CUMULATIVETIMETABLING_HPP

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


#define OLD_EXPL


namespace tempo {


template <typename T = int>
struct ProfileDelta {
    T time{-Constant::Infinity<T>};
    T delta{0};
    
    
    std::ostream &display(std::ostream &os) const {
        os << std::right << std::setw(3) << delta << " @" << std::left << std::setw(3) << time;
        return os;
    }
};


template <typename T>
std::ostream &operator<<(std::ostream &os, const ProfileDelta<T> &x) {
    return x.display(os);
}



template<typename T>
class Solver;

template <typename T> class CumulativeTimetabling : public Constraint<T> {
private:
    Solver<T> &solver;
    NumericVar<T> capacity;
    
    std::vector<Interval<T>> task;
    std::vector<NumericVar<T>> demand;
    bool sign{true};
    
    typedef std::vector<ProfileDelta<T>> Profile;
    
    std::vector<int> mandatory_parts;
    Profile profile_unique_time_;
    Profile profile_non_unique_time_;
    std::vector<int> est_order;
    
    
    std::vector<Literal<T>> pruning;
    std::vector<std::vector<Literal<T>>> explanation;
    std::vector<std::vector<Literal<T>>> alt_explanation;
    Reversible<size_t> num_explanations;
    
    
public:
    
    template <typename ItTask, typename ItNVar>
    CumulativeTimetabling(Solver<T> &solver,
                          const NumericVar<T> capacity, const ItTask beg_task,
                          const ItTask end_task, const ItNVar beg_dem);
    virtual ~CumulativeTimetabling();
    
    // helpers
    T est(const unsigned i) const;
    T lst(const unsigned i) const;
    T ect(const unsigned i) const;
    T lct(const unsigned i) const;
    T minduration(const unsigned i) const;
    T maxduration(const unsigned i) const;
    T mindemand(const unsigned i) const;
    T maxdemand(const unsigned i) const;
    T energy(const unsigned i) const;
    
    bool notify(const Literal<T>, const int rank) override;
    void post(const int idx) override;
    void propagate() override;
    
    void buildProfile();
    void computeAdjustments();
    void pushTask(const int i, int profile_index, const T usage);
    void doPruning();
    
    hint newExplanation();
    void computeExplanation(const hint h, const int i, const T l, const T u);
//    void computeLBExplanation(const hint h, const int i, const T l, const T u);
//    void computeUBExplanation(const hint h, const int i, const T l, const T u);
    
    // function used in explanation
    void xplain(const Literal<T> l, const hint h,
                std::vector<Literal<T>> &Cl) override;
    
    std::ostream &display(std::ostream &os) const override;
    
    std::ostream &print_reason(std::ostream &os, const hint h) const override;
    
    std::string prettyTask(const int i) const;
    std::string asciiArt(const int i) const;
    void printProfile();
};


template <typename T>
std::string CumulativeTimetabling<T>::prettyTask(const int i) const {
    std::stringstream ss;
    ss << "t" << task[i].id() << ": [" << (sign == bound::lower ? est(i) : -lct(i))  << ".." << (sign == bound::lower ? lct(i) : -est(i)) << "] (p="
    << minduration(i) << ", c=" << mindemand(i) << ")";
    return ss.str();
}

template <typename T>
std::string CumulativeTimetabling<T>::asciiArt(const int i) const {
    std::stringstream ss;
    ss << std::setw(3) << std::right << mindemand(i) << "x" << std::setw(3)
    << std::left << minduration(i) << " " << std::right ;
//    std::cout.flush();
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

template <typename T> T CumulativeTimetabling<T>::est(const unsigned i) const {
    return (sign==bound::lower ? task[i].getEarliestStart(solver) : -task[i].getLatestEnd(solver));
}

template <typename T> T CumulativeTimetabling<T>::lst(const unsigned i) const {
    return (sign==bound::lower ? task[i].getLatestStart(solver) : -task[i].getEarliestEnd(solver));
}

template <typename T> T CumulativeTimetabling<T>::ect(const unsigned i) const {
    return (sign==bound::lower ? task[i].getEarliestEnd(solver) : -task[i].getLatestStart(solver));
}

template <typename T> T CumulativeTimetabling<T>::lct(const unsigned i) const {
    return (sign==bound::lower ? task[i].getLatestEnd(solver) : -task[i].getEarliestStart(solver));
}

template <typename T>
T CumulativeTimetabling<T>::minduration(const unsigned i) const {
    return task[i].minDuration(solver);
}

template <typename T>
T CumulativeTimetabling<T>::maxduration(const unsigned i) const {
    return task[i].maxDuration(solver);
}

template <typename T>
T CumulativeTimetabling<T>::mindemand(const unsigned i) const {
    return demand[i].min(solver);
}

template <typename T>
T CumulativeTimetabling<T>::maxdemand(const unsigned i) const {
    return demand[i].max(solver);
}

template <typename T>
T CumulativeTimetabling<T>::energy(const unsigned i) const {
    return mindemand(i) * minduration(i);
}

//template <typename T>
//T CumulativeTimetabling<T>::getFixedPart(const unsigned i) const {
//    return std::max(0, ect(i) - lst(i));
//}


template <typename T>
template <typename ItTask, typename ItNVar>
CumulativeTimetabling<T>::CumulativeTimetabling(
                                                Solver<T> &solver, const NumericVar<T> cap,
                                                const ItTask beg_task, const ItTask end_task, const ItNVar beg_dem)
: solver(solver), num_explanations(0, &(solver.getEnv())) {
    capacity = cap;
    
    
    
    
    Constraint<T>::priority = Priority::Medium;
    
    auto dp{beg_dem};
    for (auto jp{beg_task}; jp != end_task; ++jp) {
        task.push_back(*jp);
        demand.push_back(*dp);
        ++dp;
    }
    //    est_buffer.resize(task.size());
    //    lst_buffer.resize(task.size());
    //    ect_buffer.resize(task.size());
    //    lct_buffer.resize(task.size());
    
    est_order.resize(task.size());
    std::iota(est_order.begin(), est_order.end(), 0);
}

template <typename T> CumulativeTimetabling<T>::~CumulativeTimetabling() {}


template <typename T> void CumulativeTimetabling<T>::post(const int idx) {
    
    Constraint<T>::cons_id = idx;
    Constraint<T>::idempotent = true;
    
#ifdef DBG_CEDGEFINDING
    if (DBG_CEDGEFINDING) {
        std::cout << "post " << *this << std::endl;
    }
#endif
    
    
//    std::cout << "post TT" << this->id() << ", demands =\n";
    
    
    for (size_t i{0}; i < task.size(); ++i) {
        solver.wake_me_on(lb<T>(task[i].start.id()), this->id());
        solver.wake_me_on(ub<T>(task[i].end.id()), this->id());
        solver.wake_me_on(lb<T>(demand[i].id()), this->id());
        
//        std::cout << demand[i].min(solver) << ".." << demand[i].max(solver) << std::endl;
    }
}

template <typename T>
bool CumulativeTimetabling<T>::notify(const Literal<T>, const int) {
    return true;
}

template <typename T> void CumulativeTimetabling<T>::buildProfile() {
    
    mandatory_parts.clear();
    profile_non_unique_time_.clear();
    auto n{task.size()};
    for(unsigned i{0}; i<n; ++i) {
        if(task[i].mustExist(solver) and lst(i) < ect(i) and mindemand(i) > 0) {
            profile_non_unique_time_.emplace_back(lst(i), +mindemand(i));
            profile_non_unique_time_.emplace_back(ect(i), -mindemand(i));
            mandatory_parts.push_back(i);
            
            #ifdef DBG_TT
                if (DBG_TT) {
                    std::cout << "mandatory part on t" << task[i].id() << ": [" << lst(i) << ".." << ect(i) << "]\n";
                }
            #endif
        }
    }
    
    // Sort
    std::sort(profile_non_unique_time_.begin(), profile_non_unique_time_.end(),
              [](const auto& a, const auto& b) {return a.time < b.time;} );
    
    // Build profile with unique times
    profile_unique_time_.clear();
    profile_unique_time_.emplace_back(-Constant::Infinity<T>, 0);
    
    T usage{0};
    for (const ProfileDelta<T>& step : profile_non_unique_time_) {
        if (step.time == profile_unique_time_.back().time) {
            profile_unique_time_.back().delta += step.delta;
        } else {
            profile_unique_time_.push_back(step);
        }
        // Update usage.
        usage += step.delta;
    }
    // Check final usage to be 0.
    assert(usage == 0);
    
    // Scan to find max usage.
    T max_usage = 0;
    T l{-Constant::Infinity<T>};
    T u{Constant::Infinity<T>};
    T prev_time{l};
    T max_cap{capacity.max(solver)};
    for (const ProfileDelta<T>& step : profile_unique_time_) {
//        if(usage > max_cap) {
//            l = prev_time;
//            u = step.time;
//        }
        
        usage += step.delta;
        if (usage > max_usage) {
            max_usage = usage;
//            l = step.time;
        }
        if(usage > max_cap) {
            l = step.time;
            u = step.time+1;
        }
        
        
//        else if(usage < max_usage and prev_time == l) {
//            u = step.time;
//        }
        prev_time = step.time;
        
        
        

#ifdef DBG_TT
    if (DBG_TT) {
        std::cout << " - profile step " << std::setw(3) << step << " (U=" << usage << ", max=" << max_usage << ")\n";
    }
#endif

    }

    
    
    assert(usage == 0);
 
//    if(capacity.min(solver) < max_usage) {
//        
//        std::cout << "capacity lower bound (" << capacity.min(solver) << ") pruning (max usage=" << max_usage << " in [" << l << ".." << u << "]): " << solver.pretty(capacity.greaterThanOrEqual(max_usage)) << " #cp=" << solver.num_choicepoints << std::endl;
//        
//        printProfile();
//        
//        solver.set(capacity.greaterThanOrEqual(max_usage), {this, computeExplanation(-1,l,u)});
//    }
  
 
    if(max_cap < max_usage) {

//        hint h;
//        if(sign == bound::upper) {
//            h=newExplanation();
//            sign = bound::lower;
//            computeExplanation(h,-1,-u,-l);
//            sign = bound::upper;
//        } else {
//            h=newExplanation();
//            computeExplanation(h,-1,l,u);
//        }
        hint h{newExplanation()};
        computeExplanation(h,-1,l,u);
        
        throw Failure<T>({this, h});
        
    }
    
    // Add a sentinel.
    profile_unique_time_.emplace_back(Constant::Infinity<T>, 0);
}

// Update the start min for all tasks. Runs in O(n^2) and Omega(n).
template <typename T> void CumulativeTimetabling<T>::computeAdjustments() {
    std::sort(est_order.begin(), est_order.end(),
              [&](const int a, const int b) {return est(a) < est(b);});
    T usage = 0;
    int profile_index = 0;
    for (auto i : est_order) {
        auto est_i{est(i)};
        
        if(est_i == lst(i) and lct(i) == ect(i))
            continue;
        
//        std::cout << "est(i)=" << est_i << " >? " << profile_unique_time_[profile_index].time << std::endl;
        
        while (est_i > profile_unique_time_[profile_index].time) {
            assert(static_cast<size_t>(profile_index) < profile_unique_time_.size());
            ++profile_index;
            usage += profile_unique_time_[profile_index].delta;
            
//            std::cout << "YES est(i)=" << est_i << " >? " << profile_unique_time_[profile_index].time << std::endl;
        }
        pushTask(i, profile_index, usage);
    }
}

// Push the given task to new_start_min, defined as the smallest integer such
// that the profile usage for all tasks, excluding the current one, does not
// exceed capacity_ - task->demand on the interval
// [new_start_min, new_start_min + task->interval->DurationMin() ).
template <typename T> void CumulativeTimetabling<T>::pushTask(const int i, int profile_index, const T u) {
    
#ifdef DBG_TT
    if (DBG_TT) {
//        std::cout << " push task t" << std::left << std::setw(3) << task[i].id() << ": " << profile_unique_time_[profile_index] << " (u=" << u << ")\n";
        std::cout << " push " << prettyTask(i) << ": " << profile_unique_time_[profile_index] << " (u=" << u << ")\n";
    }
#endif
    
    
    // Init
    if (mindemand(i) == 0) {  // Demand can be null, nothing to propagate.
        return;
    }
    
    T usage{u};
    
    hint h{Constant::NoHint};
//    bool beenthere{false};
    
    const T residual_capacity = (capacity.max(solver) - mindemand(i));
    const T duration = minduration(i);
    
    const ProfileDelta<T>& first_prof_delta = profile_unique_time_[profile_index];
    
    T new_start_min = est(i);
    
    assert(first_prof_delta.time >= est(i));
    // The check above is with a '>='. Let's first treat the '>' case
    if (first_prof_delta.time > new_start_min) {
        // There was no profile delta at a time between interval->StartMin()
        // (included) and the current one.
        // As we don't delete delta's of 0 value, this means the current task
        // does not contribute to the usage before:
        
#ifdef DBG_TT
        if (DBG_TT) {
            std::cout << "lst(i)=" << lst(i) << ", ect(i)=" << ect(i) << " step.time=" << first_prof_delta.time << std::endl;
        }
#endif
        
        assert((lst(i) >= first_prof_delta.time) or
               (lst(i) >= ect(i)));
        // The 'usage' given in argument is valid at first_prof_delta.time. To
        // compute the usage at the start min, we need to remove the last delta.
        const T usage_at_start_min = usage - first_prof_delta.delta;
        if (usage_at_start_min > residual_capacity) {
            h = newExplanation();
            
#ifndef OLD_EXPL
            auto prev_start_min{new_start_min};
#endif
            
            new_start_min = profile_unique_time_[profile_index].time;
            
#ifndef OLD_EXPL
            computeExplanation(h, i, prev_start_min, new_start_min);
            if(sign == bound::lower) {
                pruning.push_back(task[i].start.after(new_start_min));
                explanation[h].push_back(task[i].start.after(prev_start_min));
            } else {
                pruning.push_back(task[i].end.before(-new_start_min));
                explanation[h].push_back(task[i].end.before(-prev_start_min));
            }
#endif
            
            
////            if(ect(i) < )
//            computeExplanation(h, i, est(i), new_start_min);
        }
        
#ifdef DBG_TT
        if (DBG_TT) {
            std::cout << " usage @est(" << task[i].id() << ")=" << est(i)
            << " is " << usage_at_start_min << ", residual capacity is "
            << residual_capacity << "\n";
        }
#endif
        
    }
 
    
    // Influence of current task
    const T start_max = lst(i);
    const T end_min = ect(i);
    const T demand_min = (start_max < end_min ? mindemand(i) : 0);
    
    
    
//    std::cout << h << std::endl;
//    exit(1);
    
    while (profile_unique_time_[profile_index].time <
           (duration + new_start_min)) {
        
#ifdef DBG_TT
    if (DBG_TT) {
        std::cout << " t" << task[i].id() << " cannot fit @" << new_start_min << " because it would finish at "
        << (duration + new_start_min) << ", i.e., later than next chunck (" << profile_unique_time_[profile_index].time << ") => check if can be put in parallel\n";
    }
#endif
        
        const ProfileDelta<T>& profile_delta = profile_unique_time_[profile_index];
        assert(static_cast<size_t>(profile_index) < profile_unique_time_.size());
        // Compensate for current task
        if (profile_delta.time == start_max) {
            usage -= demand_min;
        }
        if (profile_delta.time == end_min) {
            usage += demand_min;
        }
        // Increment time
        ++profile_index;
        assert(static_cast<size_t>(profile_index) < profile_unique_time_.size());
        // Does it fit?
        
        
        
        if (usage > residual_capacity) {
            
#ifdef DBG_TT
    if (DBG_TT) {
        std::cout << " t" << task[i].id() << " cannot be put in parallel: usage=" << usage << ", residual capacity="
         << residual_capacity << "\n";
    }
#endif
            
            if(h == Constant::NoHint) {
                h = newExplanation();
            }
            
//            beenthere = true;
            
//            computeExplanation(h, i, ect(i)-Gap<T>::epsilon(), new_start_min);
            
            
//            computeExplanation(h, i, std::max(new_start_min,ect(i))-Gap<T>::epsilon(), new_start_min);
            
#ifndef OLD_EXPL
            h = newExplanation();
            auto prev_start_min{new_start_min};
#endif
            
            new_start_min = profile_unique_time_[profile_index].time;
            
#ifndef OLD_EXPL
            computeExplanation(h, i, prev_start_min, new_start_min);
            if(sign == bound::lower) {
                pruning.push_back(task[i].start.after(new_start_min));
                explanation[h].push_back(task[i].start.after(prev_start_min));
            } else {
                pruning.push_back(task[i].end.before(-new_start_min));
                explanation[h].push_back(task[i].end.before(-prev_start_min));
            }
#endif
            
            
            
            
            assert(new_start_min > est(i));
            
        }
        usage += profile_unique_time_[profile_index].delta;
    }
    
    
#ifdef DBG_TT
    if (DBG_TT) {
        std::cout << " t" << task[i].id() << " can fit @" << new_start_min << "\n";
    }
#endif
    
//    if(beenthere != (est(i) < new_start_min)) {
//        std::cout << "wtf?\n";
//        exit(1);
//    }
    
#ifdef OLD_EXPL
    if(est(i) < new_start_min) {
        
//        std::cout << h << " / " << explanation.size() << std::endl;
        
//        h = newExplanation();
        
//        std::cout << est(i) << " -> " << (std::max(new_start_min,ect(i))-Gap<T>::epsilon()) << " to " << new_start_min << std::endl;
        
//        computeExplanation(h, i, std::max(new_start_min,ect(i))-Gap<T>::epsilon(), new_start_min);
        computeExplanation(h, i, est(i), new_start_min);
        
        if(sign == bound::lower) {
            pruning.push_back(task[i].start.after(new_start_min));
//            computeLBExplanation(h, i, est(i), new_start_min);
//            computeExplanation(h, i, est(i), new_start_min);
            explanation[h].push_back(task[i].start.after(est(i)));
        } else {
            pruning.push_back(task[i].end.before(-new_start_min));
//            computeExplanation(h, i, est(i), new_start_min);
            explanation[h].push_back(task[i].end.before(-est(i)));
//            computeUBExplanation(h, i, est(i), new_start_min);
        }
        
        
    }
#endif
    
    
 
//        computeExplanation(h, i, std::max(new_start_min,ect(i))-Gap<T>::epsilon(), new_start_min);
//        if(sign == bound::lower) {
//            pruning.push_back(task[i].start.after(new_start_min));
//            explanation[h].push_back(task[i].start.after(est(i)));
//        } else {
//            pruning.push_back(task[i].end.before(-new_start_min));
//            explanation[h].push_back(task[i].end.before(-est(i)));
//        }
        
        
//    }

}


//template <typename T> void CumulativeTimetabling<T>::computeLBExplanation(const hint h, const int x, const T l, const T u) {
//    
////    bool flag{task[x].start.id() == 26};
//    
//#ifdef DBG_EXPLCTT
//    if(DBG_EXPLCTT) {
//        std::cout << "compute expl for " << task[x].start.after(u)
//        << " demand(" << task[x].id() << ") = " << mindemand(x) << ", dur(" << task[x].id() << ")= " << minduration(x) << ")\n";
//    }
//#endif
//    
//    computeExplanation(h,x,l,u);
////    explanation[h].push_back(task[x].start.after(est(x)));
//        
//#ifdef DBG_EXPLCTT
//        if(DBG_EXPLCTT) {
//            std::cout << " + " << explanation[h].back() << std::endl;
//            for(auto l : explanation[h]) {
//                std::cout << " * " << l << std::endl;
//            }
//            std::cout << " -> " << pruning.back() << std::endl;
//        }
//#endif
//        
//    
//}
//
//
//template <typename T> void CumulativeTimetabling<T>::computeUBExplanation(const hint h, const int x, const T l, const T u) {
//    
////    bool flag{task[x].start.id() == 26};
//    
//#ifdef DBG_EXPLCTT
//    if(DBG_EXPLCTT) {
//        std::cout << "compute expl for " << task[x].end.before(-u)
//        << " demand(" << task[x].id() << ") = " << mindemand(x) << ", dur(" << task[x].id() << ")= " << minduration(x) << ")\n";
//    }
//#endif
//    
//    sign = bound::lower;
//    computeExplanation(h,x,-u,-l);
//    explanation[h].push_back(task[x].end.before(lct(x)));
//    sign = bound::upper;
//    
//
//#ifdef DBG_EXPLCTT
//        if(DBG_EXPLCTT) {
//            std::cout << " + " << explanation[h].back() << std::endl;
//            for(auto l : explanation[h]) {
//                std::cout << " * " << l << std::endl;
//            }
//            std::cout << " -> " << pruning.back() << std::endl;
//        }
//#endif
//    
//}
    
template <typename T> hint CumulativeTimetabling<T>::newExplanation() {
    auto e_idx{num_explanations};
    if (explanation.size() <= e_idx) {
        explanation.resize(e_idx + 1);
    } else {
        explanation[e_idx].clear();
    }
//    hint h{static_cast<hint>(num_explanations)};
    ++num_explanations;
    return static_cast<hint>(e_idx);
}


template <typename T> void CumulativeTimetabling<T>::computeExplanation(const hint h, const int x, const T l, const T u) {
    
    if(l >= u) {
        std::cout << "bug\n";
        exit(1);
    }
    
    
//    bool flag{x == -1 ? false : task[x].start.id() == 26};
    
    bool swap{false};
    auto lb{l};
    auto ub{u};
    if(sign == bound::upper) {
        swap = true;
        sign = bound::lower;
        lb = -u;
        ub = -l;
    }
        
//    computeExplanation(h,-1,-u,-l);
        
 
    
#ifdef DBG_EXPLCTT
    if(DBG_EXPLCTT) {
        std::cout << " b/c profile overload in [" << l << ".." << u << "] capacity = " << capacity.max(solver) << std::endl;
//        std::cout << "(" << pruning.back() << ")\n";
    }
#endif
    
    
    for(auto i : mandatory_parts) if(i != x) {
        
#ifdef DBG_EXPLCTT
        if(DBG_EXPLCTT) {
            std::cout << prettyTask(i) << ": ";
        }
#endif
        
        if(lst(i) <= ub and ect(i) > lb) {
            explanation[h].push_back(task[i].start.before(std::max(lst(i),lb)));
            explanation[h].push_back(task[i].end.after(std::min(ect(i),ub)));
            
#ifdef DBG_EXPLCTT
            if(DBG_EXPLCTT) {
                std::cout << task[i].start.before(std::max(lst(i),lb)) << " and " << task[i].end.after(std::min(ect(i),ub)) << std::endl;
            }
#endif
            
        }
#ifdef DBG_EXPLCTT
        else if(DBG_EXPLCTT) {
            std::cout << "does not contributes...\n";
        }
#endif
    }
    
    
    if(swap) {
        sign = bound::upper;
        lb = -u;
        ub = -l;
    }
    
}

//


template <typename T> void CumulativeTimetabling<T>::doPruning() {
    
#ifdef DBG_TT
            if (DBG_TT) {
                if(not pruning.empty())
                    std::cout << "apply pruning" << std::endl;
            }
#endif
    
//    auto h{static_cast<hint>(num_explanations)};
//    while(not pruning.empty()) {
//        
//#ifdef DBG_TT
//            if (DBG_TT) {
//                std::cout << "pruning (" << this->id() << "/" << solver.num_cons_propagations
//                << "): " << pruning.back() << std::endl;
//            }
//#endif
//        
//        solver.set(pruning.back(), {this, --h});
//        pruning.pop_back();
//    }
    auto h{static_cast<hint>(num_explanations-pruning.size())};
    for(auto p : pruning) {
        #ifdef DBG_TT
                    if (DBG_TT) {
                        std::cout << "pruning (" << this->id() << "/" << solver.num_cons_propagations
                        << "): " << p << std::endl;
                    }
        #endif
        
                solver.set(p, {this, h++});
    }
    pruning.clear();
//    while(not pruning.empty()) {
//        
//#ifdef DBG_TT
//            if (DBG_TT) {
//                std::cout << "pruning (" << this->id() << "/" << solver.num_cons_propagations
//                << "): " << pruning.back() << std::endl;
//            }
//#endif
//        
//        solver.set(pruning.back(), {this, --h});
//        pruning.pop_back();
//    }
}


template <typename T> void CumulativeTimetabling<T>::propagate() {
    sign = bound::lower;
    pruning.clear();
    
//    assert(mindemand(0) > -Constant::Infinity<T>);
    
//    if(solver.num_cons_propagations >= 403)
//        exit(1);
    
#ifdef DBG_TT
    if (DBG_TT) {
        std::cout << "\nstart TT" << this->id() << "'s " << solver.num_cons_propagations
        << "-th propagation at lvl " << solver.level() << " (C=" << capacity.max(solver) << ")\n";
        std::sort(est_order.begin(), est_order.end(),
                  [&](const int a, const int b) {return est(a) < est(b);});
        
        for (auto j : est_order) {
//            std::cout << prettyTask(j) << std::endl;
            std::cout << "task " << std::setw(3) << task[j].id() << ": " << asciiArt(j) << std::endl;
        }
    }
#endif
    
//    if(mindemand(0) == -Constant::Infinity<T>) {
//        
//        for(auto d : demand) {
//            std::cout << d.id() << std::endl;
//        }
//        
////        std::cout << solver << std::endl;
//        
//        exit(1);
//    }
    
    
    do {
  
//        std::cout << (sign == bound::lower ? "forward\n" : "backward\n");
        
#ifdef DBG_TT
    if (DBG_TT) {
        std::cout << (sign == bound::lower ? "forward\n" : "backward\n");
    }
#endif
        
        
        buildProfile();
        
//        assert(mindemand(0) > -Constant::Infinity<T>);
        
#ifdef DBG_TT
    if (DBG_TT and sign == bound::lower) {
        std::cout << "profile=\n";
        for(auto p : profile_unique_time_) {
            std::cout << p << std::endl;
        }
        printProfile();
    }
#endif
        
        
                
        computeAdjustments();
        
//        assert(mindemand(0) > -Constant::Infinity<T>);
        
        doPruning();
        
        sign ^= 1;
    } while(sign!=bound::lower);
    
//    for(sign = bound::lower; sign!=bound::upper; sign ^= 1) {
//        
////        assert(mindemand(0) > -Constant::Infinity<T>);
//        
//        std::cout << (sign == bound::lower ? "forward\n" : "backward\n");
//        
//#ifdef DBG_TT
//    if (DBG_TT) {
//        std::cout << (sign == bound::lower ? "forward\n" : "backward\n");
//    }
//#endif
//        
//        
//        buildProfile();
//        
////        assert(mindemand(0) > -Constant::Infinity<T>);
//        
//#ifdef DBG_TT
//    if (DBG_TT) {
//        std::cout << "profile=\n";
//        for(auto p : profile_unique_time_) {
//            std::cout << p << std::endl;
//        }
//        printProfile();
//    }
//#endif
//        
//        
//                
//        computeAdjustments();
//        
////        assert(mindemand(0) > -Constant::Infinity<T>);
//        
//        doPruning();
//        
////        assert(mindemand(0) > -Constant::Infinity<T>);
//    }
    
//    
//    if(solver.num_cons_propagations == 513)
//        exit(1);
}


template <typename T>
void CumulativeTimetabling<T>::xplain(const Literal<T>, const hint h,
                                      std::vector<Literal<T>> &Cl) {
    
    for(auto p : explanation[h]) {
        Cl.push_back(p);
    }
}

 template <typename T>
void CumulativeTimetabling<T>::printProfile() {
    T time{0};
    T usage{0};
    int W{4};
    for(auto step : profile_unique_time_) {
        if(step.time < time)
            continue;
        if(step.time == Constant::Infinity<T>) {
            for(auto i{0}; i<5; ++i) {
                std::cout << std::setw(W) << "|";
                for(auto u{0}; u<usage; ++u) {
                    std::cout << " ";
                }
                std::cout << "|\n";
            }
        } else if(step.delta > 0) {
            while(time++ < step.time) {
                std::cout << std::right << std::setw(W) << "|";
                for(auto u{0}; u<usage; ++u) {
                    std::cout << " ";
                }
                std::cout << "|\n";
            }
            std::cout << std::left << std::setw(W-1) << step.time << std::right << "|";
            for(auto u{0}; u<usage; ++u) {
                std::cout << " ";
            }
            std::cout << "+";
            for(auto u{0}; u<step.delta-1; ++u) {
                std::cout << "-";
            }
            usage += step.delta;
            std::cout << "+ " << usage << "\n";
        } else if(step.delta < 0) {
            while(time++ < step.time) {
                std::cout << std::right << std::setw(W) << "|";
                for(auto u{0}; u<usage; ++u) {
                    std::cout << " ";
                }
                std::cout << "|\n";
            }
            usage += step.delta;
            std::cout << std::left << std::setw(W-1) << step.time << std::right << "|";
            for(auto u{0}; u<usage; ++u) {
                std::cout << " ";
            }
            std::cout << "+";
            for(auto u{0}; u<-step.delta-1; ++u) {
                std::cout << "-";
            }
            std::cout << "+ " << usage << "\n";
        } else {
            while(time++ <= step.time) {
                std::cout << std::right << std::setw(W) << "|";
                for(auto u{0}; u<usage; ++u) {
                    std::cout << " ";
                }
                std::cout << "|\n";
            }
        }
    }
}

template <typename T>
std::ostream &CumulativeTimetabling<T>::display(std::ostream &os) const {
  os << "Cumulative Time-Tabling";

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
std::ostream &CumulativeTimetabling<T>::print_reason(std::ostream &os, const hint) const {
  os << "cumulative-time-tabling";
  return os;
}

} // namespace tempo

#endif
