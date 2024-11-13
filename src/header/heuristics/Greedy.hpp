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
#include <algorithm>
#include <cstddef>
#include <limits>
#include <numeric>
#include <Global.hpp>
#include <vector>


namespace tempo {


/**********************************************
 * Greedy Primal Heuristic : ScheduleGenerationScheme
 **********************************************/


template<typename T=int>
struct InflectionPoint {
    InflectionPoint() {}
    InflectionPoint(const T d, const T c) : date(d), capacity(c) {}
    
    T date{0}; // time point
    T capacity{0}; // remaining capacity until next timepoint
    
    std::ostream &display(std::ostream &os) const {
        os << capacity << "@" << date;
        return os;
    }
};


template<typename T=int>
struct InsertionPoint {
    T date;
    int index;
};


template<typename T=int>
class ResourceProfile {
  
public:
    ResourceProfile() = default;
    
    void initialise(const T C);
    
    // return the insertion point of the earliest element in profile whose date >= lb and such that the capacity remains positive or null
    InsertionPoint<T> earliest(const T demand, const T dur, const InsertionPoint<T> i={0,0}) ;
    void add_after(const InsertionPoint<T> i, const T demand, const T dur) ;
    void insert(const T demand, const T dur) { add_after(earliest(demand, dur), demand, dur); }
    void clear() { profile.clear(); }
    
    std::ostream &display(std::ostream &os) ;
    
private:
    T capacity{0};
    List<InflectionPoint<T>> profile;
    
};

template<typename T>
class ScheduleGenerationScheme {
    
public:
    
    ScheduleGenerationScheme(Solver<T> &solver,
                             std::vector<std::vector<size_t>>& tasks_requirements,
                             std::vector<std::vector<T>>& task_demands,
                             std::vector<T>& resource_capacities,
                             std::vector<Interval<T>>& intervals,
                             std::vector<std::pair<int,int>>& precedences);
    
    void clear();
    T run();
    void load();
    
    T best_makespan{Constant::Infinity<T>};
  
private:
    
    // Solver
    Solver<T> &solver;
    
    // data
    std::vector<std::vector<size_t>>& tasks_requirements;
    std::vector<std::vector<T>>& task_demands;
    std::vector<T>& resource_capacities;
    std::vector<Interval<T>>& intervals;
    std::vector<std::pair<int,int>>& precedences;
    
    BacktrackEnvironment env;
    
    
    //helpers
    std::vector<int> sources;
    SparseSet<> tasks;
    DirectedGraph<int> precedence_graph;
    std::vector<ResourceProfile<T>> profile;
    std::vector<InsertionPoint<T>> event;
    std::vector<T> start_time;
    std::vector<T> best_start_time;
    SparseSet<> stack;
    std::vector<T> dem;
    
    
public:
    // statistics
    unsigned long num_insertions{0};
    
};


template<typename T>
ScheduleGenerationScheme<T>::ScheduleGenerationScheme(Solver<T> &solver,
                         std::vector<std::vector<size_t>>& tasks_requirements,
                         std::vector<std::vector<T>>& task_demands,
                         std::vector<T>& resource_capacities,
                         std::vector<Interval<T>>& intervals,
                         std::vector<std::pair<int,int>>& precedences) :
solver(solver),
tasks_requirements(tasks_requirements),
task_demands(task_demands),
resource_capacities(resource_capacities),
intervals(intervals),
precedences(precedences)
, precedence_graph(&env)
{
    
    tasks.reserve(intervals.size());
    tasks.fill();
    
    
    precedence_graph.resize(intervals.size());
    
    int n{static_cast<int>(intervals.size())};
    for(auto prec : precedences) {
        auto x{prec.first};
        auto y{prec.second};
        if(x >= 0 and y<n) {
            if(tasks.has(y))
                tasks.remove_back(y);
            precedence_graph.add(x,y);
        }
    }
    
    for(auto t : tasks) {
        sources.push_back(t);
    }
    

    profile.resize(resource_capacities.size());
    for(size_t i{0}; i<resource_capacities.size(); ++i) {
        profile[i].initialise(resource_capacities[i]);
    }
    
    event.resize(resource_capacities.size());
    start_time.resize(intervals.size(), 0);
    
    stack.reserve(resource_capacities.size());
    dem.resize(resource_capacities.size(), 0);
}

template<typename T>
void ScheduleGenerationScheme<T>::clear() {
    for(auto &pr : profile)
        pr.clear();
    
    for(size_t i{0}; i<resource_capacities.size(); ++i) {
        profile[i].initialise(resource_capacities[i]);
    }
    
    std::fill(start_time.begin(), start_time.end(), 0);
    tasks.clear();
    for(auto t : sources)
        tasks.add(t);
}

template<typename T>
void ScheduleGenerationScheme<T>::load() {
    
    assert(not best_start_time.empty());
    
    auto s{solver.saveState()};
    int i{0};
    for(auto I : intervals) {
        
//        std::cout << i << "/" << best_start_time.size() << std::endl;

#ifdef DBG_RPROF
std::cout << "post " << I.start << " <= " << best_start_time[i] << std::endl;
#endif

solver.post(I.start <= best_start_time[i]);

#ifdef DBG_RPROF
std::cout << "post " << I.start << " >= " << best_start_time[i] << std::endl;
#endif

solver.post(I.start >= best_start_time[i]);
++i;
    }

#ifdef DBG_RPROF
    std::cout << "propagate\n";
#endif

    solver.propagate();
    solver.saveSolution();
    solver.restoreState(s);
    
    
//    std::cout << solver << std::endl;
}

template<typename T>
T ScheduleGenerationScheme<T>::run() {
    
    env.save();
    
 
    int makespan{0};
    while(not tasks.empty()) {
//        auto j{tasks.front()};
        auto j{tasks.any()};
        ++num_insertions;
        
#ifdef DBG_RPROF
        std::cout << "process " << j << ", dur = " << intervals[j].maxDuration(solver) <<  " min start-time = " << start_time[j] ; //<< std::endl;
#endif
        
        stack.clear();
        for(size_t i{0}; i<tasks_requirements[j].size(); ++i) {
            
            auto m{tasks_requirements[j][i]};
            event[m] = {start_time[j],0};
            
#ifdef DBG_RPROF
            std::cout << " req " << task_demands[j][i] << " of resource " << m ;
#endif
            
            stack.add(m);
            dem[m] = task_demands[j][i];
        }
        
#ifdef DBG_RPROF
        std::cout << std::endl;
#endif
        

        while(not stack.empty()) {
            auto m{stack.front()};
            auto d{dem[m]};
            
            event[m].date = std::max(start_time[j], event[m].date);
            event[m] = profile[m].earliest(d, intervals[j].maxDuration(solver), event[m]);
            
#ifdef DBG_RPROF
            std::cout << " insertion @" << event[m].date << " in resource " << m << "?\n";
#endif
            
            if(event[m].date > start_time[j]) {
                stack.setStart(0); // revise all resource
                start_time[j] = event[m].date;
            }
            
            stack.remove_front(m);
        }
            
        for(size_t i{0}; i<tasks_requirements[j].size(); ++i) {
            
            auto m{tasks_requirements[j][i]};
            auto d{dem[m]};
            
#ifdef DBG_RPROF
            std::cout << " -> insertion @" << event[m].date << " in resource " << m << "!\n";
#endif
            
            profile[m].add_after(event[m], d, intervals[j].maxDuration(solver));
        }
            
            
        auto finish_time_j{start_time[j] + intervals[j].maxDuration(solver)};
        
        makespan = std::max(makespan, finish_time_j);
        
            for(auto e : precedence_graph[j]) {
                start_time[e] = std::max(start_time[e], finish_time_j);
                if(precedence_graph.indegree(e) == 1 and tasks.isback(e)) {
                    tasks.add(e);
                }
            }
            precedence_graph.remove(j);
            tasks.remove_front(j);
        
#ifdef DBG_RPROF
        for(auto &pr : profile) {
            std::cout << pr << std::endl;
        }
        std::cout << tasks << std::endl;
#endif
    }
    
//    for(auto &pr : profile) {
//        std::cout << pr << std::endl;
//    }
    
    if(makespan < best_makespan) {
        
//        std::cout << "improved makespan " << best_makespan << " -> " << makespan << std::endl;
        
        best_makespan = makespan;
        best_start_time = start_time;
    }
    
    env.restore(0);
    
    return makespan;
}


template<typename T>
std::ostream &operator<<(std::ostream &os, const InflectionPoint<T> &x) {
    return x.display(os);
}

template<typename T>
std::ostream &operator<<(std::ostream &os, ResourceProfile<T> &x) {
    return x.display(os);
}

template<typename T>
void ResourceProfile<T>::initialise(const T C) {
    capacity = C;
    auto init{profile.create_element(0,capacity)};
    profile.add_front(init);
    auto horizon{profile.create_element(Constant::Infinity<T>,capacity)};
    profile.add_after(init,horizon);
    
}

template<typename T>
std::ostream &ResourceProfile<T>::display(std::ostream &os) {
    
#ifdef DBG_RPROF
    std::cout << profile << std::endl;
#endif
    
    auto l{profile.rbegin()};
    
//    std::cout << *l << std::endl;
    
    ++l;
    int horizon{std::max(10,l->date+1)};
    
//    std::cout << "horizon = " << horizon << std::endl;
    
    std::vector<std::vector<bool>> image(horizon);
    for(auto &row : image) {
        row.resize(capacity, false);
    }
    auto ti{profile.first()};
    while(ti != List<InflectionPoint<T>>::tail) {
        auto n{profile.next(ti)};
                if(profile[ti].capacity < capacity) {
                    
//                    std::cout << profile[ti].capacity << std::endl;
                    
                    auto e{std::min(horizon, profile[n].date)};
                    for(auto x{profile[ti].date}; x<e; ++x) {
                        
//                        if(capacity < profile[ti].capacity) {
//                            std::cout << "bug " << profile[ti].capacity << " > " << capacity << "\n";
//                            exit(1);
//                        }
                        
                        
//                        std::cout << x << "/" << image.size() << std::endl;
                        
                        
                        for(auto y{0}; y<(capacity - profile[ti].capacity); ++y) {
                            image[x][image[x].size() - y - 1] = true;
                        }
                    }
                }
        ti = n;
    }
    os << "+";
    for(auto x{0}; x<horizon; ++x) {
        os << "-";
    }
    os << "+\n";
    for(auto y{0}; y<capacity; ++y) {
        os << "|";
        for(auto x{0}; x<horizon; ++x) {
            os << (image[x][y] ? "X" : " ");
        }
        os << "|" << std::endl;
    }
    os << "+";
    for(auto x{0}; x<horizon; ++x) {
        os << "-";
    }
    os << "+\n";
    return os;
}

template<typename T>
InsertionPoint<T> ResourceProfile<T>::earliest(const T demand, const T dur, const InsertionPoint<T> ip) {

    auto lb{ip.date};
    auto i{ip.index};
    auto ti{profile.begin()};
    if(i > 0)
        ti = profile.at(i);
    auto succ{ti};
    while(ti!=profile.end()) {
        
#ifdef DBG_RPROF
        std::cout << *ti << "?";
#endif
      
        ++succ;
        
        if(succ->date > lb) {
            if(ti->capacity >= demand) {
                auto finish{std::max(lb,ti->date) + dur};
                auto next{ti};
                bool fits{true};
                ++next;
                while(fits and next->date < finish) {
                    fits = (next->capacity >= demand);
                    ++next;
                }
                if(fits) {
#ifdef DBG_RPROF
                    std::cout << " yes!\n";
#endif
                    return {std::max(lb,ti->date), ti.index};
                }
#ifdef DBG_RPROF
                else {
                    std::cout << " no\n";
                }
#endif
            }
        }
        
        ti = succ;
    }
    return {0,0};
}

template<typename T>
void ResourceProfile<T>::add_after(const InsertionPoint<T> i, const T demand, const T duration) {
    
    
    auto t{i.index};
    auto lb{i.date};
    
    
#ifdef DBG_RPROF
    std::cout << "add dem=" << demand << " dur=" << duration << " after " << profile[t] << std::endl;
#endif
    
        auto finish{profile[t].date + duration};
        auto tp{profile.at(t)};
        auto next{tp};
    
    if(lb <= tp->date) {
        auto prev{tp};
        --prev;
        if(prev != profile.end() and prev->capacity == tp->capacity - demand) {
            
#ifdef DBG_RPROF
            std::cout << " remove " << *tp << std::endl;
#endif
            
            ++next;
            profile.remove(tp.index);
        }
    } else {
        
#ifdef DBG_RPROF
            std::cout << " start at " << lb << " instead of " << tp->date << std::endl;
#endif
        
            ++next;
            finish = lb + duration;
        
            auto j{profile.create_element(lb, tp->capacity - demand)};
            profile.add_after(tp.index, j);
    }
    
        while(next->date < finish) {
            
#ifdef DBG_RPROF
            std::cout << " reduce capa of " << *next << " by " << demand << std::endl;
#endif
            
            next->capacity -= demand;
            ++next;
        }
        
#ifdef DBG_RPROF
        std::cout << " stop before " << *next << std::endl;
#endif
    
        if(next->date > finish) {
            --next;
            
#ifdef DBG_RPROF
            std::cout << " create a new event with the capa of " << *next << std::endl;
#endif
            
            auto j{profile.create_element(finish, next->capacity + demand)};
            profile.add_after(next.index, j);
        } else {
            if(next->date < finish) {
                std::cout << "bug!\n";
                exit(1);
            }
            
            auto prev{next};
            --prev;
            if(next->capacity == prev->capacity) {
                
#ifdef DBG_RPROF
                std::cout << " remove " << *next << std::endl;
#endif
                
                profile.remove(next.index);
            }
        }
}



/**********************************************
 * Greedy Primal Heuristic for Disjunctive Problems
 **********************************************/

template <typename T> class Greedy  {
public:
    
    Greedy(Solver<T>& s) : solver(s) {
        precedences.resize(solver.numeric.size());        
        prec_map_ptrs.resize(solver.numeric.size());
    }

    void addIntervals(std::vector<Interval<T>> &J) {
      Intervals = J;
      unscheduled_Intervals.reserve(Intervals.size());
      unscheduled_Intervals.fill();
    }


    void addResource(const std::vector<BooleanVar<T>>::iterator bx,
                     const std::vector<BooleanVar<T>>::iterator ex) {

      for (auto xi{bx}; xi != ex; ++xi) {          
          addDisjunct(*xi);
      }
    }
    
    void addDisjunct(const BooleanVar<T> b) {
        addVar(b.id());
    }

    void addVar(const var_t b) {

        // std::cout << "--" << b << "--" << std::endl;
        auto l{solver.boolean.getLiteral(true, b)};
        // std::cout << l << std::endl;
        // std::cout << solver.boolean.getEdge(l) << std::endl;
        auto pc{solver.boolean.getEdge(l)};
        auto nc{solver.boolean.getEdge(~l)};

        // std::cout << precedences.size() << " - " << pc.to << std::endl;
        // std::cout << "("<< solver.numeric.size() << ")" << std::endl;

        if (pc != Constant::NoEdge<T>) {
            precedences[pc.to].push_back(l);
        }

        if(nc != Constant::NoEdge<T>) {
            precedences[nc.to].push_back((~l));
        }


        // std::cout << "prec maps " <<  pc.to << " - " << pc.from << std::endl;

        if(prec_map_ptrs[pc.to].size() <= pc.from){
            prec_map_ptrs[pc.to].resize(pc.from+1, nullptr);
        }
        prec_map_ptrs[pc.to][pc.from] = &precedences[pc.to].back();
        // std::cout << *prec_map_ptrs[pc.to][pc.from] << " - " << precedences[pc.to].back() << std::endl;

        // std::cout << "-----------" << std::endl;
        // for(auto k{0u}; k < prec_map_ptrs.size(); ++k){
        //     for(auto l{0u}; l < prec_map_ptrs[k].size(); ++l){
        //         if(prec_map_ptrs[k][l] != nullptr){
        //             std::cout << *prec_map_ptrs[k][l] << " ";
        //         }else{
        //             std::cout << prec_map_ptrs[k][l] << " ";
        //         }
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout << "-----------" << std::endl;
    }

    bool runEarliestStart();
    bool runLatestEnd();
    bool runLex();
    bool runOrienteering(std::vector<int> profits, std::vector<std::vector<T>> transition);
    
private:
    Solver<T>& solver;
    std::vector<Interval<T>> Intervals;
    SparseSet<> unscheduled_Intervals;
    std::vector<std::list<Literal<T>>> precedences;
    std::vector<std::vector<Literal<T>*>> prec_map_ptrs;

};


template <typename T>
bool Greedy<T>::runLex() {
 
    solver.propagate();


    while (not unscheduled_Intervals.empty()) {
      ++solver.num_choicepoints;

      int next = unscheduled_Intervals.front();
      unscheduled_Intervals.pop_front();

      try {

        unscheduled_Intervals.remove_back(next);
        for (auto p : precedences[Intervals[next].start.id()]) {
          if (solver.boolean.isUndefined(p.variable())) {
            //            std::cout << " -> " << solver.pretty(p) << std::endl;
            solver.set(p);
          }
        }
          for (auto p : precedences[Intervals[next].end.id()]) {
            if (solver.boolean.isUndefined(p.variable())) {
              //              std::cout << " -> " << solver.pretty(p) <<
              //              std::endl;
              solver.set(p);
            }
          }

        solver.propagate();

        //                    std::cout << solver << std::endl;

      } catch (Failure<T> &f) {
        //                    std::cout << "FAILED!\n";
        //          exit(1);
        break;
      }
    }

    bool r{unscheduled_Intervals.empty()};
    unscheduled_Intervals.fill();
    return r;
}



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

        //                    std::cout << std::endl << "next="<<
        //                    Intervals[next] << std::endl;

        unscheduled_Intervals.remove_back(next);

          for (auto p : precedences[Intervals[next].start.id()]) {
            if (solver.boolean.isUndefined(p.variable())) {

              //                std::cout << " -> " << solver.pretty(p) <<
              //                std::endl;

              solver.set(p);
            }
          }
          for (auto p : precedences[Intervals[next].end.id()]) {
            if (solver.boolean.isUndefined(p.variable())) {

              //                std::cout << " -> " << solver.pretty(p) <<
              //                std::endl;

              solver.set(p);
            }
          }
        if (unscheduled_Intervals.backsize() == 1)
          solver.set(Intervals[next].end.before(
              Intervals[next].getEarliestEnd(solver)));

        solver.propagate();

        //            std::cout << solver << std::endl;

      } catch (Failure<T> &f) {
//                    std::cout << "FAILED!\n";
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

          for (auto p : precedences[Intervals[next].start.id()]) {
            if (solver.boolean.isUndefined(p.variable())) {
              solver.set(p);
            }
          }
          for (auto p : precedences[Intervals[next].end.id()]) {
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

template<typename T>
bool Greedy<T>::runOrienteering(std::vector<int> profits, std::vector<std::vector<T>> transition) {    
    // std::cout << "In greedy" << std::endl;
    // Sort intervalls by profit
    assert(profits.size() == Intervals.size());
    std::vector<T> releases;
    releases.reserve(profits.size());
    std::vector<T> dues;
    dues.reserve(profits.size());
    std::vector<T> _est;
    _est.resize(profits.size());
    for(auto elt : Intervals) {
        releases.push_back(elt.getEarliestStart(solver));
        dues.push_back(elt.getLatestEnd(solver));
    }

    std::vector<index_t> idxs(profits.size());
    std::iota(idxs.begin(), idxs.end(), 0);

    std::sort(idxs.begin(), idxs.end(), [&profits](const int &a, const int &b){return profits[a] > profits[b];});

    std::vector<index_t> next(idxs.size()+2);
    std::vector<index_t> prev(idxs.size()+2);
    std::vector<T> cumul_gap_after(idxs.size()+2, 0);
    std::vector<T> min_marge_after(idxs.size()+2, 0);    
    const index_t FIRST{static_cast<index_t>(next.size()-1)};
    const index_t LAST{static_cast<index_t>(next.size()-2)};
    
    for(auto  i{0u}; i < next.size()-2; ++i) {
        next[i] = i;
        prev[i] = i;
    }
    next[FIRST] = LAST;
    next[LAST] = FIRST;
    prev[FIRST] = LAST;
    prev[FIRST] = LAST;
    min_marge_after[FIRST] = std::numeric_limits<T>::max();
    min_marge_after[LAST] = std::numeric_limits<T>::max();

    auto est = [&_est] (index_t idx) { return _est[idx]; };
    auto ect = [&_est, this] (index_t idx) { return _est[idx] + Intervals[idx].minDuration(solver); };
    auto slack = [&next, &dues, &ect] (index_t idx) { return next[idx] != idx ? dues[idx] - ect(idx) : 0; };
    auto gap = [&next, &prev, &FIRST, &LAST, &transition, &est, &ect] (index_t idx) { 
        if(next[idx] == idx || prev[idx] == FIRST || idx == FIRST || idx == LAST) {
            return 0;
        }
        //std::cout << " Gap of " << idx << " : " << est(idx) << " - " << ect(prev[idx]) << " - " << transition[prev[idx]][idx] << " = " << (est(idx) - ect(prev[idx]) - transition[prev[idx]][idx]) << std::endl;
        return est(idx) - ect(prev[idx]) - transition[prev[idx]][idx];
    };
    auto gap_between = [&next, &cumul_gap_after] (index_t idx1, index_t idx2) {
        if(next[idx1] == idx1 ||  next[idx2] == idx2) {
            return 0;
        }
        return cumul_gap_after[idx1] - cumul_gap_after[idx2];
    };

    auto marge_insertion_after = [&next, &gap_between, &slack](index_t idx1, index_t idx2) {
        if(next[idx1] == idx1 ||  next[idx2] == idx2) {
            return 0;
        }
        return gap_between(idx1, idx2) + slack(idx2);
    };
    // std::cout << "End initialisation" << std::endl;

    // std::cout << "Insertion of first element" << std::endl;
    // exist to true + propagate and choice point        
    // first element    
    next[FIRST] = idxs.front();
    prev[LAST] = idxs.front();
    next[idxs.front()] = LAST;
    prev[idxs.front()] = FIRST;
    cumul_gap_after[idxs.front()] = 0;
    min_marge_after[FIRST] = slack(idxs.front());
    min_marge_after[idxs.front()] = min_marge_after[LAST];
    _est[idxs.front()] = releases[idxs.front()];

    // {
    //     auto k = next[FIRST];
    //     std::cout << "Sequence is" << std::endl;
    //     std::cout << "[" << cumul_gap_after[FIRST] << " - " << min_marge_after[FIRST] << "]" << std::endl;
    //     while(k != LAST){
    //         std::cout << k << "(" << est(k) << ") - " << slack(k) << std::endl;
    //         std::cout << cumul_gap_after[k] << " - " << min_marge_after[k] << std::endl;
    //         k = next[k];
    //     }
    //     std::cout << "[" << cumul_gap_after[LAST] << " - " << min_marge_after[LAST] << "]" << std::endl;

    // }
    


    // std::cout << "Start loop" << std::endl;
    // std::cout << "First = " << FIRST << std::endl;
    // std::cout << "Last = " << LAST << std::endl;

    auto i{1u};
    while(i < idxs.size() && profits[idxs[i]] > 0) {
        // find insertion
        auto idx{idxs[i]};
        std::vector<index_t> best_insertion_idxs;
        T best_insertion_value{std::numeric_limits<T>::max()};

        // std::cout << "Try to insert " << idx << std::endl;

        // Test insert front
        
        auto new_est{std::max(releases[next[FIRST]], releases[idx] + Intervals[idx].minDuration(solver) + transition[idx][next[FIRST]])};
        auto shift{new_est - releases[next[FIRST]]};
        if(shift < min_marge_after[next[FIRST]]) {
            best_insertion_value = transition[idx][next[FIRST]];
            best_insertion_idxs.push_back(FIRST);
        }
        
        

        // std::cout << "Front insertion is " << (best_insertion_value<std::numeric_limits<T>::max()) << std::endl;



        for(auto idx_before{next[FIRST]}; idx_before != LAST; ){
            bool valid{true};

            // std::cout << "Try inserting after " << idx_before << std::endl;

            auto actual_end_transition{ect(idx_before)};            
            if(next[idx_before] != LAST){ //last element
                actual_end_transition += transition[idx_before][next[idx_before]];                
            }

            // std::cout << " Actual end_transition is " << actual_end_transition << std::endl;
            auto est_inserted{std::max(releases[idx], ect(idx_before) + transition[idx_before][idx])};

            auto new_end_transition{est_inserted + Intervals[idx].minDuration(solver)};
            if(next[idx_before] != LAST){
                new_end_transition += transition[idx][next[idx_before]];
            }
            // std::cout << " New end transition is " << new_end_transition << std::endl;
            auto shift{new_end_transition - actual_end_transition};

            // std::cout << "Marge after is " <<  min_marge_after[idx_before] << std::endl;
            valid = est_inserted + Intervals[idx].minDuration(solver) < dues[idx] && shift <= min_marge_after[idx_before];
            
            if(valid) {
                // std::cout << "Is valid" << std::endl;
                auto delta_trans{transition[idx_before][idx]};
                if(next[idx_before] != LAST){
                    delta_trans += transition[idx][next[idx_before]] - transition[idx_before][next[idx_before]];
                }
                // std::cout << "Delta is " << delta_trans << std::endl;
                if(delta_trans <= best_insertion_value){
                    // std::cout << "Best" << std::endl;
                    if(delta_trans < best_insertion_value){
                        // std::cout << "  -> New best" << std::endl;
                        best_insertion_idxs.clear();
                        best_insertion_value = delta_trans;
                    }                                        
                    best_insertion_idxs.push_back(idx_before);
                }
            } 
            // else{
            //     std::cout << "Is not valid" << std::endl;
            // }
            idx_before = next[idx_before];
        }

        // insert and update
        if(best_insertion_value < std::numeric_limits<T>::max()) {            
            auto best_idx{best_insertion_idxs[random()%best_insertion_idxs.size()]};
            // std::cout << "Insertion in sequence after " << best_idx << std::endl;
            next[idx] = next[best_idx];
            next[best_idx] = idx;
            prev[idx] = prev[next[idx]];
            prev[next[idx]] = idx;

            auto k = next[FIRST];            
            _est[k] = releases[k];
            while(next[k] != LAST) {
                _est[next[k]] = std::max(releases[next[k]], ect(k) + transition[k][next[k]]);
                k = next[k];
            }


            k = prev[LAST];
            cumul_gap_after[k] = 0;
            min_marge_after[k] = std::numeric_limits<T>::max();
            while(k != FIRST) {
                cumul_gap_after[prev[k]] = cumul_gap_after[k] + gap(k);
                // std::cout << k << " : " << std::endl;
                // std::cout << "Min marg after  " << min_marge_after[k] << ", gap is " << gap(k) << ", marge insertion after " << marge_insertion_after(prev[k], k) << "(" << slack(k) << ", " << est(k) << "/" << ect(k)  <<")" << std::endl;;
                if(min_marge_after[k] == std::numeric_limits<int>::max()){
                    min_marge_after[prev[k]] = marge_insertion_after(prev[k], k);
                }else{
                    min_marge_after[prev[k]] = std::min(min_marge_after[k] + gap(k), marge_insertion_after(prev[k], k));
                }
                
                
                k = prev[k];
            }            
            // {
            //     k = next[FIRST];
            //     std::cout << "Sequence is" << std::endl;
            //     std::cout << "[" << cumul_gap_after[FIRST] << " - " << min_marge_after[FIRST] << "]" << std::endl;
            //     while(k != LAST){
            //         std::cout << k << "(" << est(k) << ", " << ect(k) << ") - " << slack(k) << std::endl;
            //         std::cout << "[" << cumul_gap_after[k] << " - " << min_marge_after[k] << "]" << std::endl;
            //         k = next[k];
            //     }
            //     std::cout << "[" << cumul_gap_after[LAST] << " - " << min_marge_after[LAST] << "]" << std::endl;
            // }
        }
        ++i;
        // std::cout << " ------- " << std::endl;
    }

    
    
    
    auto k{next[FIRST]};
    // {
    //     auto obj{0};
    //     while(k != LAST) {
    //         obj+= profits[k];
    //         std::cout << k << " [" << releases[k] << "," << dues[k] << "]-" << Intervals[k].minDuration(solver);
    //         if(next[k] != LAST)
    //             std::cout << " => " << transition[k][next[k]];
    //         std::cout << std::endl;
    //         std::cout << est(k) << std::endl;
    //         k = next[k];
    //     }
    //     std::cout << "Profit = " << obj << std::endl;
    // }

    k = next[FIRST];
    while(k != LAST) {
        // std::cout << "Interval " << k << std::endl;
        solver.set(solver.boolean.getLiteral(true, Intervals[k].exist));
        auto sid{Intervals[k].start.id()};
        // auto eid{Intervals[k].end.id()};
        auto l{next[k]};
        while(l != LAST) {
            auto sid2{Intervals[l].start.id()};
            // auto eid2{Intervals[l].end.id()};
            // std::cout << " -> " << l << std::endl;
            // std::cout << "      prec  " <<  sid << " - " << sid2  << " : ";
            // if(prec_map_ptrs[sid][sid2] != nullptr)
            //     std::cout << *prec_map_ptrs[sid][sid2] << " - " << solver.boolean.isUndefined(prec_map_ptrs[sid][sid2]->variable()) << std::endl;
            // else
            //     std::cout << "None" << std::endl;

            if(prec_map_ptrs[sid][sid2] != nullptr && solver.boolean.isUndefined(prec_map_ptrs[sid][sid2]->variable())){
                // std::cout << "      Set the precedence " << solver.boolean.getEdge(*prec_map_ptrs[sid][sid2]) << std::endl;
                // std::cout << "      " <<  Intervals[k].getEarliestStart(solver) << "-" << Intervals[k].getLatestStart(solver) << "   ---   " << Intervals[l].getEarliestStart(solver) << "-" << Intervals[l].getLatestStart(solver) << std::endl;
                solver.set(*prec_map_ptrs[sid][sid2]);
            }
            l = next[l];
        }
        k = next[k];
    }
    solver.propagate();


    // k = next[FIRST];
    // while(k != LAST) {
    //     std::cout << k <<"[" << releases[k] << ","<<dues[k] << "]-" << Intervals[k].minDuration(solver) <<"  =>   " <<  Intervals[k].getEarliestStart(solver) << "-" << Intervals[k].getLatestStart(solver) << std::endl;
    //     k = next[k];
    // }

    return true;
}


}

#endif // __GREEDY_HPP
