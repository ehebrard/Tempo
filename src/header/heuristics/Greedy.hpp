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
        
//        std::cout << "post " << I.start << " <= " << best_start_time[i] << std::endl;
        
        solver.post(I.start <= best_start_time[i]);
        
//        std::cout << "post " << I.start << " >= " << best_start_time[i] << std::endl;
        
        solver.post(I.start >= best_start_time[i]);
        ++i;
    }
    
//    std::cout << "propagate\n";
    
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
//      Interval_map.resize(solver.numeric.size(), -1);
        precedences.resize(solver.numeric.size());
    }

    void addIntervals(std::vector<Interval<T>> &J) {
      Intervals = J;
      unscheduled_Intervals.reserve(Intervals.size());
      unscheduled_Intervals.fill();
//      precedences.resize(Intervals.size());
//      int i{0};
//      for (auto &j : Intervals) {
//        Interval_map[j.start.id()] = i;
//        Interval_map[j.end.id()] = i;
//        ++i;
//      }
        
//        for (auto &j : Intervals) {
//            std::cout << j << std::endl;
//        }
    }

    //    void addResource(const std::vector<DisjunctVar<T>>::iterator bx,
    //                     const std::vector<DisjunctVar<T>>::iterator ex) {
    //
    //      for (auto xi{bx}; xi != ex; ++xi) {
    //
    //        auto l{solver.boolean.getLiteral(true, *xi)};
    //        auto c{solver.boolean.getEdge(l)};
    //
    //        precedences[Interval_map[c.to]].push_back(l);
    //        precedences[Interval_map[c.from]].push_back(~l);
    //      }
    //    }

    void addResource(const std::vector<BooleanVar<T>>::iterator bx,
                     const std::vector<BooleanVar<T>>::iterator ex) {

      for (auto xi{bx}; xi != ex; ++xi) {
          
          addDisjunct(*xi);

//        auto l{solver.boolean.getLiteral(true, *xi)};
//        auto c{solver.boolean.getEdge(l)};
//
//        precedences[Interval_map[c.to]].push_back(l);
//        precedences[Interval_map[c.from]].push_back(~l);
      }
    }
    
    void addDisjunct(const BooleanVar<T> b) {
        addVar(b.id());
    }

    void addVar(const var_t b) {
      auto l{solver.boolean.getLiteral(true, b)};
      auto pc{solver.boolean.getEdge(l)};
        auto nc{solver.boolean.getEdge(~l)};
        
//        std::cout << pc << " OR " << nc << std::endl;
        
        
//      precedences[Interval_map[c.to]].push_back(l);
//      precedences[Interval_map[c.from]].push_back(~l);
        precedences[pc.to].push_back(l);
        precedences[nc.to].push_back(~l);
        
//        exit(1);
    }

    bool runEarliestStart();
    bool runLatestEnd();
    bool runLex();
    
private:
    Solver<T>& solver;
    std::vector<Interval<T>> Intervals;
    SparseSet<> unscheduled_Intervals;
//    std::vector<int> Interval_map;

    // for each numeric var, all its conditional precedences
    std::vector<std::vector<Literal<T>>> precedences;

    //    std::vector<SparseSet<>> unscheduled_Intervals_of;
};


template <typename T>
bool Greedy<T>::runLex() {
 
    solver.propagate();

    //    std::cout << solver << std::endl;

    while (not unscheduled_Intervals.empty()) {
      ++solver.num_choicepoints;

      int next = unscheduled_Intervals.front();
      unscheduled_Intervals.pop_front();

      try {

        //          std::cout << std::endl << "next="<< Intervals[next] <<
        //          std::endl;

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
//        if (unscheduled_Intervals.backsize() == 1)
//          solver.set(Intervals[next].end.before(
//              Intervals[next].getEarliestEnd(solver)));

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
//        for (auto p : precedences[next]) {
//          if (solver.boolean.isUndefined(p.variable())) {
//            solver.set(p);
//          }
//        }
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
//        for (auto p : precedences[next]) {
//          if (solver.boolean.isUndefined(p.variable())) {
//            solver.set(p);
//          }
//        }
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




}

#endif // __GREEDY_HPP
