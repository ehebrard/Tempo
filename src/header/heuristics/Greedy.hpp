
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
        job_map.resize(solver.numeric.size(), -1);
    }
    
    void addJobs(std::vector<Job<T>>& J) {
        jobs = J;
        unscheduled_jobs.reserve(jobs.size());
        unscheduled_jobs.fill();
        precedences.resize(jobs.size());
        int i{0};
        for(auto& j : jobs) {
            job_map[j.start.id()] = i;
            job_map[j.end.id()] = i;
            ++i;
        }
    }
    
    void addResource(std::vector<BooleanVar<>>::iterator bx, std::vector<BooleanVar<>>::iterator ex) {
        for(auto xi{bx}; xi!=ex; ++xi) {
            auto l{solver.boolean.getLiteral(true, *xi)};
            auto c{solver.boolean.getEdge(l)};
            precedences[job_map[c.to]].push_back(l);
            precedences[job_map[c.from]].push_back(~l);
        }
    }
    
    bool run();
    
private:
    Solver<T>& solver;
    std::vector<Job<T>> jobs;
    SparseSet<> unscheduled_jobs;
    std::vector<int> job_map;
    std::vector<std::vector<Literal<T>>> precedences;
    
};


template <typename T>
bool Greedy<T>::run() {
 
    solver.propagate();
    
    while(not unscheduled_jobs.empty()) {
        ++solver.num_choicepoints;
        
        int next{-1};
        for(auto j : unscheduled_jobs) {
            if(next == -1) {
                next = j;
            } else if(jobs[next].getEarliestStart(solver) > jobs[j].getEarliestStart(solver)) {
                next = j;
//            } else if(jobs[next].getEarliestStart(solver) == jobs[j].getEarliestStart(solver) and (random() % 2) == 1) {
//                next = j;
//            }
            } else if(jobs[next].getEarliestStart(solver) == jobs[j].getEarliestStart(solver) and jobs[next].getLatestEnd(solver) > jobs[j].getLatestEnd(solver)) {
                next = j;
            } else if(jobs[next].getEarliestStart(solver) == jobs[j].getEarliestStart(solver) and jobs[next].getLatestEnd(solver) == jobs[j].getLatestEnd(solver) and (random() % 2) == 1) {
                next = j;
            }
        }
        
        try {
            unscheduled_jobs.remove_back(next);
            for(auto p : precedences[next]) {
                if(solver.boolean.isUndefined(p.variable())) {
                    solver.set(p);
                }
            }

            solver.propagate();
            
        } catch(Failure<T>& f) {
            std::cout << "FAILED!\n";
            break;
        }
    }
    
    bool r{unscheduled_jobs.empty()};
    unscheduled_jobs.fill();
    return r;
}

}

#endif // __GREEDY_HPP
