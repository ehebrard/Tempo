/**
* @author Emmanuel Hebrard
* @date 12.09.24
* @brief
*/

#ifndef TEMPO_LRBHEAP_HPP
#define TEMPO_LRBHEAP_HPP


//#include <algorithm>

#include "Solver.hpp"
#include "heuristics/impl/LearningRateMap.hpp"
#include "heuristics/impl/ActivityMap.hpp"
#include "util/Heap.hpp"
#include "util/SparseSet.hpp"
#include "util/traits.hpp"
#include "util/SubscribableEvent.hpp"

#include "heuristic_interface.hpp"
#include "ReversibleObject.hpp"

//#define DBG_LearningRate true

namespace tempo {
template <typename T> class Solver;
}

namespace tempo::heuristics {
/**
 * @brief Random variable selection strategy
 * @details @copybrief Randomly chooses a variable from the remenaing search
 * variables
 * @note Right now, only binary variables are selected
 */
template<typename T>
class LearningRateHeap : public ReversibleObject, public BaseVariableHeuristic<T> {
    
public:
//    template <concepts::scalar T>
    LearningRateHeap(Solver<T> &solver)
    : ReversibleObject(&(solver.getEnv()))
    , m_solver(solver)
    , conflictCallback(solver.ConflictExtracted.subscribe_handled(
                                                              [this](const auto &arg) { this->incrementActivity(arg); }))
    , propagationCallback(solver.PropagationCompleted.subscribe_handled(
                                                              [this](const auto &arg) { this->resetActivity(arg); }))
    , boolean_rate(solver.getBooleanLearningRate())
    , numeric_rate(solver.getNumericLearningRate())
    {
        
        var_heap.resize(solver.boolean.size());
        std::iota(var_heap.begin(), var_heap.end(), 0);
        index.resize(solver.boolean.size(), 0);
        trail.push_back(var_heap.size());
        
        //@TODO: better initialisation
        for (size_t i{0}; i < var_heap.size(); ++i) {
            std::swap(var_heap[i], var_heap[i + random() % (var_heap.size() - i)]);
            index[var_heap[i]] = static_cast<index_t>(i);
        }
        
        boolean_rate.resize(solver.boolean.size(), impl::ActivityMap::baseIncrement);
        numeric_rate.resize(solver.numeric.size(), impl::ActivityMap::baseIncrement);

        assert(checkVars(solver));
    }
    
    LearningRateHeap(const LearningRateHeap &) = delete;
    LearningRateHeap(LearningRateHeap &&) = delete;
    LearningRateHeap &operator=(const LearningRateHeap &) = delete;
    LearningRateHeap &operator=(LearningRateHeap &&) = delete;
    ~LearningRateHeap() override = default;
    
    /**
     * Heuristic interface
     * @tparam T timing type
     * @param solver solver for which to select the variable
     * @return randomly selected variable
     */
    auto nextVariable(const Solver<T> &solver) noexcept -> VariableSelection override {
        
//#ifdef DBG_LearningRate
//        std::cout << "beg next var\n";
//#endif
        
        const concepts::same_template<SparseSet> auto &variables =
        solver.getBranch();
        
//#ifdef DBG_LearningRate
//        std::cout << "pick var from " << variables << " @lvl " << solver.level() << " (self lvl = " << trail.size() << ")\nheap:";
//        for (auto v{var_heap.begin()}; v != (var_heap.begin() + trail.back()); ++v) {
//            std::cout << " " << *v;
//        }
//        std::cout << std::endl;
//#endif
        
        assert(not variables.empty());
        
        
#ifdef DBG_LearningRate
        if(static_cast<int>(trail.size()) < solver.level())
            std::cout << " - catch up from lvl " << trail.size() << " to lvl " << solver.level() << std::endl;
#endif
        
//        std::cout << " lvls " << trail.size() << " " << solver.level() << " " << this->local_env->level() << std::endl;
//        std::cout << static_cast<size_t>(&(solver.getEnv())) << " " << static_cast<size_t>(env) << std::endl;
        
        while (static_cast<int>(trail.size()) < solver.level())
            trail.push_back(trail.back());
        
        assert(checkVars(solver));
        assert(checkHeap(0));
        
        auto last{trail.back()};
        var_t x;
        do {
            x = pickBest(last);
        } while (not variables.has(x));
        
        
//        std::cout << x << ": " << boolean_rate[x] << " / " << (boolean_rate[x] - boolean_rate[var_heap[0]]) << std::endl;
        
        trail.push_back(last);
        this->save();
     
        assert(checkHeap(0));
        
        return {x, VariableType::Boolean};
    }
    
    void undo() override {
        
        
        
//        ++backtracks;
        
#ifdef DBG_LearningRate
        std::cout << "(r) backtrack (" << this->local_env->level() << ")\n";
#endif
        
        updateRate(m_solver);
        
    }
    
    
protected:
    
//    template <concepts::scalar T>
    void resetActivity(const Solver<T> &solver) {
        
//#ifdef DBG_LearningRate
//        std::cout << "beg reset\n";
//#endif
        
        bool_buffer.clear();
        num_buffer.clear();
        
        while(litPointer < solver.numLiteral()) {
            
//            std::cout << litPointer << " / " << solver.numLiteral() << std::endl;
            
            
            auto l{solver.getLiteral(litPointer++)};
            
//            std::cout << l << std::endl;
            
            
            if(l.isNumeric()) {
                num_buffer.push_back(l.variable());
            } else {
                bool_buffer.push_back(l.variable());
            }
            
//            std::cout << bool_buffer.size() << " | " << num_buffer.size() << std::endl;
        }
        
        
//#ifdef DBG_LearningRate
//        std::cout << "mid reset\n";
//#endif
        
//
//        auto &variables{solver.getBranch()};
//        
//        
//#ifdef DBG_LearningRate
//            std::cout << " - (before reset) var ptr " << varPointer << " / " << variables.start_idx() << std::endl;
//#endif
//        
//        while(varPointer < variables.start_idx()) {
//            bool_buffer.push_back(variables[varPointer++]);
//            
//#ifdef DBG_LearningRate
//            std::cout << " - reset rate of " << bool_buffer.back() << std::endl;
//#endif
//            
//        }
//        
//        
//#ifdef DBG_LearningRate
//            std::cout << " - (after reset) var ptr " << varPointer << " / " << variables.start_idx() << std::endl;
//#endif
        
//        10 -- x <= 100
//        50 -- x <= 70           3     3/40
//        90 -- x <= 25           10  + 7/40
//        140 -- undo x <= 25     12  + 2/50
//        200 -- undo x <= 70     12  + 0
//        250 -- x <= 30          20  + 8/50
//        400 -- undo x <= 30     30
//        500 -- undo x <= 100    33
        
        
        
        numeric_rate.update(num_buffer, solver.num_fails);
        numeric_rate.resetActivity(num_buffer, solver.num_fails);
        
#ifdef DBG_LearningRate
        for(auto v : num_buffer) {
            std::cout << " - update rate of x" << v << " (" << numeric_rate.participated[v] << "/" << (solver.num_fails - numeric_rate.assigned_at[v]) << ") and reset" << std::endl;
        }
#endif
        
#ifdef DBG_LearningRate
        for(auto v : bool_buffer)
            std::cout << " - reset rate of b" << v << std::endl;
#endif

        boolean_rate.resetActivity(bool_buffer, solver.num_fails);
        
//#ifdef DBG_LearningRate
//        std::cout << "end reset\n";
//#endif
    }
    

    void updateRate(const Solver<T> &solver) {
        
#ifdef DBG_LearningRate
        std::cout << "beg update after backtrack " << litPointer << "/" << solver.numLiteral() << "\n";
#endif
        
//        bool_buffer.clear();
        num_buffer.clear();
        
        auto n{trail.back()};
        while(litPointer > solver.numLiteral()) {
            auto l{solver.getLiteral(--litPointer)};
            
            
//#ifdef DBG_LearningRate
//            std::cout << " - lit " << l << std::endl;
//#endif
            
            if(l.isNumeric()) {
                num_buffer.push_back(l.variable());
            } else {
                auto v{l.variable()};
//                bool_buffer.push_back(v);
                
                auto r{boolean_rate[v]};
                
#ifdef DBG_LearningRate
            std::cout << " - update rate of b" << v << std::endl;
#endif
                
                boolean_rate.update(v, solver.num_fails);
    
                if(index[v] < n) {
                    if(r < boolean_rate[v]) {
                        heap::percolate_up(var_heap.begin(), index[v], index,
                                           [&](const var_t x, const var_t y) {
                            return boolean_rate[x] > boolean_rate[y];
                        });
                    } else {
                        heap::percolate_down(var_heap.begin(), var_heap.begin()+n, index[v], index,
                                           [&](const var_t x, const var_t y) {
                            return boolean_rate[x] > boolean_rate[y];
                        });
                    }
                }
            }
        }
        
#ifdef DBG_LearningRate
        for(auto v : num_buffer)
            std::cout << " - update rate of x" << v << " (" << numeric_rate.participated[v] << "/" << (solver.num_fails - numeric_rate.assigned_at[v]) << ") and reset" << std::endl;
#endif
        
        numeric_rate.update(num_buffer, solver.num_fails);
        numeric_rate.resetActivity(num_buffer, solver.num_fails);
        
        
//        auto &variables{solver.getBranch()};
//        
//#ifdef DBG_LearningRate
//            std::cout << " - (before update) var ptr " << varPointer << " / " << variables.start_idx() << std::endl;
//#endif
//      
//        auto n{trail.back()};
//        while(varPointer > variables.start_idx()) {
//            auto v{variables[--varPointer]};
//            bool_buffer.push_back(v);
//            
//#ifdef DBG_LearningRate
//            std::cout << " - update rate of " << bool_buffer.back() << std::endl;
//#endif
//            
//            auto r{boolean_rate[v]};
//            boolean_rate.update(v, solver.num_fails);
//            
//            if(index[v] < n) {
//                if(r < boolean_rate[v]) {
//                    heap::percolate_up(var_heap.begin(), index[v], index,
//                                       [&](const var_t x, const var_t y) {
//                        return boolean_rate[x] > boolean_rate[y];
//                    });
//                } else {
//                    heap::percolate_down(var_heap.begin(), var_heap.begin()+n, index[v], index,
//                                       [&](const var_t x, const var_t y) {
//                        return boolean_rate[x] > boolean_rate[y];
//                    });
//                }
//            }
//        }
        
        
//        
//#ifdef DBG_LearningRate
//            std::cout << " - (after update) var ptr " << varPointer << " / " << variables.start_idx() << std::endl;
//#endif
 
        trail.pop_back();

        
        while (n < trail.back()) {
            
#ifdef DBG_LearningRate
            std::cout << " - restore and percolate var " << var_heap[n] << "\n";
//            assert(variables.has(var_heap[n]));
#endif
            
            heap::percolate_up(var_heap.begin(), n, index,
                               [&](const var_t x, const var_t y) {
                return boolean_rate[x] > boolean_rate[y];
            });
            ++n;
        }
        
//        
//#ifdef DBG_LearningRate
//            std::cout << " end backtrack" << std::endl;
//#endif
        
    }
    
    
//    template <concepts::scalar T>
    void incrementActivity(const Solver<T> &solver) {
        
//#ifdef DBG_LearningRate
//        std::cout << "beg increment\n";
//#endif
        
        bool_buffer.clear();
        num_buffer.clear();
        for (auto l : solver.lastLearnt()) {
            if (l.isNumeric()) {
                num_buffer.push_back(l.variable());
            } else {
                bool_buffer.push_back(l.variable());
            }
            if (weight_reasons) {
                for(auto p : solver.getExplanation(~l)) {
                    if (p.isNumeric()) {
                        num_buffer.push_back(p.variable());
                    } else {
                        bool_buffer.push_back(p.variable());
                    }
                }
            }
        }
        
        for (auto i : solver.cut.cached_) {
            auto l{solver.getLiteral(i)};
            if (l.isNumeric()) {
                num_buffer.push_back(l.variable());
            } else {
                bool_buffer.push_back(l.variable());
            }
        }
        
        boolean_rate.incrementActivity(bool_buffer);
        numeric_rate.incrementActivity(num_buffer);
        
        assert(checkHeap(0));
    }
    
//    template <concepts::scalar T>
    bool checkVars(const Solver<T> &solver) {
        auto last{trail.back()};
        for (auto x : solver.getBranch()) {
            bool notin{true};
            for (auto v{var_heap.begin()}; notin and v != (var_heap.begin() + last);
                 ++v) {
                if (*v == x) {
                    notin = false;
                }
            }
            if (notin) {
                std::cout << x << " is not in the heap!\n";
                if(index.size() > static_cast<size_t>(x)) {
                    std::cout << "(@index "<< index[x] << "/" << last << ")\n";
                }
                return false;
            }
        }
        return true;
    }
    
    bool checkHeap(const int i) {
        auto lc{heap::left(i)};
        auto rc{heap::right(i)};
        
        auto n{static_cast<int>(trail.back())};
        
        bool ok_left =
        ((lc >= n) or
         ((boolean_rate[var_heap[i]] >= boolean_rate[var_heap[lc]]) and checkHeap(lc)));
        
        if (not ok_left) {
            std::cout << lc << "|" << var_heap[lc] << ":(" << boolean_rate[var_heap[lc]]
            << ") is the left child of " << i << "|" << var_heap[i] << ":("
            << boolean_rate[var_heap[i]] << ")\n";
        }
        
        if (ok_left) {
            auto ok_right =
            ((rc >= n) or ((boolean_rate[var_heap[i]] >= boolean_rate[var_heap[rc]]) and
                           checkHeap(rc)));
            
            if (not ok_right) {
                std::cout << rc << "|" << var_heap[rc] << ":(" << boolean_rate[var_heap[rc]]
                << ") is the right child of " << i << "|" << var_heap[i]
                << ":(" << boolean_rate[var_heap[i]] << ")\n";
            }
            
            return ok_right;
        }
        
        return false;
    }
    
    var_t pickBest(size_t &_end_) {
        var_t x{var_heap[0]};
        heap::remove_min(var_heap.begin(), var_heap.begin() + _end_, index,
                         [&](const var_t x, const var_t y) {
            return boolean_rate[x] > boolean_rate[y];
        });
        --_end_;
        
#ifdef DBG_LearningRate
        std::cout << " - pick " << x << " (" << _end_ << ")\n";
#endif
        
        return x;
    }
    
    Solver<T> &m_solver;
    
    // the successive sizes of the heap -> decreasing at each call to
    // nextVariable() |trail| should therefore be equal to solver.level(), if not,
    // there has been a backtrack then the values in trail can be used to
    // re-integrate the necessary variables
    std::vector<size_t> trail;
    
    // position of each variable in the heap (so that we can percolate)
    std::vector<index_t> index;
    
    std::vector<var_t> var_heap;
    
    SubscriberHandle conflictCallback;
    SubscriberHandle propagationCallback;
//    SubscriberHandle backtrackCallback;
    
    std::vector<var_t> bool_buffer;
    std::vector<var_t> num_buffer;
    
    impl::LearningRateMap& boolean_rate;
    impl::LearningRateMap& numeric_rate;
//    impl::ActivityMap& num_activity; // learning rate for numeric vars is sketchy
   
//    Reversible<index_t> trailPointer;
    
//    size_t varPointer{0};
    size_t litPointer{1};
//    int backtracks{0};
    
    bool weight_reasons{false};
};

struct LearningRateHeapFactory : MakeVariableHeuristicFactory<LearningRateHeapFactory> {
    LearningRateHeapFactory();

    template<concepts::scalar T>
    auto build_impl(Solver<T>& solver) const -> VariableHeuristic<T> {
        return std::make_unique<LearningRateHeap<T>>(solver);
    }
};

}

#endif // TEMPO_LearningRateHEAP_HPP
