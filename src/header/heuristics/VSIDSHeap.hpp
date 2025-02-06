/**
* @author Emmanuel Hebrard
* @date 12.09.24
* @brief
*/

#ifndef TEMPO_VSIDSHEAP_HPP
#define TEMPO_VSIDSHEAP_HPP


//#include <algorithm>

#include "util/traits.hpp"
#include "util/SparseSet.hpp"
#include "util/Heap.hpp"

#include "heuristic_interface.hpp"
#include "ReversibleObject.hpp"

namespace tempo {
    template<typename T>
    class Solver;
}

namespace tempo::heuristics {
    /**
     * @brief Random variable selection strategy
     * @details @copybrief Randomly chooses a variable from the remenaing search variables
     * @note Right now, only binary variables are selected
     */
    struct VSIDSHeap {

        template<concepts::scalar T>
        VSIDSHeap(Solver<T>& solver) :
        handlerToken(solver.ClauseAdded.subscribe_handled(
                                                                             [this](const auto &arg) { this->updateActivity(arg); }))
        , decay(solver.getOptions().vsids_decay)
        {
            
            var_heap.resize(solver.boolean.size());
            std::iota(var_heap.begin(), var_heap.end(), 0);
            var_activity.resize(solver.boolean.size(), base_epsilon);
            index.resize(solver.boolean.size(),0);

            trail.push_back(var_heap.size());
            
            //@TODO: better initialisation
            for(size_t i{0}; i<var_heap.size(); ++i) {
                std::swap(var_heap[i], var_heap[i+random()%(var_heap.size()-i)]);
                index[var_heap[i]] = static_cast<index_t>(i);
            }
        }
        
        VSIDSHeap(const VSIDSHeap&) = delete;
        VSIDSHeap(VSIDSHeap &&) = delete;
        VSIDSHeap &operator=(const VSIDSHeap&) = delete;
        VSIDSHeap &operator=(VSIDSHeap &&) = delete;
        ~VSIDSHeap() = default;
        
        /**
         * Heuristic interface
         * @tparam T timing type
         * @param solver solver for which to select the variable
         * @return randomly selected variable
         */
        template<concepts::scalar T>
        auto nextVariable(const Solver<T> &solver) noexcept -> VariableSelection {
            const concepts::same_template<SparseSet> auto &variables = solver.getBranch();
            
            assert(not variables.empty());
            assert(checkHeap(0));
            
            auto n{trail.back()};
            while(static_cast<int>(trail.size()) > solver.level()) {
                trail.pop_back();
            }
            
            while(n < trail.back()) {
                                heap::percolate_up(var_heap.begin(), n, index, [&](const var_t x, const var_t y) {return var_activity[x] > var_activity[y];});
                ++n;
            }
    
            auto last{trail.back()};
            var_t x;
            do {
                x = pickBest(last);
            } while(not variables.has(x));
            
            trail.push_back(last);
            
//            for(auto y : variables) {
//                if(var_activity[y] > var_activity[x]) {
//                    std::cout << "bug, b" << y << " has higher activity than b" << x << ": " << var_activity[y] << " vs. " << var_activity[x] << std::endl;
//                    exit(1);
//                }
//            }
            
            assert(checkHeap(0));
            
            return {x, VariableType::Boolean};
        }
        
        template<class Iterable>
        void updateActivity(const Iterable& literals) {
            bool need_rescaling{false};
            for(auto l : literals) {
                if(not l.isNumeric()) {
                    var_activity[l.variable()] += epsilon;
                    need_rescaling |= (var_activity[l.variable()] > max_activity);
                    if(index[l.variable()] < trail.back()) {
                        heap::percolate_up(var_heap.begin(), index[l.variable()], index, [&](const var_t x, const var_t y) {return var_activity[x] > var_activity[y];});
                    }
                }
            }
            
            if(need_rescaling) {
                
//                std::cout << "rescale";
//                for(auto a : var_activity) {
//                    std::cout << " " << a;
//                }
                
                auto [lp, up] = std::ranges::minmax_element(var_activity.begin()+1, var_activity.end());
                double l{std::numeric_limits<double>::max()};
                if(lp != var_activity.end())
                    l = *lp;
                
                double u{-std::numeric_limits<double>::max()};
                if(up != var_activity.end())
                    u = *up;

                auto factor{base_gap / (u - l)};
//                for(auto &a : var_activity) {
//                    a = base_epsilon + (a - l) * factor;
//                }
                for(auto a{var_activity.begin()+1}; a!=var_activity.end(); ++a) {
                    *a = base_epsilon + (*a - l) * factor;
                }
                
//                std::cout << "\n ----> ";
//                for(auto a : var_activity) {
//                    std::cout << " " << a;
//                }
//                std::cout << std::endl;
//                std::cout << std::endl;
                
                epsilon = base_epsilon;
            } else {
                epsilon /= decay;
            }
            
//            std::cout <<" epsilon = " << epsilon << std::endl;
            
            assert(checkHeap(0));
        }
        
        bool checkHeap(const int i) {
            auto lc{heap::left(i)};
            auto rc{heap::right(i)};
            
            auto n{static_cast<int>(trail.back())};
            
            bool ok_left = ((lc >= n) or ((var_activity[var_heap[i]] >= var_activity[var_heap[lc]]) and checkHeap(lc)));
            
            if(ok_left) {
                auto ok_right = ((rc >= n) or ((var_activity[var_heap[i]] >= var_activity[var_heap[rc]]) and checkHeap(rc)));
                return ok_right;
            }
            
            return false;
        }
        
        var_t pickBest(size_t& _end_) {
            var_t x{var_heap[0]};
            heap::remove_min(var_heap.begin(), var_heap.begin()+_end_, index, [&](const var_t x, const var_t y) {return var_activity[x] > var_activity[y];});
            --_end_;
            return x;
        }
 
        // the successive sizes of the heap -> decreasing at each call to nextVariable()
        // |trail| should therefore be equal to solver.level(), if not, there has been a backtrack
        // then the values in trail can be used to re-integrate the necessary variables
        std::vector<size_t> trail;
        
        // position of each variable in the heap (so that we can percolate)
        std::vector<index_t> index;
        
        
        std::vector<var_t> var_heap;
        std::vector<double> var_activity;
        
        SubscriberHandle handlerToken;
        
        constexpr static const double base_epsilon{1e-6};
        const double decay{.999};
        constexpr static const double max_activity{1e12};
        constexpr static const double base_gap{1 - base_epsilon};
        
        double epsilon{base_epsilon};
    };
}

#endif //TEMPO_RANDOMVARIABLESELECTION_HPP
