/**
* @author Tim Luchterhand
* @date 08.11.24
* @brief
*/

#ifndef TEMPO_SOLUTIONPOOL_HPP
#define TEMPO_SOLUTIONPOOL_HPP

#include <vector>
#include <ranges>
#include <stdexcept>

#include "Global.hpp"
#include "Constant.hpp"
#include "util/traits.hpp"
#include "util/random.hpp"
#include "Literal.hpp"
#include "Solver.hpp"
#include "Solution.hpp"

namespace tempo::lns {

    /**
     * @brief Basically a stack that also supports removal from random positions
     * @tparam T element type
     */
    template<typename T>
    struct Pool : protected std::vector<T> {
        using std::vector<T>::size;
        using std::vector<T>::empty;
        using std::vector<T>::emplace_back;
        using std::vector<T>::clear;
        using std::vector<T>::begin;
        using std::vector<T>::end;
        using std::vector<T>::cbegin;
        using std::vector<T>::cend;

        const T &peekLast() const noexcept {
            return this->back();
        }

        T &peekLast() noexcept {
            return this->back();
        }

        T popLast() noexcept(std::is_nothrow_move_constructible_v<T>) {
            auto elem = std::move(this->back());
            this->pop_back();
            return elem;
        }

        T &peekRandom() noexcept {
            return random_select(*this);
        }

        const T&peekRandom() const noexcept {
            return const_cast<Pool*>(this)->peekRandom();
        }

        T popRandom() {
            auto it = random_select_iterator(*this);
            auto elem = std::move(*it);
            this->erase(it);
            return elem;
        }
    };

    /**
     * @brief Holds solutions (literal polarities)
     * @tparam T problem objective type
     */
    template<concepts::scalar T>
    class SolutionPool : public Pool<Solution<T>> {
        NumericVar<T> makespan;
        T makespanValue = Constant::Infinity<T>;
    public:
        explicit SolutionPool(const NumericVar<T> &makespanVar) noexcept : makespan(makespanVar) {}

        void addSolution(const Solver<T> &solver) {
            if (not solver.boolean.hasSolution() or not solver.numeric.hasSolution()) {
                throw std::runtime_error("solver hasn't found a solution yet");
            }

            this->emplace_back(solver);
            makespanValue = this->back().numeric.lower(makespan);
        }

        T bestMakespan() const noexcept {
            return makespanValue;
        }

        T getMakespan() const noexcept {
            return this->back().numeric.lower(makespan);
        }

        T getMakespan(const Solution<T> &sol) const noexcept {
            return sol.numeric.lower(makespan);
        }
    };

}

#endif //TEMPO_SOLUTIONPOOL_HPP
