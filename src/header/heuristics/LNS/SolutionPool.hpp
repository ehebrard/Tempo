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
#include <optional>
#include <limits>

#include "Constant.hpp"
#include "util/traits.hpp"
#include "util/random.hpp"
#include "Solver.hpp"
#include "Solution.hpp"

namespace tempo::lns {

    /**
     * @brief Ringbuffer that supports stack operations as well as random element retrival
     * @tparam T element type
     * @note Element order is LiFo
     */
    template<typename T>
    class Pool {
        std::vector<T> buffer;
        std::size_t head = 0;
        std::size_t currSize = 0;
        std::optional<std::size_t> maxSize;

        void addElem() noexcept {
            head = currSize == 0 ? 0 : (head + 1) % capacity();
            currSize = std::min(currSize + 1, capacity());
        }

        void removeElem() noexcept {
            --currSize;
            if (head == 0) {
                head = capacity() - 1;
            } else {
                --head;
            }
        }
    public:
        /**
         * Ctor
         * @param size optional max size parameter. If not given, pool has no size limit
         */
        constexpr explicit Pool(std::optional<std::size_t> size = {}) : buffer(size.value_or(0)), maxSize(size)  {}

        /**
         * Maximum number of elements the pool can hold
         * @return max number of elements
         */
        [[nodiscard]] constexpr std::size_t capacity() const noexcept {
            return maxSize.value_or(std::numeric_limits<std::size_t>::max());
        }

        /**
         * Number of elements contained in the pool
         * @return number of elements contained in the pool
         */
        [[nodiscard]] constexpr std::size_t size() const noexcept {
            return currSize;
        }

        /**
         * Whether pool is empty
         * @return true if pool is empty, false otherwise
         */
        [[nodiscard]] constexpr bool empty() const noexcept {
            return currSize == 0;
        }

        /**
         * @brief Add an element to the back of the pool
         * @tparam Args element ctor argument types
         * @param args arguments to element ctor
         */
        template<typename... Args>
        constexpr void emplace_back(Args &&... args) {
            addElem();
            if (currSize > buffer.size()) {
                buffer.emplace_back(std::forward<Args>(args)...);
            } else {
                std::construct_at(buffer.data() + head, std::forward<Args>(args)...);
            }
        }

        /**
         * @brief Retrieve last element without removing it
         * @return reference to last element
         */
        constexpr const T &peekLast() const noexcept {
            return buffer[head];
        }

        /**
         * @copydoc peekLast
         */
        constexpr T &peekLast() noexcept {
            return buffer[head];
        }

        /**
         * @brief Remove last element and retrieve it
         * @return element retrieved from the back of the pool
         */
        constexpr T popLast() noexcept(std::is_nothrow_move_constructible_v<T>) {
            auto elem = std::move(buffer[head]);
            removeElem();
            return elem;
        }

        /**
         * Retrieve a random element without removing it
         * @return reference to random element
         */
        constexpr T &peekRandom() noexcept {
            return random_select(buffer, currSize);
        }

        /**
         * @copydoc peekRandom
         */
        constexpr const T&peekRandom() const noexcept {
            return const_cast<Pool*>(this)->peekRandom();
        }

        /**
         * @brief Remove random element and retrieve it
         * @return element retrieved from the back of the pool
         * @note complexity is linear in the size of the pool
         */
        constexpr T popRandom() {
            auto it = random_select_iterator(buffer, currSize);
            auto elem = std::move(*it);
            buffer.erase(it);
            removeElem();
            return elem;
        }

        /**
         * remove all elements
         */
        constexpr void clear() noexcept {
            currSize = 0;
            head = 0;
            buffer.clear();
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
        explicit SolutionPool(const NumericVar<T> &makespanVar, std::optional<std::size_t> maxCapacity = {}) noexcept
            : Pool<Solution<T>>(maxCapacity), makespan(makespanVar) {}

        void addSolution(const Solver<T> &solver) {
            if (not solver.boolean.hasSolution() or not solver.numeric.hasSolution()) {
                throw std::runtime_error("solver hasn't found a solution yet");
            }

            this->emplace_back(solver);
            makespanValue = this->peekLast().numeric.lower(makespan);
        }

        T bestMakespan() const noexcept {
            return makespanValue;
        }

        T getMakespan() const noexcept {
            return this->peekLast().numeric.lower(makespan);
        }

        T getMakespan(const Solution<T> &sol) const noexcept {
            return sol.numeric.lower(makespan);
        }
    };

}

#endif //TEMPO_SOLUTIONPOOL_HPP
