/**
* @author Tim Luchterhand
* @date 14.11.24
* @brief contains random generators and utilities
*/

#ifndef TEMPO_RANDOM_HPP
#define TEMPO_RANDOM_HPP

#include <limits>
#include <ranges>
#include <optional>

namespace tempo {

    /**
     * set the random seed
     * @param x_ new x value
     * @param y_ new y value
     * @param x_ new z value
     */
    void seed(unsigned long x_, unsigned long y_ = 362436069, unsigned long z_ = 521288629) noexcept;

    /**
     * create a random integer
     * @return random integer number
     * @note period of the generator is period 2^96-1
     */
    unsigned long random() noexcept;

    /**
     * helper method for triggering random events
     * @tparam Res probability resolution (default 4 digits)
     * @param probability probability of the event
     * @return true if event occurred, false otherwise
     */
    template<unsigned long Res = 10000>
    bool random_event_occurred(double probability) noexcept {
        return random() % Res < static_cast<unsigned long>(probability * Res);
    }

    /**
     * @brief Random number generator that can be used in stl algorithms like std::shuffle
     */
    struct RNG {
        using result_type = decltype(random());

        static constexpr auto min() { return std::numeric_limits<result_type>::min(); }
        static constexpr auto max() { return std::numeric_limits<result_type>::max(); }
        auto operator()() const noexcept -> unsigned long;
    };

    /**
     * Selects a random element from a range
     * @tparam R range type
     * @param range the range
     * @param maxIdx optional maximum index (exclusive). UB if greater than size of range!
     * @return random element from the range
     */
    template<std::ranges::random_access_range R> requires(std::ranges::sized_range<R>)
    decltype(auto) random_select(R &&range, std::optional<std::size_t> maxIdx = {}) {
        auto max = maxIdx.value_or(std::ranges::size(range));
        return std::forward<R>(range)[random() % max];
    }

    /**
     * selects a random iterator to an element in a range
     * @tparam R range type
     * @param range the range
     * @param maxIdx optional maximum index (exclusive). UB if greater than size of range!
     * @return random iterator to an element in the range
     */
    template<std::ranges::random_access_range R> requires(std::ranges::sized_range<R>)
    auto random_select_iterator(R &&range, std::optional<std::size_t> maxIdx = {}) {
        auto idx = random() % maxIdx.value_or(std::ranges::size(range));
        return std::forward<R>(range).begin() + idx;
    }
}

#endif //TEMPO_RANDOM_HPP
