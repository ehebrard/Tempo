/**
* @author Tim Luchterhand
* @date 14.11.24
* @brief contains random generators and utilities
*/

#ifndef TEMPO_RANDOM_HPP
#define TEMPO_RANDOM_HPP

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
    bool randomEventOccurred(double probability) noexcept {
        return random() % Res < static_cast<unsigned long>(probability * Res);
    }
}

#endif //TEMPO_RANDOM_HPP
