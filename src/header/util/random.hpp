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
#include <vector>
#include <concepts>
#include <numeric>
#include <algorithm>
#include <stdexcept>

#include "util/traits.hpp"

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
     * selects a random iterator to an element in a range
     * @tparam R range type
     * @param range the range
     * @param maxIdx optional maximum index (exclusive). UB if greater than size of range!
     * @return random iterator to an element in the range
     */
    template<std::ranges::random_access_range R> requires(std::ranges::sized_range<R>)
    auto random_select_iterator(R &&range, std::optional<std::size_t> maxIdx = {}) {
        auto idx = random() % maxIdx.value_or(std::ranges::size(range));
        return std::ranges::begin(std::forward<R>(range)) + idx;
    }


    /**
     * Selects a random element from a range
     * @tparam R range type
     * @param range the range
     * @param maxIdx optional maximum index (exclusive). UB if greater than size of range!
     * @return random element from the range
     */
    template<std::ranges::random_access_range R> requires(std::ranges::sized_range<R>)
    decltype(auto) random_select(R &&range, std::optional<std::size_t> maxIdx = {}) {
        return *random_select_iterator(std::forward<R>(range), maxIdx);
    }


    namespace detail {
        template<typename R>
        concept scalar_range = std::ranges::range<R> and concepts::scalar<std::ranges::range_value_t<R>>;
    }

    class ReplacementDistributionSampler ;

    /**
     * @brief Samples from an arbitrary distribution function using inverse transform sampling.
     * @details @copybrief
     * This sampler is defined over an arbitrary PDF. It works the following way. Given the following PDF
     *
     * ^
     * |            |
     * |    |       |     |
     * |    | |     | |   |
     * |  | | | | | | | | | |
     * |--0-1-2-3-4-5-6-7-8-9--> x
     * the sampler selects values on the x-axis. In this example, the x-value 5 is the most likely
     */
    class DistributionSampler {
        std::vector<unsigned long> cdf;
        friend class ReplacementDistributionSampler;
        static constexpr auto FPTolerance = 1.05;

        static bool conversionNeeded(auto val) noexcept {
            if constexpr (std::integral<decltype(val)>) {
                return false;
            } else {
                return val / static_cast<double>(static_cast<long>(val)) > FPTolerance;
            }
        }

        template<concepts::scalar T>
        static double smoothVal(T val, double mean, double lambda) noexcept {
            return static_cast<double>((1 - lambda) * val + lambda * mean);
        }

        template<detail::scalar_range R>
        void accumulate(const R&values, bool fpConversion, unsigned long fpPrecision) {
            cdf.reserve(std::ranges::size(values));
            unsigned long sum = 0;
            auto converted = values | std::ranges::views::transform([fpConversion, fpPrecision](auto val) {
                if constexpr (std::integral<decltype(val)>) {
                    return val;
                } else {
                    if (fpConversion) {
                        return val * fpPrecision;
                    }

                    return val;
                }
            });

            for (auto val : converted) {
                if (val < 0) {
                    throw std::invalid_argument("negative probability mass");
                }

                sum += static_cast<unsigned long>(val);
                cdf.emplace_back(sum);
            }
        }

    protected:
        [[nodiscard]] auto getCDF() const noexcept -> const std::vector<unsigned long>&;
    public:
        /**
         * CTor
         * @param cdf Cumulative density function. Does not need to sum up to one but needs to be a valid CDF
         * otherwise, i.e. ordered in strictly ascending order. If not, this may lead to UB!
         */
        explicit DistributionSampler(std::vector<unsigned long> cdf);

        /**
         * Ctor
         * @tparam R pdf range type
         * @param pdf Probability density function. Does not need to sum to 1 but must not contain negative values
         * @param smoothingFactor optional smoothing factor that brings pdf values closer to their mean value. 0
         * corresponds to no smoothing, 1 creates a uniform distribution
         * @param fpPrecision fixed point precision for float to integer conversion
         */
        template<detail::scalar_range R> requires(std::ranges::sized_range<R>)
        explicit DistributionSampler(const R &pdf, double smoothingFactor = 0, unsigned long fpPrecision = 10000) {
            using T = std::ranges::range_value_t<R>;
            if (std::ranges::empty(pdf)) {
                throw std::invalid_argument("empty pdf");
            }

            if (smoothingFactor < 0 or smoothingFactor > 1) {
                throw std::invalid_argument("Invalid smoothing factor");
            }

            if (smoothingFactor == 1) {
                auto values = std::views::iota(1ul, std::ranges::size(pdf) + 1);
                cdf.assign(values.begin(), values.end());
                return;
            }

            if (smoothingFactor == 0) {
                bool convert = std::floating_point<T> and std::ranges::any_of(pdf, [](auto val) {
                    return conversionNeeded(val);
                });
                accumulate(pdf, convert, fpPrecision);
                return;
            }

            bool fpConversionNeeded = false;
            std::vector<double> smoothed;
            smoothed.reserve(std::ranges::size(pdf));
            const auto mean = std::accumulate(std::ranges::begin(pdf), std::ranges::end(pdf), T(0)) /
                              static_cast<double>(std::ranges::size(pdf));
            for (auto v : pdf) {
                smoothed.emplace_back(smoothVal(v, mean, smoothingFactor));
                fpConversionNeeded |= conversionNeeded(smoothed.back());
            }

            accumulate(smoothed, fpConversionNeeded, fpPrecision);
        }

        /**
         * Gets the x-index selected according to the probability distribution
         * @return
         */
        [[nodiscard]] std::size_t random() const noexcept;

        /**
         * Selects a random iterator to an element in the given range according to the PDF
         * @tparam R range type
         * @param range the range
         * @return iterator to element in the range
         * @note The range must contain as many elements as the PDF
         */
        template<std::ranges::random_access_range R> requires(std::ranges::sized_range<R>)
        auto randomSelectIterator(R &&range) {
            if (std::ranges::size(std::forward<R>(range)) != cdf.size()) {
                throw std::runtime_error("range needs to have as many elements as PDF");
            }

            return std::ranges::begin(std::forward<R>(range)) + this->random();
        }

        /**
         * Selects a random element from the given range according to the PDF
         * @tparam R range type
         * @param range the range
         * @return element in the range
         * @note The range must contain as many elements as the PDF
         */
        template<std::ranges::random_access_range R> requires(std::ranges::sized_range<R>)
        decltype(auto) randomSelect(R &&range) const {
            return *this->randomSelectIterator(range);
        }

    };

    /**
     * @brief Samples WITH REPLACEMENT from an arbitrary distribution function using inverse transform sampling.
     * @details @copybrief
     * This class behaves similarly to DistributionSampler with the exceptions that samples are drawn WITH
     * REPLACEMENT. This means that the distribution changes with each sample drawn
     */
    class ReplacementDistributionSampler {
        std::vector<std::size_t> indices;
        DistributionSampler sampler;
    public:
        /**
         * Ctor. Same signatures as DistributionSampler
         * @tparam Args argument types
         * @param args arguments to DistributionSampler
         */
        template<typename ...Args>
        explicit ReplacementDistributionSampler(Args &&...args): sampler(std::forward<Args>(args)...) {
            std::ranges::iota_view values(static_cast<std::size_t>(0), sampler.cdf.size());
            indices.assign(values.begin(), values.end());
        }

        /**
         * Whether all values have been sampled
         * @return true if all values have been drawn, false otherwise
         */
        [[nodiscard]] bool exhausted() const noexcept;

        /**
         * Draws an x-value according to the probability distribution. The same x value cannot be drawn twice.
         * @return selected x-value
         */
        std::size_t random();

        /**
         * Selects a random iterator to an element in the given range according to the PDF. Each iterator can only be
         * drawn once
         * @tparam R range type
         * @param range the range
         * @return iterator to element in the range
         * @note The range must contain at least as many elements as the PDF. You should not remove elements from your
         * range. The sampler keeps track of the indices that have already been drawn
         */
        template<std::ranges::random_access_range R> requires(std::ranges::sized_range<R>)
        auto randomSelectIterator(R &&range) {
            if (exhausted()) {
                throw std::runtime_error("sampler is exhausted");
            }

            if (std::ranges::size(std::forward<R>(range)) < indices.back()) {
                throw std::runtime_error("range has too few elements");
            }

            return std::ranges::begin(std::forward<R>(range)) + this->random();
        }

        /**
         * Selects a random element from the given range according to the PDF. Each element can only be chosen once.
         * @tparam R range type
         * @param range the range
         * @return element in the range
         * @note The range must contain at least as many elements as the PDF. You should not remove elements from your
         * range. The sampler keeps track of the indices that have already been drawn
         */
        template<std::ranges::random_access_range R> requires(std::ranges::sized_range<R>)
        decltype(auto) randomSelect(R &&range) {
            return *this->randomSelectIterator(range);
        }
    };
}

#endif //TEMPO_RANDOM_HPP
