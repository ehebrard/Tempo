/**
* @author Tim Luchterhand
* @date 14.11.24
* @brief
*/

//#include <iostream>
#include <algorithm>
#include <cassert>
#include <numeric>
#include <Iterators.hpp>

#include "util/random.hpp"

namespace tempo {

    static unsigned long x = 123456789, y = 362436069, z = 521288629;

    void seed(unsigned long x_, unsigned long y_, unsigned long z_) noexcept {
        x = x_;
        y = y_;
        z = z_;
    }

    unsigned long random() noexcept {
        unsigned long t;
        x ^= x << 16;
        x ^= x >> 5;
        x ^= x << 1;

        t = x;
        x = y;
        y = z;
        z = t ^ x ^ y;
        
//        std::cout << z << std::endl;

        return z;
    }

    auto RNG::operator()() const noexcept -> unsigned long { return random(); }

    UniformReplacementSampler::UniformReplacementSampler(std::size_t populationSize) {
        std::ranges::iota_view iota(0ul, populationSize);
        indices.reserve(populationSize);
        std::ranges::copy(iota, std::back_inserter(indices));
        std::ranges::shuffle(indices, RNG{});
    }

    std::size_t UniformReplacementSampler::random() noexcept {
        auto ret = indices.back();
        indices.pop_back();
        return ret;
    }

    bool UniformReplacementSampler::exhausted() const noexcept {
        return indices.empty();
    }

    std::size_t DistributionSampler::random() const noexcept {
        auto rVal = tempo::random() % cdf.back();
        auto res = std::ranges::upper_bound(cdf, rVal);
        assert(res != cdf.end());
        return std::distance(cdf.begin(), res);
    }

    DistributionSampler::DistributionSampler(std::vector<unsigned long> cdf): cdf(std::move(cdf)) {
        if (this->cdf.empty()) {
            throw std::runtime_error("empty CDF");
        }
    }

    auto DistributionSampler::getCDF() const noexcept -> const std::vector<unsigned long> & {
        return cdf;
    }

    std::size_t ReplacementDistributionSampler::random() {
        using namespace std::views;
        if (exhausted()) {
            throw std::runtime_error("sample is exhausted");
        }

        const auto idx = sampler.random();
        unsigned long decrement = sampler.cdf[idx];
        if (idx > 0) {
            decrement -= sampler.cdf[idx - 1];
        }

        const auto ret = indices[idx];
        auto source =
                iterators::zip(sampler.cdf | transform([decrement](auto val) { return val - decrement; }), indices) |
                drop(idx + 1);
        auto destination = iterators::zip_i(sampler.cdf.begin() + idx, indices.begin() + idx);
        std::ranges::copy(source, destination);
        sampler.cdf.pop_back(); indices.pop_back();
        return ret;
    }

    bool ReplacementDistributionSampler::exhausted() const noexcept {
        return indices.empty();
    }
}
