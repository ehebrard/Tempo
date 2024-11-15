/**
* @author Tim Luchterhand
* @date 15.11.24
* @brief
*/

#include <gtest/gtest.h>
#include <vector>
#include <unordered_set>

#include "util/random.hpp"
#include "testing.hpp"

struct TestSampler : public tempo::DistributionSampler {
    using tempo::DistributionSampler::DistributionSampler;
    using tempo::DistributionSampler::getCDF;
};

TEST(random, sampler_ctor_cdf) {
    std::vector cdf{1ul, 2ul, 3ul};
    TestSampler sampler(cdf);
    EXPECT_EQ(sampler.getCDF(), cdf);
}

TEST(random, sampler_ctor_pdf) {
    std::vector pdf{1, 2, 4, 1};
    TestSampler sampler(pdf);
    EXPECT_EQ(sampler.getCDF(), (std::vector<unsigned long>{1, 3, 7, 8}));
    sampler = TestSampler(std::vector{0.5, 0.75, 0.25}, 100);
    EXPECT_EQ(sampler.getCDF(), (std::vector<unsigned long>{50, 125, 150}));
    EXPECT_THROW(TestSampler(std::vector{1, 4, -1, 2, 9}), std::runtime_error);
}

template<tempo::detail::integer_range PDF, typename T, typename DrawFun>
void testRandomDraw(const PDF &pdf, const std::unordered_set<T> elements, const DrawFun &drawFun) {
    using namespace tempo;
    for (int i = 0; i < 1000; ++i) {
        auto elems = elements;
        seed(tempo::testing::random_int(0, 10000));
        ReplacementDistributionSampler sampler(pdf);
        while (not sampler.exhausted()) {
            auto elem = drawFun(sampler);
            ASSERT_TRUE(elems.contains(elem));
            elems.erase(elem);
        }

        EXPECT_TRUE(elems.empty());
    }
}

TEST(random, replacement_sampler) {
    using namespace tempo;
    testRandomDraw(std::vector{1, 2, 4, 1}, std::unordered_set{0, 1, 2, 3},
                   [](auto &sampler) { return sampler.random(); });
}

TEST(random, replacement_sampler_from_range) {
    using namespace tempo;
    testRandomDraw(std::vector{1, 2, 4, 1}, std::unordered_set{9, 2, 5, 1},
                   [elems = {9, 2, 5, 1}](auto &sampler) { return sampler.randomSelect(elems); });
}