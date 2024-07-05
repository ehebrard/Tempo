/**
 * @author Tim Luchterhand
 * @date 27.04.23.
 */

#ifndef SCHEDCL_TESTING_HPP
#define SCHEDCL_TESTING_HPP

#include <filesystem>
#include <random>
#include <concepts>
#include <utility>

#include "util/Matrix.hpp"
#include "util/serialization.hpp"
#include "Global.hpp"
#include "util/SchedulingProblemHelper.hpp"
#include "Model.hpp"

namespace tempo::testing {
    struct TestData {
        static constexpr auto TestProblem = __TEST_DATA_DIR__ "/test_problem.json";
        static constexpr auto ExtendedTestProblem = __TEST_DATA_DIR__ "/extended_test_problem.json";
        static constexpr auto GraphBuilderConfig = __TEST_DATA_DIR__ "/graph_builder_config.json";
        static constexpr auto TestNN = __TEST_DATA_DIR__ "/Identity_export.pt";
    };

    struct Resource : public std::vector<Interval<int>> {
        std::vector<int> demands;
        int capacity;
        Resource(int capacity, std::vector<Interval<int>> tasks, std::vector<int> demands);
        [[nodiscard]] int getDemand(unsigned taskId) const;
        [[nodiscard]] int resourceCapacity() const;
    };

    using ProblemInstance = tempo::SchedulingProblemHelper<int, Resource>;

    auto createTestProblem() -> ProblemInstance;

    auto createTasks(std::initializer_list<std::pair<int, int>> durations) -> std::vector<Interval<int>>;

    auto createDummyTasks(unsigned numberOfTasks) -> std::vector<Interval<int>>;


    /**
     * Generates a random integer value in [min, max]
     * @tparam T integral type of value
     * @param min lower bound
     * @param max upper bound
     * @return random value
     */
    template<std::integral T>
    T random_int(T min, T max) {
        static std::random_device rd;
        static std::default_random_engine el(rd());
        std::uniform_int_distribution<T> dist(min, max);
        return dist(el);
    }

    /**
     * Generates a random float value in [min, max)
     * @tparam T floating point type of value
     * @param min lower bound
     * @param max upper bound
     * @return random value
     */
    template<std::floating_point T>
    T random_float(T min, T max) {
        static std::random_device rd;
        static std::default_random_engine el(rd());
        std::uniform_real_distribution<T> dist(min, max);
        return dist(el);
    }

    struct TaskSpec {
        int minDur, maxDur, release, deadline;
    };

}

#endif //SCHEDCL_TESTING_HPP
