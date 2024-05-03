/**
 * @author Tim Luchterhand
 * @date 27.04.23.
 */

#ifndef SCHEDCL_TESTING_HPP
#define SCHEDCL_TESTING_HPP

#include <filesystem>
#include <random>
#include <concepts>

#include "util/Matrix.hpp"
#include "util/serialization.hpp"
#include "Global.hpp"

namespace tempo::testing {
    struct TestData {
        static constexpr auto TestProblem = __TEST_DATA_DIR__ "/test_problem.json";
        static constexpr auto ExtendedTestProblem = __TEST_DATA_DIR__ "/extended_test_problem.json";
        static constexpr auto GraphBuilderConfig = __TEST_DATA_DIR__ "/graph_builder_config.json";
        static constexpr auto TestNN = __TEST_DATA_DIR__ "/Identity_export.pt";
    };

    template<typename T>
    void setTaskDistance(task from, task to, T dist, Matrix<T> &distMat) {
        using namespace tempo;
        distMat.at(START(from), END(to)) = dist;
    }

    template<typename T>
    void setTaskDurations(task t, T minDuration, T maxDuration, T earliestStart, T latestCompletion, Matrix<T> &distMat) {
        using namespace tempo;
        distMat.at(START(t), END(t)) = maxDuration;
        distMat.at(END(t), START(t)) = -minDuration;
        distMat.at(START(t), ORIGIN) = earliestStart;
        distMat.at(HORIZON, END(t)) = latestCompletion;
    }

    template<typename T>
    void setUpperBound(T bound, tempo::Matrix<T> &distMat) {
        using namespace tempo;
        distMat.at(ORIGIN, HORIZON) = bound;
    }

    template<typename T>
    void setLowerBound(T bound, tempo::Matrix<T> &distMat) {
        using namespace tempo;
        distMat.at(HORIZON, ORIGIN) = -bound;
    }

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
}

#endif //SCHEDCL_TESTING_HPP
