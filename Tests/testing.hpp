/**
 * @author Tim Luchterhand
 * @date 27.04.23.
 */

#ifndef TEMPO_TESTING_HPP
#define TEMPO_TESTING_HPP

#include <filesystem>
#include <random>
#include <concepts>
#include <utility>
#include <tuple>

#include "util/Matrix.hpp"
#include "util/serialization.hpp"
#include "Global.hpp"
#include "util/SchedulingProblemHelper.hpp"
#include "Model.hpp"

namespace tempo::testing {
    struct TestData {
        static constexpr auto GraphBuilderConfig = __TEST_DATA_DIR__ "/graph_builder_config.json";
        static constexpr auto TestNN = __TEST_DATA_DIR__ "/Identity_export.pt";
    };

    struct Task : public Interval<int> {
        Task(NumericVar<int> start, NumericVar<int> end, NumericVar<int> duration) :
                Interval<int>(start, end, duration) {}
    };

    struct Resource : public std::vector<Interval<int>> {
        std::vector<int> demands;
        int capacity;
        Resource(int capacity, std::vector<Interval<int>> tasks, std::vector<int> demands);
        [[nodiscard]] int getDemand(unsigned taskId) const;
        [[nodiscard]] int resourceCapacity() const;
    };

    using ProblemInstance = tempo::SchedulingProblemHelper<int, Resource>;

    class BoundProvider {
        std::vector<int> u, l;
    public:
        BoundProvider(std::vector<int> upper, std::vector<int> lower);
        [[nodiscard]] int upper(tempo::var_t var) const;
        [[nodiscard]] int lower(tempo::var_t var) const;
    };

    struct DummyScheduler {
        tempo::testing::BoundProvider numeric;
        template<typename ...Args>
        explicit DummyScheduler(Args &&...args): numeric(std::forward<Args>(args)...) {}
    };


    auto createTestProblem() -> std::pair<ProblemInstance, DummyScheduler>;

    auto createExtendedTestProblem() -> std::pair<ProblemInstance, DummyScheduler>;

    auto createRandomProblem(std::size_t numTasks, std::size_t numResources,
                             double precedenceChance = 0.3) -> std::tuple<ProblemInstance, DummyScheduler, Matrix<int>>;

    struct TaskSpec {
        int minDur, maxDur, earliestStart{0}, latestDeadline{0};
    };

    auto createTasks(const std::vector<TaskSpec> &specs) -> std::pair<std::vector<Interval<int>>, DummyScheduler>;

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
}

#endif //TEMPO_TESTING_HPP
