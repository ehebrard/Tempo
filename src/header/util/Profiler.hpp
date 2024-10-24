/**
* @author Tim Luchterhand
* @date 09.08.2024
* @brief Timing and profiling utilities
*/

#ifndef TEMPO_PROFILER_HPP
#define TEMPO_PROFILER_HPP

#include <string>
#include <cmath>
#include <chrono>
#include <unordered_map>
#include <vector>
#include <ostream>
#include <iomanip>
#include <algorithm>
#include <bit>
#include <ranges>
#include <Iterators.hpp>

#include "heuristics/heuristic_interface.hpp"

namespace tempo::util {

    namespace detail {
        template<typename T>
        struct timing_symbol {};

        template<>
        struct timing_symbol<std::chrono::milliseconds> {
            static constexpr auto symbol = "ms";
        };

        template<>
        struct timing_symbol<std::chrono::microseconds> {
            static constexpr auto symbol = "µs";
        };

        template<>
        struct timing_symbol<std::chrono::seconds> {
            static constexpr auto symbol = "s";
        };

        using TP = std::chrono::time_point<std::chrono::high_resolution_clock>;

        template<std::integral T>
        int printLen(T val) {
            return static_cast<int>(std::log(val + 1) / std::log(10)) + 1;
        }
    }

    /**
     * @brief Represents an event with a duration.
     * @details @copybrief
     */
    struct TimingEvent {
        detail::TP start; ///< start point of the event
        detail::TP end; ///< end point of the event

        /**
         * Ctor
         * @param start start point of the event
         * @param end end point of the event
         */
        TimingEvent(const detail::TP &start, const detail::TP &end) noexcept;

        /**
         * gets the duration of the event
         * @tparam T duration type
         * @return duration of the event
         */
        template<typename T>
        auto duration() const noexcept {
            return std::chrono::duration_cast<T>(end - start).count();
        }
    };

    /**
     * @brief Profiling result
     * @details @copybrief
     * @tparam T duration type
     */
    template<typename T>
    struct Result {
        T min, max, avg, stddev, med, sum;
    };

    /**
     * @brief Profiler that manages multiple events.
     * @details @copybrief
     */
    class Profiler {
        std::unordered_map<std::string, std::vector<TimingEvent>> events;

    public:

        /**
         * Adds an event ot the profiler
         * @param event event to be added
         * @param name event name
         */
        void addEvent(const TimingEvent &event, const std::string &name);

        /**
         * Adds an event ot the profiler
         * @param start start point of the even
         * @param end end point of the event
         * @param name event name
         */
        void addEvent(detail::TP start, detail::TP end, const std::string &name);

        /**
         * gets the profiling result for an event
         * @tparam T duration type
         * @param eventName name of the event
         * @return profiling result
         */
        template<typename T>
        auto getResult(const std::string &eventName) const {
            using Res = decltype(std::declval<TimingEvent>().duration<T>());
            Res sum = 0;
            Res sqSum = 0;
            Res max = std::numeric_limits<Res>::min();
            Res min = std::numeric_limits<Res>::max();
            const auto &eventList = events.at(eventName);
            std::vector<Res> results;
            results.reserve(eventList.size());
            for (const auto &evt : eventList) {
                auto val = evt.duration<T>();
                sum += val;
                results.emplace_back(val);
                sqSum += val * val;
                max = std::max(max, val);
                min = std::min(min, val);
            }

            Res mean = sum / eventList.size();
            std::ranges::sort(results);
            Res med;
            if (results.empty()) {
                med = 0;
            } else if (results.size() % 2 == 0) {
                med = (results[results.size() / 2] + results[results.size() / 2 - 1]) / 2;
            } else {
                med = results[results.size() / 2];
            }

            auto stddev = static_cast<Res>(std::sqrt(sqSum / eventList.size() - mean * mean));
            return Result<Res>{.min=min, .max=max, .avg=mean, .stddev=stddev, .med = med, .sum=sum};

        }

        /**
         * prints a profiling result to an out stream
         * @tparam T timing type
         * @param eventName name of the event
         * @param os out stream
         * @param nameWidth optional width formatting parameter for the event name
         * @param valWidth optional width formatting parameter for the result values
         */
        template<typename T>
        void print(const std::string &eventName, std::ostream &os, int nameWidth = 0, int valWidth = 0) const {
            auto [min, max, avg, stddev, med, sum] = getResult<T>(eventName);
            constexpr auto s = detail::timing_symbol<T>::symbol;
            os << "-- " << std::setw(nameWidth) << std::left << eventName << ": \tmin: " << std::setw(valWidth)
               << std::left << min << s
               << ", max: " << std::setw(valWidth) << std::left << max << s << ", avg : " << std::setw(valWidth)
               << std::left << avg << s
               << ", std: " << std::setw(valWidth) << std::left << stddev << s << ", median: " << std::setw(valWidth)
               << std::left << med
               << s << ", total: " << sum << s << "\n";
        }

        /**
         * Checks whether profiling data on a specific event is available
         * @param event event name
         * @return true if event occurred at least once, false otherwise
         */
        bool has(const std::string &event) const noexcept {
            return events.contains(event);
        }

        /**
         * prints all events to an out stream
         * @tparam T timing type
         * @param os out stream
         */
        template<typename T>
        void printAll(std::ostream &os) const {
            using ResT = decltype(getResult<T>("").min);
            constexpr auto NumFields = sizeof(Result<ResT>) / sizeof(T);
            constexpr auto s = detail::timing_symbol<T>::symbol;
            int nameWidth = 0;
            std::array<int, NumFields> valWidths{0};
            std::vector<Result<ResT>> results;
            for (const auto &[name, _] : events) {
                nameWidth = std::max(nameWidth, static_cast<int>(name.length()) + 1);
                results.emplace_back(getResult<T>(name));
                const auto fields = std::bit_cast<std::array<ResT, NumFields>>(results.back());
                for (auto [f, fWidth] : iterators::zip(fields, valWidths)) {
                    fWidth = std::max(fWidth, static_cast<int>(detail::printLen(f) + 1));
                }
            }

            for (auto [name, res] : iterators::zip(events | std::views::keys, results)) {
                os << "-- " << std::setw(nameWidth)
                   << std::left << name << ": \tmin: " << std::setw(valWidths[0]) << std::left << res.min << s
                   << ", max: " << std::setw(valWidths[1]) << std::left << res.max << s
                   << ", avg : " << std::setw(valWidths[2]) << std::left << res.avg << s
                   << ", std: " << std::setw(valWidths[3]) << std::left << res.stddev << s
                   << ", median: " << std::setw(valWidths[4]) << std::left << res.med << s
                   << ", total: " << std::setw(valWidths[5]) << std::left << res.sum << s << "\n";
            }
        }
    };

    /**
     * @brief Used to measure time between start a start and a stop event.
     * @details @copybrief
     */
    class StopWatch {
    protected:
        detail::TP startTp;
    public:

        /**
         * Ctor. Starts the timer
         */
        StopWatch();

        /**
         * Resets the start event to the current time
         */
        void start();

        /**
         * Get the timing event since the last call to start and the time of the call
         * @return TimingEvent with start point set to the last call to start() and end
         * event the current time
         */
        TimingEvent getTiming() const;
    };


    /**
     * @brief Stop watch that automatically adds a timing event to a profiler at destruction
     */
    class ScopeWatch: protected StopWatch {
        Profiler &profiler;
        std::string name;
    public:
        /**
         * CTor. Starts the stop watch
         * @param profiler profiler where the event is registered
         * @param eventName name of the event
         */
        ScopeWatch(Profiler &profiler, std::string eventName);

        /**
         * DTor. Stops the stop watch and adds the event to the scheduler
         */
        ~ScopeWatch();
    };

    /**
     * @brief Wrapper class that measures the evaluation times of heuristics.
     * @details @copybrief.
     * It satisfies the heuristic interfaces and forwards all calls to the actual heuristic while measuring their
     * execution time.
     * @tparam H heuristic type
     * @tparam Resolution time resolution
     */
    template<typename H>
    class HeuristicProfiler {
        H h;
        Profiler &profiler;
    public:
        /**
         * Ctor
         * @tparam Args argument types for the actual heuristic
         * @param profiler profiler to use
         */
        template<typename ...Args>
        constexpr explicit HeuristicProfiler(Profiler &profiler, Args &&...args) :
                h(std::forward<Args>(args)...), profiler(profiler) {}

        /**
         * Variable heuristic interface
         * @tparam T timing type
         * @param solver
         * @return variable selection
         */
        template<concepts::scalar T> requires(heuristics::variable_heuristic<H, T>)
        auto nextVariable(const Solver<T> &solver) -> heuristics::VariableSelection {
            ScopeWatch sw(profiler, "variable branching");
            return h.nextVariable(solver);
        }

        /**
         * Value heuristic interface
         * @tparam T timing type
         * @param x variable selection
         * @param solver
         * @return branching literal
         */
        template<concepts::scalar T> requires(heuristics::value_heuristic<H, Solver<T>>)
        auto valueDecision(heuristics::VariableSelection x, const Solver<T> &solver) -> Literal<T> {
            ScopeWatch sw(profiler, "value branching");
            return h.valueDecision(x, solver);
        }

        /**
         * Heuristic interface
         * @tparam T timing type
         * @param solver
         * @return branching literal
         */
        template<concepts::scalar T> requires(heuristics::heuristic<H, T>)
        auto branch(const Solver<T> &solver) -> Literal<T> {
            ScopeWatch sw(profiler, "branching heuristic");
            return h.branch(solver);
        }

        /**
         * Access to the actual heuristic
         * @return pointer to the stored heuristic
         */
        constexpr auto operator->() const noexcept {
            return &h;
        }
    };

    /**
     * Movable wrapper around HeuristicProfiler
     * @tparam H heuristic type
     * @tparam Resolution time resolution
     */
    template<typename H>
    using ProfiledHeuristic = heuristics::MovableHeuristic<HeuristicProfiler<H>>;

}

#endif //TEMPO_PROFILER_HPP
