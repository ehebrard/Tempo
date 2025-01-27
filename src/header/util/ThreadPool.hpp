/**
 * @author Tim Luchterhand
 * @date 09.08.21.
 * @brief Thread pool implementation
 */

#ifndef TEMPO_THREADPOOL_HPP
#define TEMPO_THREADPOOL_HPP

#include <thread>
#include <vector>
#include <deque>
#include <condition_variable>
#include <memory>
#include <functional>
#include <future>

namespace tempo {
    /**
     * @brief Simple thread pool implementation
     */
    class ThreadPool {
    public:
        using JobType = std::function<void()>;

        ThreadPool(const ThreadPool &) = delete;

        ThreadPool(ThreadPool &&) = delete;

        ThreadPool &operator=(const ThreadPool &) = delete;

        ThreadPool &operator=(ThreadPool &&) = delete;

        /**
         * Ctor
         * @param nWorkers number of threads to use
         */
        explicit ThreadPool(unsigned int nWorkers);

        /**
         * Submit a job to the threadpool
         * @tparam FUN job function type
         * @param f function to execute in the pool
         * @return future result
         */
        template<typename FUN>
        auto submit(FUN &&f) -> std::future<decltype(f())> {
            using T = decltype(f());
            std::packaged_task<T()> task(std::forward<FUN>(f));
            auto ret = task.get_future();
            std::lock_guard lock(mutex);
            jobs.emplace_back([t = std::make_shared<decltype(task)>(std::move(task))]() mutable { (*t)(); });
            cv.notify_one();
            return ret;
        }

        /**
         * Number of workers used by the thread pool
         * @return
         */
        [[nodiscard]] unsigned numWorkers() const noexcept;

        /**
         * Dtor.
         * Waits for all running workers to finish and cancels all other queued jobs
         */
        ~ThreadPool();

    private:
        void doWork();

        std::vector<std::thread> workers;
        std::mutex mutex;
        std::condition_variable cv;
        std::deque<JobType> jobs;
        unsigned int active = 0;
        bool finished = false;
    };
}

#endif //TEMPO_THREADPOOL_HPP
