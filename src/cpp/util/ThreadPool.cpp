/**
* @author Tim Luchterhand
* @date 27.01.25
* @file ThreadPool.cpp
* @brief
*/

#include "util/ThreadPool.hpp"

namespace tempo {
    ThreadPool::ThreadPool(unsigned int nWorkers): workers(nWorkers) {
        for (auto &w: workers) {
            w = std::thread([this]() { doWork(); });
        }
    }

    unsigned ThreadPool::numWorkers() const noexcept {
        return workers.size();
    }

    ThreadPool::~ThreadPool() { {
            std::lock_guard lock(mutex);
            finished = true;
            cv.notify_all();
        }

        for (auto &w: workers) {
            if (w.joinable()) {
                w.join();
            }
        }
    }

    void ThreadPool::doWork() {
        while (true) {
            std::unique_lock lock(mutex);
            cv.wait(lock, [this]() { return !jobs.empty() || (finished && active == 0); });
            if (jobs.empty()) {
                break;
            }

            auto job = std::move(jobs.front());
            jobs.pop_front();
            active++;
            lock.unlock();
            job();
            lock.lock();
            active--;
        }

        cv.notify_all();
    }
}
