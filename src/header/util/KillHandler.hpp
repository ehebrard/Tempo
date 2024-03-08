//
// Created by tim on 19.09.23.
//

#ifndef TEMPO_SIGNALHANDLER_HPP
#define TEMPO_SIGNALHANDLER_HPP
#include <atomic>

namespace tempo {
    /**
     * Singleton class that intercepts the signals SIGINT and SIGTERM if used
     */
    class KillHandler {
    public:

        /**
         * Singleton instance access
         * @return reference to the instance
         */
        static KillHandler& instance() noexcept;

        /**
         * Get the status of the Handler, i.e. the signal number of the last signal intercepted since call to reset()
         * @return
         */
        int getStatus() const noexcept;

        /**
         * Whether a signal has been received since the last call to reset()
         * @return true if signal has been intercepted, false otherwise
         */
        bool signalReceived() const noexcept;

        /**
         * Resets the intercepted signal status
         */
        void reset() noexcept;

        /**
         * Sets the internal signal flag to SIGTERM to signalize termination request
         */
        void kill() noexcept;
    private:
        static void handler(int signal) noexcept;
        static void configureHandler(int signal);
        KillHandler();
        static inline std::atomic_int status;
    };
}


#endif //TEMPO_SIGNALHANDLER_HPP
