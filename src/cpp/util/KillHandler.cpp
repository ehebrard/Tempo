//
// Created by tim on 19.09.23.
//

#include <csignal>
#include <stdexcept>

#include "util/KillHandler.hpp"

namespace tempo {

    KillHandler &KillHandler::instance() noexcept {
        static KillHandler s;
        return s;
    }

    void KillHandler::reset() noexcept {
        KillHandler::status = 0;
    }

    int KillHandler::getStatus() const noexcept {
        return KillHandler::status;
    }

    bool KillHandler::signalReceived() const noexcept {
        return KillHandler::status != 0;
    }

    KillHandler::KillHandler() {
        KillHandler::status = 0;
        KillHandler::configureHandler(SIGINT);
        KillHandler::configureHandler(SIGTERM);
    }

    void KillHandler::handler(int signal) noexcept {
        KillHandler::status = signal;
    }

    void KillHandler::configureHandler(int signal) {
        struct sigaction sa = {};
        sa.sa_handler = KillHandler::handler;
        if (sigaction(signal, &sa, nullptr) < 0) {
            throw std::runtime_error("unable to set kill signal handler");
        }
    }

    void KillHandler::kill() noexcept {
        KillHandler::status = SIGTERM;
    }
}
