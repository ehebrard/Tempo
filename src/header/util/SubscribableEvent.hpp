//
// Created by tim on 16.11.22.
//

#ifndef TEMPO_SUBSCRIBABLEEVENT_HPP
#define TEMPO_SUBSCRIBABLEEVENT_HPP

#include <functional>
#include <vector>
#include <ranges>
#include <memory>
#include <algorithm>

namespace tempo {

    namespace impl {
        struct EventStatus {
            [[nodiscard]] constexpr bool isAlive() const noexcept {
                return alive;
            }

            constexpr void invalidate() noexcept {
                alive = false;
            }
        private:
            bool alive = true;
        };

        using EventStatusPtr = std::shared_ptr<EventStatus>;
        using cEventStatusPtr = std::shared_ptr<const EventStatus>;
    }

    /**
     * @brief Identifies a function handler subscribed to a SubscribableEvent. Can be used to unsubscribe the function
     * handler.
     * @details @copybrief
     * This handle automatically unsubscribes the associated event handler at destruction. It avoids undefined behaviour
     * if the event is already destroyed. However, this mechanism is not thread-safe!
     */
    class SubscriberHandle {
        template<typename ...Args>
        friend class SubscribableEvent;
    public:
        using Id = unsigned long long;
        SubscriberHandle(const SubscriberHandle &) = delete;
        SubscriberHandle &operator=(const SubscriberHandle &) = delete;
        SubscriberHandle(SubscriberHandle &&other) noexcept;
        SubscriberHandle &operator=(SubscriberHandle &&other) noexcept;

        /**
         * mark the handle as disposed. It can then no longer be used to unregister the event handler from the event.
         * Automatic unsubscription at destruction is also disabled
         */
        void dispose() noexcept;

        /**
         * Whether the handle has been manually disposed
         * @return
         */
        [[nodiscard]] bool isDisposed() const noexcept;

        /**
         * Manually unregister the associated event handler
         */
        void unregister();

        /**
         * DTor. Automatically unregisters the associated event handler if not already disposed.
         */
        ~SubscriberHandle();

        /**
         * swaps the contents of two handles
         * @param lhs left hand side
         * @param rhs right hand side
         */
        friend void swap(SubscriberHandle &lhs, SubscriberHandle &rhs) noexcept;

    protected:
        template<typename Fun>
        SubscriberHandle(Fun &&deregister, impl::cEventStatusPtr status, Id id, bool disposed):
                deregister(std::forward<Fun>(deregister)), eventStatus(std::move(status)), disposed(disposed), id(id) {}

    private:
        SubscriberHandle() noexcept;
        void invokeDeregister();

        std::function<void(Id)> deregister;
        impl::cEventStatusPtr eventStatus;
        bool disposed = false;
        Id id{};
    };


    /**
     * @brief C# inspired event class that manages a list of event handlers that can be invoked together
     * @tparam Args Arguments types of the event handler functions
     */
    template<typename ...Args>
    class SubscribableEvent {
        using handler = std::function<void(Args...)>;
        using HandlerId = SubscriberHandle::Id;
    public:
        SubscribableEvent() noexcept: eventStatus(std::make_shared<impl::EventStatus>()) {}

        SubscribableEvent(const SubscribableEvent &) = delete;
        SubscribableEvent &operator=(const SubscribableEvent &) = delete;
        SubscribableEvent(SubscribableEvent &&) noexcept = default;
        SubscribableEvent &operator=(SubscribableEvent &&) noexcept = default;

        ~SubscribableEvent() {
            eventStatus->invalidate();
        }

        /**
         * Adds a functor to the event handler list. The functor is called when trigger is invoked
         * @tparam Handler type of handler functor
         * @param handlerFunction functor to be subscribed to the event
         */
        template<typename Handler>
        void subscribe_unhandled(Handler &&handlerFunction) {
            subscribe(std::forward<Handler>(handlerFunction), true);
        }

        /**
         * @copydoc subscribe_unhandled
         * @return handle to the event that is used to unsubscribe the handler (see SubscriberHandle)
         * @note if the returned value is discarded the handler will be unregistered immediately after the call
         */
        template<typename Handler>
        [[nodiscard]]SubscriberHandle subscribe_handled(Handler &&handlerFunction) {
            return subscribe(std::forward<Handler>(handlerFunction), false);
        }

        /**
         * Triggers the event. All subscribed event handlers are invoked with the provided arguments
         * @param args arguments to the event handlers
         */
        template<typename ...InvokeArgs>
        void trigger(InvokeArgs&&... args) const {
            for (const auto &[handler, _] : handlers) {
                handler(std::forward<InvokeArgs>(args)...);
            }
        }

    private:
        template<typename Handler>
        SubscriberHandle subscribe(Handler &&handlerFunction, bool discardHandler) {
            static_assert(std::is_invocable_r_v<void, Handler, Args...>, "invalid event handler signature");
            handlers.emplace_back(std::forward<Handler>(handlerFunction), handlerId);
            return SubscriberHandle([this](auto id) { unsubscribe(id);}, eventStatus, handlerId++, discardHandler);
        }

        void unsubscribe(HandlerId id) {
            auto res = std::ranges::find_if(handlers, [id](const auto &p) { return p.second == id; });
            if (res != handlers.end()) {
                std::swap(*res, handlers.back());
                handlers.pop_back();
            }
        }

        std::vector<std::pair<handler, HandlerId>> handlers{};
        impl::EventStatusPtr eventStatus;
        HandlerId handlerId{};
    };
}

#endif //SCHEDCL_SUBSCRIBABLEEVENT_HPP
