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
        struct Token{};
    public:
        SubscriberHandle(const SubscriberHandle &) = default;
        SubscriberHandle &operator=(const SubscriberHandle &) = default;
        SubscriberHandle(SubscriberHandle &&other) noexcept = default;
        SubscriberHandle &operator=(SubscriberHandle &&other) = default;

        /**
         * Manually unregister the associated event handler
         */
        void unregister();

        /**
         * DTor. Automatically unregisters the associated event handler if not already disposed.
         */
        ~SubscriberHandle();

        /**
         * Ctor. Creates an empty handle
         */
        constexpr SubscriberHandle() noexcept = default;

        /**
         * Whether the handle refers to a valid event
         * @return true if handle refers to an active event, false otherwise
         */
        [[nodiscard]] bool isSubscribed() const noexcept;

    protected:

        static constexpr Token Subscribe{};

        /**
         * Ctor
         * creates a handle with specified id
         * @param id id value of the handle
         */
        SubscriberHandle(Token);

    private:
        static constexpr char Alive = 1;
        std::shared_ptr<char> alive;
    };


    /**
     * @brief C# inspired event class that manages a list of event handlers that can be invoked together
     * @tparam Args Arguments types of the event handler functions
     */
    template<typename ...Args>
    class SubscribableEvent {
        using handler = std::function<void(Args...)>;
    public:
        constexpr SubscribableEvent() noexcept = default;

        SubscribableEvent(const SubscribableEvent &) = delete;
        SubscribableEvent &operator=(const SubscribableEvent &) = delete;
        SubscribableEvent(SubscribableEvent &&) noexcept = default;
        SubscribableEvent &operator=(SubscribableEvent &&) noexcept = default;

        ~SubscribableEvent() {
            for (auto &[_, handle]: handlers) {
                handle.unregister();
            }
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
            auto it = handlers.begin();
            auto end = handlers.end();
            while (not handlers.empty() and it != end) {
                if (it->second.isSubscribed()) {
                    it->first(std::forward<InvokeArgs>(args)...);
                    ++it;
                } else {
                    std::swap(*it, handlers.back());
                    handlers.pop_back();
                    end = handlers.end();
                }
            }
        }

    private:
        template<typename Handler>
        SubscriberHandle subscribe(Handler &&handlerFunction, bool discardHandler) {
            static_assert(std::is_invocable_r_v<void, Handler, Args...>, "invalid event handler signature");
            handlers.emplace_back(std::forward<Handler>(handlerFunction), SubscriberHandle(SubscriberHandle::Subscribe));
            return discardHandler ? SubscriberHandle() : handlers.back().second;
        }

        mutable std::vector<std::pair<handler, SubscriberHandle>> handlers{};
    };
}

#endif //SCHEDCL_SUBSCRIBABLEEVENT_HPP
