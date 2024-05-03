/**
 * @author Tim Luchterhand
 * @date 14.04.23.
 */

#include <gtest/gtest.h>
#include <optional>

#include "util/SubscribableEvent.hpp"

class TestHandle : public tempo::SubscriberHandle {
public:
    template<typename Fun>
    TestHandle(Fun &&deregister, tempo::impl::cEventStatusPtr status, Id id, bool disposed):
            SubscriberHandle(std::forward<Fun>(deregister), std::move(status), id, disposed) {}
};

TEST(util, EventStatus_basic) {
    using namespace tempo;
    impl::EventStatus status;
    EXPECT_TRUE(status.isAlive());
    status.invalidate();
    EXPECT_FALSE(status.isAlive());
}

TEST(util, SubscriberHandle_basic) {
    TestHandle handle([](auto) {}, {}, 3, false);
    EXPECT_FALSE(handle.isDisposed());
    handle = TestHandle([](auto) {}, {}, 3, true);
    EXPECT_TRUE(handle.isDisposed());
}

TEST(util, SubscriberHandle_dispose) {
    TestHandle handle([](auto) {}, {}, 3, false);
    handle.dispose();
    EXPECT_TRUE(handle.isDisposed());
}

TEST(util, SubscriberHandle_move) {
    TestHandle handle([](auto) {}, {}, 3, false);
    TestHandle handle1 = std::move(handle);
    EXPECT_FALSE(handle1.isDisposed());
    EXPECT_TRUE(handle.isDisposed());
    handle = std::move(handle1);
    EXPECT_FALSE(handle.isDisposed());
    EXPECT_TRUE(handle1.isDisposed());
    TestHandle handle2([](auto) {}, {}, 4, false);
    handle = std::move(handle2);
    EXPECT_TRUE(handle2.isDisposed());
    EXPECT_FALSE(handle.isDisposed());
}

TEST(util, SubscribableEvent_basic) {
    using namespace tempo;
    SubscribableEvent<bool> event;
    bool success1 = false;
    bool success2 = false;
    bool success3 = false;
    event.subscribe_unhandled([&success1](bool first) { if (first) { success1 = true; }});
    event.subscribe_unhandled([&success2](bool first) { if (not first) { success2 = true; }});
    event.subscribe_unhandled([&success3](bool) { success3 = true; });
    event.trigger(true);
    EXPECT_TRUE(success1 and success3);
    EXPECT_FALSE(success2);
    event.trigger(false);
    EXPECT_TRUE(success2);
}


TEST(util, SubscribableEvent_unsubscribe) {
    using namespace tempo;
    int variable = 0;
    int variable1 = 0;
    SubscribableEvent<int> event;
    auto handle = event.subscribe_handled([&variable](int val) { variable = val; });
    event.subscribe_unhandled([&variable1](int val) { variable1 = val; });
    EXPECT_FALSE(handle.isDisposed());
    event.trigger(1);
    EXPECT_EQ(variable, 1);
    EXPECT_EQ(variable1, 1);
    EXPECT_FALSE(handle.isDisposed());
    event.trigger(17);
    EXPECT_EQ(variable, 17);
    EXPECT_EQ(variable1, 17);
    handle.unregister();
    event.trigger(29);
    EXPECT_EQ(variable, 17);
    EXPECT_EQ(variable1, 29);
}

TEST(util, SubscribableEvent_auto_unsubscribe) {
    using namespace tempo;
    int variable = 0;
    SubscribableEvent<int> event;
    {
        auto handle = event.subscribe_handled([&variable](int val) { variable = val; });
        event.trigger(3);
        EXPECT_EQ(variable, 3);
    }

    event.trigger(19);
    EXPECT_EQ(variable, 3);
}

TEST(util, SubscribableEvent_unsubscribe_disposed) {
    using namespace tempo;
    int variable = 0;
    SubscribableEvent<int> event;
    auto handle = event.subscribe_handled([&variable](int val) { variable = val; });
    handle.dispose();
    EXPECT_TRUE(handle.isDisposed());
    handle.unregister();
    event.trigger(11);
    EXPECT_EQ(variable, 11);
}

TEST(util, SubscribableEvent_unsubscribe_ub_safety) {
    using namespace tempo;
    std::optional<SubscriberHandle> handle;
    {
        SubscribableEvent<> event;
        handle = event.subscribe_handled([]() {});
    }

    handle->unregister();
}
