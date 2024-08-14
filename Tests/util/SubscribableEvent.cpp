/**
 * @author Tim Luchterhand
 * @date 14.04.23.
 */

#include <gtest/gtest.h>

#include "util/SubscribableEvent.hpp"

struct TestHandle : public tempo::SubscriberHandle {
    TestHandle() = default;
    explicit TestHandle(int): tempo::SubscriberHandle(Subscribe) {}
};

TEST(util, SubscriberHandle_basic) {
    TestHandle handle;
    EXPECT_FALSE(handle.isSubscribed());
    TestHandle handle1(5);
    EXPECT_TRUE(handle1.isSubscribed());
}

TEST(util, SubscriberHandle_dispose) {
    TestHandle handle(1);
    handle.unregister();
    EXPECT_FALSE(handle.isSubscribed());
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
    EXPECT_TRUE(handle.isSubscribed());
    event.trigger(1);
    EXPECT_EQ(variable, 1);
    EXPECT_EQ(variable1, 1);
    EXPECT_TRUE(handle.isSubscribed());
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

struct TestSubscriber {
    tempo::SubscriberHandle handle;
    int value = 0;

    explicit TestSubscriber(tempo::SubscribableEvent<int> &event) : handle(
            event.subscribe_handled([this](int val) { value = val; })) {}
};

TEST(util, SubscribableEvent_auto_unsubscribe1) {
    using namespace tempo;
    SubscribableEvent<int> event;
    {
        TestSubscriber subscriber(event);
        event.trigger(17);
        EXPECT_EQ(subscriber.value, 17);
    }

    event.trigger(12);
}

TEST(util, SubscribableEvent_unsubscribe_ub_safety) {
    using namespace tempo;
    SubscriberHandle handle;
    {
        SubscribableEvent<> event;
        handle = event.subscribe_handled([]() {});
    }

    EXPECT_FALSE(handle.isSubscribed());
}

TEST(util, SubscribableEvent_move) {
    using namespace tempo;
    SubscribableEvent<int> event;
    int val1 = 0;
    int val2 = 0;
    event.subscribe_unhandled([&val1](int v) { val1 = v; });
    auto handle = event.subscribe_handled([&val2](int v) { val2 = v; });
    event.trigger(14);
    EXPECT_EQ(val1, 14);
    EXPECT_EQ(val2, 14);
    auto event1 = std::move(event);
    event.trigger(17);
    EXPECT_EQ(val1, 14);
    EXPECT_EQ(val2, 14);
    EXPECT_TRUE(handle.isSubscribed());
    event1.trigger(17);
    EXPECT_EQ(val1, 17);
    EXPECT_EQ(val2, 17);
    handle.unregister();
    handle = event1.subscribe_handled([&val2](int v) {val2 = v * 2;});
    EXPECT_TRUE(handle.isSubscribed());
    event1.trigger(4);
    EXPECT_EQ(val1, 4);
    EXPECT_EQ(val2, 8);
    event = std::move(event1);
    EXPECT_TRUE(handle.isSubscribed());
    event1.trigger(9);
    EXPECT_EQ(val1, 4);
    EXPECT_EQ(val2, 8);
    event.trigger(6);
    EXPECT_EQ(val1, 6);
    EXPECT_EQ(val2, 12);
}
TEST(util, SubscribableEvent_move1) {
    using namespace tempo;
    SubscribableEvent<int> event;
    int val1 = 0;
    int val2 = 0;
    auto handle = event.subscribe_handled([&val1](int v) { val1 = v; });
    SubscribableEvent<int> event1;
    auto handle1 = event1.subscribe_handled([&val2](int v) { val2 = v; });
    event = std::move(event1);
    EXPECT_FALSE(handle.isSubscribed());
    event1.trigger(17);
    EXPECT_EQ(val1, 0);
    EXPECT_EQ(val2, 0);
    event.trigger(4);
    EXPECT_EQ(val1, 0);
    EXPECT_EQ(val2, 4);
}
