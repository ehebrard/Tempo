/**
 * @author Tim Luchterhand
 * @date 14.04.23.
 */

#include <gtest/gtest.h>
#include <Iterators.hpp>
#include <functional>
#include <concepts>

#include "util/SubscribableEvent.hpp"

struct TestHandle : public tempo::SubscriberHandle {
    TestHandle() = default;
    explicit TestHandle(int): tempo::SubscriberHandle(Subscribe) {}
};

TEST(util, ActiveList_basic) {
    using namespace tempo;
    detail::ActiveList<int> list;
    EXPECT_TRUE(list.empty());
    EXPECT_EQ(list.begin(), list.end());
    EXPECT_EQ(list.size(), 0);
}

TEST(util, ActiveList_basic1) {
    using namespace tempo;
    detail::ActiveList<int> list{1, 2, 3};
    EXPECT_EQ(list.size(), 3);
    for (auto [gt, i] : iterators::enumerate(list, 1)) {
        EXPECT_EQ(gt, i);
    }
}

TEST(util, ActiveList_inactive) {
    using namespace tempo;
    detail::ActiveList<int> list{1, 2, 3};
    list.markInactive(list.begin() + 1);
    EXPECT_EQ(list.size(), 2);
    EXPECT_EQ(*list.begin(), 1);
    EXPECT_EQ(list.back(), 3);
    list.markInactive(list.end());
    EXPECT_EQ(list.size(), 2);
    EXPECT_EQ(*list.begin(), 1);
    EXPECT_EQ(list.back(), 3);
    list.markInactive(list.begin() - 1);
    EXPECT_EQ(list.size(), 2);
    EXPECT_EQ(*list.begin(), 1);
    EXPECT_EQ(list.back(), 3);
}

enum class Action {
    Nothing,
    Dtor,
    MoveAssign,
    MoveCtor,
    CopyAssign,
    CopyCtor
};

class Data {
    std::function<void(Action)> fun;
public:
    int val;
    template<std::invocable<Action> Fun>
    explicit Data(Fun &&f, int i) : fun(std::forward<Fun>(f)), val(i) {}

    Data(const Data &data): fun(data.fun), val(data.val) {
        fun(Action::CopyCtor);
    }

    Data(Data &&data) noexcept: fun(std::move(data.fun)), val(data.val) {
        fun(Action::MoveCtor);
    }

    Data &operator=(const Data &data) {
        if (&data != this) {
            fun = data.fun;
            val = data.val;
            fun(Action::CopyAssign);
        }

        return *this;
    }

    Data &operator=(Data &&data) noexcept {
        if(&data != this) {
            fun = std::move(data.fun);
            val = data.val;
            fun(Action::MoveAssign);
        }

        return *this;
    }

    ~Data() {
        if (fun) {
            fun(Action::Dtor);
        }
    }

};

TEST(util, ActiveList_inactive_overwrite) {
    using namespace tempo;
    using enum Action;
    Action obj1 = Nothing;
    Action obj2 = Nothing;
    detail::ActiveList<Data> list;
    list.add([&obj1](auto action) { obj1 = action; }, 1);
    EXPECT_EQ(obj1, Nothing);
    list.add([&obj2](auto action) { obj2 = action; }, 2);
    EXPECT_NE(obj1, Dtor);
    list.markInactive(list.begin());
    EXPECT_NE(obj1, Dtor);
    EXPECT_NE(obj2, Dtor);
    EXPECT_EQ(list.size(), 1);
    obj2 = Nothing;
    list.add([](auto){}, 3);
    EXPECT_EQ(obj1, Dtor);
    EXPECT_EQ(obj2, Nothing);
}

TEST(util, ActiveList_move) {
    using namespace tempo;
    using enum Action;
    Action obj1, obj2, obj3;
    detail::ActiveList<Data> list;
    list.add([&obj1](auto a) { obj1 = a; }, 1);
    list.add([&obj2](auto a) { obj2 = a; }, 2);
    list.markInactive(list.begin());
    obj1 = Nothing; obj2 = Nothing; obj3 = Nothing;
    auto moved = std::move(list);
    EXPECT_EQ(obj1, Nothing);
    EXPECT_EQ(obj2, Nothing);
    EXPECT_TRUE(list.empty());
    EXPECT_EQ(moved.size(), 1);
    EXPECT_EQ(moved.back().val, 2);
    moved.add([&obj3](auto a) { obj3 = a; }, 3);
    EXPECT_EQ(obj1, Dtor);
    EXPECT_EQ(obj2, Nothing);
    EXPECT_EQ(moved.back().val, 3);
}

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

TEST(util, SubscribableEvent_unsubscribe1) {
    using namespace tempo;
    int variable = 0;
    int variable1 = 0;
    int variable2 = 0;
    SubscribableEvent<int> event;
    event.subscribe_unhandled([&variable](int val) { variable = val; });
    auto handle1 = event.subscribe_handled([&variable1](int val) { variable1 = val; });
    event.trigger(14);
    EXPECT_EQ(variable, 14);
    EXPECT_EQ(variable1, 14);
    EXPECT_EQ(variable2, 0);
    handle1.unregister();
    event.trigger(8);
    EXPECT_EQ(variable, 8);
    EXPECT_EQ(variable1, 14);
    EXPECT_EQ(variable2, 0);
    auto handle2 = event.subscribe_handled([&variable2](int val) { variable2 = val; });
    event.trigger(17);
    EXPECT_EQ(variable, 17);
    EXPECT_EQ(variable1, 14);
    EXPECT_EQ(variable2, 17);
    handle1 = event.subscribe_handled([&variable1](int val) { variable1 = val; });
    event.trigger(3);
    EXPECT_EQ(variable, 3);
    EXPECT_EQ(variable1, 3);
    EXPECT_EQ(variable2, 3);
    handle1.unregister();
    handle2.unregister();
    event.trigger(6);
    EXPECT_EQ(variable, 6);
    EXPECT_EQ(variable1, 3);
    EXPECT_EQ(variable2, 3);
    event.subscribe_unhandled([](auto){});
    event.subscribe_unhandled([&variable1](int val) { variable1 = val; });
    event.subscribe_unhandled([&variable2](int val) { variable2 = val; });
    event.trigger(16);
    EXPECT_EQ(variable, 16);
    EXPECT_EQ(variable1, 16);
    EXPECT_EQ(variable2, 16);
}

TEST(util, SubscribableEvent_call_cycle) {
    using namespace tempo;
    SubscribableEvent<> event;
    SubscriberHandle handle;
    int counter = 0;
    handle = event.subscribe_handled([&] {
        ++counter;
        handle.unregister();
        event.trigger();
    });

    event.subscribe_unhandled([]{});
    event.trigger();
    EXPECT_EQ(counter, 1);
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
