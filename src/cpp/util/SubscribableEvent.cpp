//
// Created by tim on 22.03.23.
//

#include "util/SubscribableEvent.hpp"

namespace tempo {

    void SubscriberHandle::dispose() noexcept {
        disposed = true;
    }

    bool SubscriberHandle::isDisposed() const noexcept {
        return disposed;
    }

    void SubscriberHandle::unregister() {
        invokeDeregister();
    }

    SubscriberHandle::~SubscriberHandle() {
        invokeDeregister();
    }

    void SubscriberHandle::invokeDeregister() {
        if (not disposed and nullptr != eventStatus and eventStatus->isAlive()) {
            deregister(id);
            dispose();
        }
    }

    SubscriberHandle::SubscriberHandle(SubscriberHandle &&other) noexcept: SubscriberHandle() {
        using std::swap;
        swap(*this, other);
    }

    SubscriberHandle &SubscriberHandle::operator=(SubscriberHandle &&other) noexcept {
        using std::swap;
        SubscriberHandle tmp = std::move(other);
        swap(*this, tmp);
        return *this;
    }

    SubscriberHandle::SubscriberHandle() noexcept: deregister(), eventStatus(), disposed(true), id() {}

    void swap(SubscriberHandle &lhs, SubscriberHandle &rhs) noexcept {
        using std::swap;
        swap(lhs.deregister, rhs.deregister);
        swap(lhs.eventStatus, rhs.eventStatus);
        swap(lhs.disposed, rhs.disposed);
        swap(lhs.id, rhs.id);
    }
}
