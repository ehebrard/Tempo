#ifndef _TEMPO_CONSTRAINT_HPP
#define _TEMPO_CONSTRAINT_HPP

#include <ostream>

#include "Explanation.hpp"
#include "Global.hpp"

namespace tempo {


template <typename T> class Constraint : public Explainer<T> {

public:
  //    int id() { return cons_id; }

  virtual ~Constraint() = default;
  Constraint() = default;
  Constraint(const Constraint<T> &) = delete;
  Constraint(Constraint<T> &&) = delete;
  Constraint &operator=(const Constraint<T> &) = delete;
  Constraint &operator=(Constraint<T> &&) = delete;

  Priority priority = Priority::High;
  bool idempotent{false};

  //
  virtual void post(const int idx) = 0;
  // propagate the constraint
  virtual void propagate() = 0;
  // notify a change (with the literal and it's variable rank in the scope)
  virtual bool notify(const Literal<T>, const int) { return false; }
  //  virtual bool notify_edge(const Literal<T>, const int) { return false; }

  virtual std::ostream &display(std::ostream &os) const = 0;

  // protected:
  //   int cons_id;
};

template <typename T>
std::ostream &operator<<(std::ostream &os, const Constraint<T> &x) {
  return x.display(os);
}

} // namespace tempo

#endif // _TEMPO_CONSTRAINT_HPP
