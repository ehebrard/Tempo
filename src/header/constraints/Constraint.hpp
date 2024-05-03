#ifndef _TEMPO_CONSTRAINT_HPP
#define _TEMPO_CONSTRAINT_HPP

#include <ostream>

#include "Explanation.hpp"
#include "Global.hpp"

namespace tempo {

class Constraint : public Explainer {

public:
  //    int id() { return cons_id; }

  virtual ~Constraint() = default;
  Constraint() = default;
  Constraint(const Constraint &) = delete;
  Constraint(Constraint &&) = delete;
  Constraint &operator=(const Constraint &) = delete;
  Constraint &operator=(Constraint &&) = delete;

  Priority priority = Priority::High;
  bool idempotent{false};

  //
  virtual void post(const int idx) = 0;
  // propagate the constraint
  virtual void propagate() = 0;
  // notify a change (with the literal and it's variable rank in the scope)
  virtual bool notify_bound(const lit, const int) { return false; }
  virtual bool notify_edge(const lit, const int) { return false; }

  virtual std::ostream &display(std::ostream &os) const = 0;

  // protected:
  //   int cons_id;
};

std::ostream &operator<<(std::ostream &os, const Constraint &x);

} // namespace tempo

#endif // _TEMPO_CONSTRAINT_HPP
