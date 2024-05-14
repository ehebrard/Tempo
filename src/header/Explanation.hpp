
#ifndef _TEMPO_EXPLANATION_HPP
#define _TEMPO_EXPLANATION_HPP

#include <vector>

#include "Global.hpp"

namespace tempo {



/* Explanation system:

An Explainer is a class with the methods 'xplain(const lit l, const hint h, std::vector<lit> &Cl)' that inserts an explanation for the literal l to the clause Cl

Typically, each constraint is an explainer, and must store the info necessary to access/compute the reason 'h' can be used to store or access the explanation

TemporalGraph are also explainers, for bounds

The static object "NoReason" is the empty reason

*/
class Explainer {
protected:
  int cons_id{-1};

public:
  virtual void xplain(const lit, const hint, std::vector<lit>&) ;

  virtual std::ostream &print_reason(std::ostream &, const hint) const;

  // for introspection
  virtual int getType() const;

  virtual ~Explainer() = default;
  Explainer() = default;
  Explainer(Explainer &) = default;
  Explainer( Explainer &&) = default;
  Explainer &operator=(Explainer &) = default;
  Explainer &operator=( Explainer &&) = default;

  int id() const { return cons_id; }
};

class Explanation {

public:
  Explanation(Explainer *e, hint h);
  void explain(const lit l, std::vector<lit> &Cl) ;
  std::ostream &display(std::ostream &os) const;
  bool operator==(const Explanation &e) const;

  int getType() const;

  // private:
  Explainer *expl{NULL};
  hint the_hint{NoHint};
};

// class Constant {
// public:
//   static Explanation NoReason;
// };

std::ostream &operator<<(std::ostream &os, const Explanation &x);
//std::ostream &operator<<(std::ostream &os, Explanation &x);

} // namespace schedcl

#endif
