
#ifndef _TEMPO_EXPLANATION_HPP
#define _TEMPO_EXPLANATION_HPP

#include <vector>

#include "Global.hpp"

namespace tempo {

template <typename T> struct Literal;

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

template <typename T> class NewExplainer {
protected:
  int cons_id{-1};

public:
  virtual void xplain(const Literal<T>, const hint, std::vector<Literal<T>> &);

  virtual std::ostream &print_reason(std::ostream &, const hint) const;

  // for introspection
  virtual int getType() const;

  virtual ~NewExplainer() = default;
  NewExplainer() = default;
  NewExplainer(NewExplainer &) = default;
  NewExplainer(NewExplainer &&) = default;
  NewExplainer &operator=(NewExplainer &) = default;
  NewExplainer &operator=(NewExplainer &&) = default;

  int id() const { return cons_id; }
};

template <typename T> class NewExplanation {

public:
  NewExplanation(NewExplainer<T> *e, hint h);
  void explain(const Literal<T> l, std::vector<Literal<T>> &Cl);
  std::ostream &display(std::ostream &os) const;
  bool operator==(const NewExplanation<T> &e) const;

  int getType() const;

  // private:
  NewExplainer<T> *expl{NULL};
  hint the_hint{NoHint};
};

template <typename T>
void NewExplainer<T>::xplain(const Literal<T>, const hint,
                             std::vector<Literal<T>> &) {}

template <typename T>
std::ostream &NewExplainer<T>::print_reason(std::ostream &os,
                                            const hint) const {
  os << "no reason";
  return os;
}

template <typename T>
void NewExplanation<T>::explain(const Literal<T> l,
                                std::vector<Literal<T>> &Cl) {
  expl->xplain(l, the_hint, Cl);
}

template <typename T> int NewExplainer<T>::getType() const { return NOEXPL; }

template <typename T>
std::ostream &NewExplanation<T>::display(std::ostream &os) const {

  assert(expl != NULL);

  // os << "b/c ";
  expl->print_reason(os, the_hint);
  return os;
}

template <typename T>
bool NewExplanation<T>::operator==(const NewExplanation<T> &e) const {
  return expl == e.expl and the_hint == e.the_hint;
}

template <typename T>
NewExplanation<T>::NewExplanation(NewExplainer<T> *e, hint h)
    : expl(e), the_hint(h) {}

// Explanation Constant::NoReason = Explanation(new Explainer(), NoHint);

template <typename T> int NewExplanation<T>::getType() const {
  return expl->getType();
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const NewExplanation<T> &x) {
  return x.display(os);
}

std::ostream &operator<<(std::ostream &os, const Explanation &x);

} // namespace tempo

#endif
