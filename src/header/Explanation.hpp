/************************************************
 * Tempo Explanation.hpp
 *
 * Copyright 2024 Emmanuel Hebrard
 *
 * Tempo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 * Tempo is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Tempo.  If not, see <http://www.gnu.org/licenses/>.
 *
 ***********************************************/

#ifndef _TEMPO_EXPLANATION_HPP
#define _TEMPO_EXPLANATION_HPP

#include <vector>

#include "Global.hpp"

namespace tempo {

template <typename T> struct Literal;


//! Explainer interface
/*!
 An Explainer is a class with the methods 'xplain(const lit l, const hint h, std::vector<lit> &Cl)' that inserts an explanation for the literal l to the clause Cl

 Typically, each constraint is an explainer, and must store the info necessary to access/compute the reason 'h' can be used to store or access the explanation
 */
template <typename T> class Explainer {
protected:
  int cons_id{-1};

public:
  virtual void xplain(const Literal<T>, const hint, std::vector<Literal<T>> &);

  virtual std::ostream &print_reason(std::ostream &, const hint) const;

  // for introspection
//  virtual int getType() const;

  virtual ~Explainer() = default;
  Explainer() = default;
  Explainer(Explainer &) = default;
  Explainer(Explainer &&) = default;
  Explainer &operator=(Explainer &) = default;
  Explainer &operator=(Explainer &&) = default;

  int id() const { return cons_id; }
};

//! Explanation class
/*!
 An Explanation is a pointer to an Explainer  and a 'hint' to help it explain
  - e.g., it can hold the index of the explantion in a precomputed list, or it can give the necessary info to recompute it
 */
template <typename T> class Explanation {

public:
    Explanation() {}
  Explanation(Explainer<T> *e, hint h);
  void explain(const Literal<T> l, std::vector<Literal<T>> &Cl);
  std::ostream &display(std::ostream &os) const;
  bool operator==(const Explanation<T> &e) const;

//  int getType() const;

  // private:
  Explainer<T> *expl{NULL};
  hint the_hint{static_cast<int>(-1)};
};

template <typename T>
void Explainer<T>::xplain(const Literal<T>, const hint,
                             std::vector<Literal<T>> &) {}

template <typename T>
std::ostream &Explainer<T>::print_reason(std::ostream &os,
                                            const hint) const {
  os << "no reason";
  return os;
}

//template <typename T> int Explainer<T>::getType() const { return NOEXPL; }

template <typename T>
void Explanation<T>::explain(const Literal<T> l,
                                std::vector<Literal<T>> &Cl) {
  expl->xplain(l, the_hint, Cl);
}

template <typename T>
std::ostream &Explanation<T>::display(std::ostream &os) const {

  assert(expl != NULL);

  // os << "b/c ";
  expl->print_reason(os, the_hint);
  return os;
}

template <typename T>
bool Explanation<T>::operator==(const Explanation<T> &e) const {
  return expl == e.expl and the_hint == e.the_hint;
}

template <typename T>
Explanation<T>::Explanation(Explainer<T> *e, hint h)
    : expl(e), the_hint(h) {}

// Explanation Constant::NoReason = Explanation(new Explainer(), NoHint);

//template <typename T> int Explanation<T>::getType() const {
//  return expl->getType();
//}

template <typename T>
std::ostream &operator<<(std::ostream &os, const Explanation<T> &x) {
  return x.display(os);
}

} // namespace tempo

#endif
