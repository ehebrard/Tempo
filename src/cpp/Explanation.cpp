
#include <iostream>

#include <assert.h>

#include "Explanation.hpp"


using namespace tempo;

void Explainer::xplain(const lit, const hint, std::vector<lit> &) {}

std::ostream &Explainer::print_reason(std::ostream &os, const hint) const {
  os << "no reason";
  return os;
}

 void Explanation::explain(const lit l, std::vector<lit> &Cl) {
   expl->xplain(l, the_hint, Cl);
 }

 int Explainer::getType() const { return NOEXPL; }

 std::ostream &Explanation::display(std::ostream &os) const {

   assert(expl != NULL);

   // os << "b/c ";
   expl->print_reason(os, the_hint);
   return os;
 }

bool Explanation::operator==(const Explanation &e) const {
  return expl == e.expl and the_hint == e.the_hint;
}

Explanation::Explanation(Explainer *e, hint h) : expl(e), the_hint(h) {}

// Explanation Constant::NoReason = Explanation(new Explainer(), NoHint);

int Explanation::getType() const { return expl->getType(); }

std::ostream &tempo::operator<<(std::ostream &os, Explanation &x) {
  return x.display(os);
}
