
#include "Constraint.hpp"

std::ostream &tempo::operator<<(std::ostream &os, const tempo::Constraint &x) {
  return x.display(os);
}
