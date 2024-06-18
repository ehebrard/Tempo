

#ifndef _TEMPO_CLAUSE_HPP
#define _TEMPO_CLAUSE_HPP

#include <vector>

#include "Explanation.hpp"
#include "Global.hpp"
#include "Literal.hpp"

namespace tempo {


template <typename T> class Clause : public std::vector<Literal<T>> {

public:
  const static Clause<T> *Empty;

  Clause(const int i = -1);

  int id;
  index_t watched_index[2];

  Literal<T> watched(const bool) const;
  bool watch_rank(const Literal<T>) const;
  size_t watch_index(const bool) const;

  //  bool operator<(const Clause &) const;
  //  bool operator==(const Clause &) const;

  //    // explanation stuff
  //     void xplain(const lit l, const hint h, std::vector<lit> &Cl) ;
  //     std::ostream &print_reason(std::ostream &, const hint) const;
  //     int getType() const;

  std::ostream &display(std::ostream &) const;
};

template <typename T>
Clause<T>::Clause(const int i) : std::vector<Literal<T>>(), id(i) {
  watched_index[0] = 0;
  watched_index[1] = 1;
}


template <typename T> Literal<T> Clause<T>::watched(const bool r) const {
  return this->operator[](watched_index[r]);
}

template <typename T> bool Clause<T>::watch_rank(const Literal<T> l) const {
  //  return this->operator[](watched_index[1]) == l;
  //    auto p{this->operator[](watched_index[1])};
  //    return (p.isNumeric() == l.isNumeric()) and (p.variable() ==
  //    l.variable());
    auto p{this->operator[](watched_index[1])};
  return l.sameVariable(p) and l.sign() != p.sign();
}

template <typename T> size_t Clause<T>::watch_index(const bool r) const {
  return watched_index[r];
}

template <typename T>
std::ostream &Clause<T>::display(std::ostream &os) const {
  os << id << ":";
  if (this->size() == 0)
    os << "()";
  else {
    os << "(" << this->operator[](0);
    for (size_t i{1}; i < this->size(); ++i) {
      os << ", " << this->operator[](i);
    }
      os << ")";
//      os << ") w=" << watched(0) << " & " << watched(1) << "]";
  }
  return os;
}

template <typename T>
const Clause<T> *Clause<T>::Empty = new Clause<T>(-1);

template <typename T>
std::ostream &operator<<(std::ostream &os, const Clause<T> &x) {
  return x.display(os);
}
}

#endif // _TEMPO_CLAUSE_HPP

