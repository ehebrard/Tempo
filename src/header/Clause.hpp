

#ifndef _TEMPO_CLAUSE_HPP
#define _TEMPO_CLAUSE_HPP

#include <vector>

#include "Explanation.hpp"
#include "Global.hpp"
#include "Literal.hpp"

namespace tempo {

class Clause : public std::vector<lit> { //}, public Explainer {

public:

  const static Clause* Empty;

  Clause(const int i = -1);
  virtual ~Clause();

  int id;
  index_t watcher_index[2];

  lit watcher(const bool) const;
  bool watch_rank(const lit) const;
size_t watch_index(const bool) const;

  bool operator<(const Clause &) const;
  bool operator==(const Clause &) const;
    
//    // explanation stuff
//     void xplain(const lit l, const hint h, std::vector<lit> &Cl) ;
//     std::ostream &print_reason(std::ostream &, const hint) const;
//     int getType() const;

  std::ostream& display(std::ostream &) const;

};

std::ostream &operator<<(std::ostream &, const Clause &);

template <typename T> class NewClause : public std::vector<Literal<T>> {

public:
  const static NewClause<T> *Empty;

  NewClause(const int i = -1);

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
NewClause<T>::NewClause(const int i) : std::vector<Literal<T>>(), id(i) {
  watched_index[0] = 0;
  watched_index[1] = 1;
}

// bool NewClause<T>::operator<(const NewClause &c) const {
//   if (size() == c.size()) {
//     for (size_t i{0}; i < size(); ++i) {
//       if ((*this)[i] < c[i])
//         return true;
//       else if ((*this)[i] > c[i])
//         return false;
//     }
//   }
//   return size() < c.size();
// }
//
// bool NewClause<T>::operator==(const Clause &c) const {
//   if (size() == c.size()) {
//     for (size_t i{0}; i < size(); ++i) {
//       if ((*this)[i] != c[i])
//         return false;
//     }
//     return true;
//   }
//   return false;
// }

template <typename T> Literal<T> NewClause<T>::watched(const bool r) const {
  return this->operator[](watched_index[r]);
}

template <typename T> bool NewClause<T>::watch_rank(const Literal<T> l) const {
  //  return this->operator[](watched_index[1]) == l;
  //    auto p{this->operator[](watched_index[1])};
  //    return (p.isNumeric() == l.isNumeric()) and (p.variable() ==
  //    l.variable());
  return l.sameVariable(this->operator[](watched_index[1]));
}

template <typename T> size_t NewClause<T>::watch_index(const bool r) const {
  return watched_index[r];
}

template <typename T>
std::ostream &NewClause<T>::display(std::ostream &os) const {
  if (this->size() == 0)
    os << "()";
  else {
    os << "(" << this->operator[](0);
    for (size_t i{1}; i < this->size(); ++i) {
      os << ", " << this->operator[](i);
    }
    os << ")";
  }
  return os;
}

template <typename T>
const NewClause<T> *NewClause<T>::Empty = new NewClause<T>(-1);

template <typename T>
std::ostream &operator<<(std::ostream &os, const NewClause<T> &x) {
  return x.display(os);
}
}

#endif // _TEMPO_CLAUSE_HPP

