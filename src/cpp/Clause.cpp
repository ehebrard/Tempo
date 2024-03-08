

#include <iostream>
#include <assert.h>

#include "Clause.hpp"


tempo::Clause::Clause(const int i) : std::vector<lit>(), id(i) {
  watcher_index[0] = 0;
  watcher_index[1] = 1;
}

tempo::Clause::~Clause() {}

bool tempo::Clause::operator<(const Clause &c) const {
  if (size() == c.size()) {
    for (size_t i{0}; i < size(); ++i) {
      if ((*this)[i] < c[i])
        return true;
      else if ((*this)[i] > c[i])
        return false;
    }
  }
  return size() < c.size();
}

bool tempo::Clause::operator==(const Clause &c) const {
  if (size() == c.size()) {
    for (size_t i{0}; i < size(); ++i) {
      if ((*this)[i] != c[i])
        return false;
    }
    return true;
  }
  return false;
}

tempo::lit tempo::Clause::watcher(const bool r) const {
  return this->operator[](watcher_index[r]);
}

bool tempo::Clause::watch_rank(const lit l) const {
  return this->operator[](watcher_index[1]) == l;
}

size_t tempo::Clause::watch_index(const bool r) const {
  return watcher_index[r];
}

std::ostream &tempo::Clause::display(std::ostream &os) const {
  if(size()==0)
    os << "()";
  else {
    os << "(" << this->operator[](0);
    for(size_t i{1}; i<size(); ++i) {
      os << ", " << this->operator[](i);
    }
    os << ")";
  }
  return os;
}

//void tempo::Clause::xplain(const lit l, const hint h, std::vector<lit> &Cl) {
////    assert((*this)[h] == l);
////    for(auto p : *this)
////        if(p != l)
////            Cl.push_back(p);
//    auto  n{size()};
//    for(size_t i=1; i<n; ++i) {
//        Cl.push_back(this->operator[]((h+i)%n));
//    }
//}
//
//std::ostream &tempo::Clause::print_reason(std::ostream &os, const hint) const {
//    return display(os);
//}
//
//int tempo::Clause::getType() const {
//    return CLAUSEEXPL;
//}

const tempo::Clause *tempo::Clause::Empty = new tempo::Clause(-1);


std::ostream &tempo::operator<<(std::ostream &os, const tempo::Clause &x) {
  return x.display(os);
}

