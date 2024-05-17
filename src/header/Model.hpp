
#ifndef __TEMPO_MODEL_HPP
#define __TEMPO_MODEL_HPP


#include "Literal.hpp"

using namespace std;

namespace tempo {

template<typename T> class Solver;

/**********************************************
 * Job
 **********************************************/

template<typename T=int>
class NumericVar {

public:
  NumericVar(const var_t i) : _id_(i) {}

  T min(Solver<T> &sc) const;
  T max(Solver<T> &sc) const;

  Literal<T> operator>=(const T t) const;
  Literal<T> operator>(const T t) const;
  Literal<T> operator<=(const T t) const;
  Literal<T> operator<(const T t) const;

  var_t id() const { return _id_; }

  operator var_t() const { return _id_; }

  std::ostream &display(std::ostream &os) const;

  static bool isNumeric() { return true; }

protected:
  const var_t _id_;
};


template<typename T=int>
class BooleanVar {

public:
  BooleanVar(const var_t i) : _id_(i) {}

  Literal<T> operator==(const bool t) const;

  var_t id() const { return _id_; }

  operator var_t() const { return _id_; }

  std::ostream &display(std::ostream &os) const;

  static bool isNumeric() { return false; }

protected:
  const var_t _id_;
};

template<typename T>
Literal<T> BooleanVar<T>::operator==(const bool t) const {
  return Literal<T>(t, _id_, Constant::NoSemantic);
}


template<typename T=int>
class DisjunctVar : public BooleanVar<T> {

public:
  DisjunctVar(const var_t i, const info_t d) : BooleanVar<T>(i), _edge_id_(d) {}

  Literal<T> operator==(const bool t) const;

  std::ostream &display(std::ostream &os) const;

  static bool isNumeric() { return false; }

private:
  const info_t _edge_id_;
};

template<typename T>
Literal<T> DisjunctVar<T>::operator==(const bool t) const {
  return Literal<T>(t, BooleanVar<T>::_id_, _edge_id_ + t);
}

template<typename T=int>
class TemporalVar : public NumericVar<T> {

public:
  TemporalVar(const var_t i, const T o = 0) : NumericVar<T>(i), _offset(o) {}

  T earliest(Solver<T> &) const;
  T latest(Solver<T> &) const;

  Literal<T> after(const T t) const;
  Literal<T> before(const T t) const;

  DistanceConstraint<T> after(const TemporalVar<T> &e, const T t = 0) const;
  DistanceConstraint<T> before(const TemporalVar<T> &e, const T t = 0) const;

  T offset() const { return _offset; }

  std::ostream &display(std::ostream &os) const;

  static bool isNumeric() { return true; }

private:
  const T _offset;
};

template<typename T>
T NumericVar<T>::min(Solver<T>& s) const {
  return s.numeric.lower(_id_);
}

template<typename T>
T NumericVar<T>::max(Solver<T>& s) const {
  return s.numeric.upper(_id_);
}

template<typename T>
Literal<T> NumericVar<T>::operator<=(const T t) const {
  return leq<T>(_id_, t);
}

template<typename T>
Literal<T> NumericVar<T>::operator>=(const T t) const {
  return geq<T>(_id_, t);
}

template<typename T>
Literal<T> NumericVar<T>::operator<(const T t) const {
  return lt<T>(_id_, t);
}

template<typename T>
Literal<T> NumericVar<T>::operator>(const T t) const {
  return gt<T>(_id_, t);
}


template<typename T>
T TemporalVar<T>::earliest(Solver<T>& s) const {
  return NumericVar<T>::min(s) + _offset;
}

template<typename T>
T TemporalVar<T>::latest(Solver<T>& s) const {
  return NumericVar<T>::max(s) + _offset;
}

template<typename T>
Literal<T> TemporalVar<T>::after(const T t) const {
  return geq<T>(NumericVar<T>::_id_, t - _offset);
}

template<typename T>
Literal<T> TemporalVar<T>::before(const T t) const {
  return leq<T>(NumericVar<T>::_id_, t - _offset);
}

template<typename T>
DistanceConstraint<T> TemporalVar<T>::after(const TemporalVar<T>& e, const T t) const {
  return e.before(*this, t);
}

template<typename T>
DistanceConstraint<T> TemporalVar<T>::before(const TemporalVar<T>& e, const T t) const {
  return {e.id(), NumericVar<T>::_id_, e.offset() - _offset - t};
}

template<typename T>
std::ostream &BooleanVar<T>::display(std::ostream &os) const {
  os << "b" << id();
  return os;
}

template<typename T>
std::ostream &NumericVar<T>::display(std::ostream &os) const {
  os << "x" << id();
  return os;
}

template<typename T>
std::ostream &TemporalVar<T>::display(std::ostream &os) const {
  os << "x" << NumericVar<T>::id();
  return os;
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const BooleanVar<T> &x) {
  return x.display(os);
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const NumericVar<T> &x) {
  return x.display(os);
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const TemporalVar<T> &x) {
  return x.display(os);
}




template<typename T=int>
class Job {
public:
  Job(Solver<T> &s, const T mindur, const T maxdur);
    
    Job(const Job<T>&) = default;

    T getEarliestStart(Solver<T> &s) const;
    T getLatestStart(Solver<T> &s) const;
    T getEarliestEnd(Solver<T> &s) const;
    T getLatestEnd(Solver<T> &s) const;

    bool mustExist(Solver<T> &s) const;
    bool cannotExist(Solver<T> &s) const;

    T minDuration() const;
    T maxDuration() const;

    event getStart() const;
    event getEnd() const;

    int id() const;
    bool operator==(const Job<T> &t) const;

    std::ostream &display(std::ostream &os) const;

    TemporalVar<T> start;
    TemporalVar<T> end;

  private:
    T min_duration;
    T max_duration;
    BooleanVar<T> optional;
};

template <typename T=int> class DisjunctiveResource : public vector<Job<T>> {
public:
  using vector<Job<T>>::vector;
  template <typename C>
  void createOrderVariables(Solver<T> &solver, C &container);
};

template <typename T>
template <typename C>
void DisjunctiveResource<T>::createOrderVariables(Solver<T> &solver,
                                                  C &container) {
  for (auto a{this->begin()}; a != this->end(); ++a) {
    for (auto b{a + 1}; b != this->end(); ++b) {
      //                  container.insert(container.end(),
      //                  solver.newDisjunct(a->end.before(b->start),
      //                  b->end.before(a->start)));
      container.push_back(
          solver.newDisjunct(a->end.before(b->start), b->end.before(a->start)));
    }
  }
}

template <typename T>
Job<T>::Job(Solver<T> &solver, const T mindur, const T maxdur)
    : start(solver.newTemporal()),
      end((mindur == maxdur ? TemporalVar(start.id(), mindur)
                            : solver.newTemporal())),
      optional(Constant::NoVarx) {
  min_duration = mindur;
  max_duration = maxdur;

  if (start.id() != end.id()) {
    solver.set(start.before(end, min_duration));
    solver.set(end.before(start, -max_duration));
  }
}

template <typename T> int Job<T>::id() const { return start.id(); }

template<typename T>
bool Job<T>::operator==(const Job<T>& t) const {
  return id() == t.id();
}

template <typename T> T Job<T>::getEarliestStart(Solver<T> &solver) const {
  return start.earliest(solver);
}

template <typename T> T Job<T>::getLatestStart(Solver<T> &solver) const {
  return start.latest(solver);
}

template <typename T> T Job<T>::getEarliestEnd(Solver<T> &solver) const {
  return end.earliest(solver);
}

template <typename T> T Job<T>::getLatestEnd(Solver<T> &solver) const {
  return end.latest(solver);
}

template <typename T> bool Job<T>::mustExist(Solver<T> &) const { return true; }

template <typename T> bool Job<T>::cannotExist(Solver<T> &) const {
  return false;
}

template <typename T> T Job<T>::minDuration() const { return min_duration; }

template <typename T> T Job<T>::maxDuration() const { return max_duration; }

template <typename T> event Job<T>::getStart() const { return start.id(); }

template <typename T> event Job<T>::getEnd() const { return end.id(); }

template <typename T> ostream &Job<T>::display(ostream &os) const {
  os << "t" << id(); //<< ": [" << start.earliest(solver) << ".." <<
                     // end.latest(solver) << "]";
  return os;
}

template <typename T> ostream &operator<<(ostream &os, const Job<T> &x) {
  return x.display(os);
}
}

#endif // __MODEL_HPP
