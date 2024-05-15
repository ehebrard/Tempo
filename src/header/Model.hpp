
#ifndef __TEMPO_MODEL_HPP
#define __TEMPO_MODEL_HPP


#include "Literal.hpp"


namespace tempo {

template<typename T> class Solver;

/**********************************************
 * Job
 **********************************************/

template<typename T>
class NumericVar {

public:
  NumericVar(const var_t i) : _id_(i) {}

  NumericVar(NumericVar<T> &) = default;
  NumericVar(NumericVar<T> &&) = default;
  NumericVar<T> &operator=(NumericVar<T> &) = default;
  NumericVar<T> &operator=(NumericVar<T> &&) = default;

  T min(Solver<T> &sc) const;
  T max(Solver<T> &sc) const;

  Literal<T> operator>=(const T t) const;
  Literal<T> operator>(const T t) const;
  Literal<T> operator<=(const T t) const;
  Literal<T> operator<(const T t) const;

  var_t id() const { return _id_; }

  operator var_t() const { return _id_; }

  std::ostream &display(std::ostream &os) const;

protected:
  const var_t _id_;
};


template<typename T>
class BooleanVar {

public:
  BooleanVar(const var_t i) : _id_(i) {}

  BooleanVar(BooleanVar<T> &) = default;
  BooleanVar(BooleanVar<T> &&) = default;
  BooleanVar<T> &operator=(BooleanVar<T> &) = default;
  BooleanVar<T> &operator=(BooleanVar<T> &&) = default;

  Literal<T> operator==(const bool t) const;

  var_t id() const { return _id_; }

  operator var_t() const { return _id_; }

  std::ostream &display(std::ostream &os) const;

protected:
  const var_t _id_;
};

template<typename T>
Literal<T> BooleanVar<T>::operator==(const bool t) const {
  return Literal<T>(t, _id_, Constant::NoSemantic);
}


template<typename T>
class DisjunctVar : public BooleanVar<T> {

public:
  DisjunctVar(const var_t i, const info_t d) : BooleanVar<T>(i), _edge_id_(d) {}

  DisjunctVar(DisjunctVar<T> &) = default;
  DisjunctVar(DisjunctVar<T> &&) = default;
  DisjunctVar<T> &operator=(DisjunctVar<T> &) = default;
  DisjunctVar<T> &operator=(DisjunctVar<T> &&) = default;

  Literal<T> operator==(const bool t) const;

  std::ostream &display(std::ostream &os) const;

private:
  const info_t _edge_id_;
};

template<typename T>
Literal<T> DisjunctVar<T>::operator==(const bool t) const {
  return Literal<T>(t, BooleanVar<T>::_id_, _edge_id_ + t);
}

template<typename T>
class TemporalVar : public NumericVar<T> {

public:
  TemporalVar(const var_t i, const T o = 0) : NumericVar<T>(i), _offset(o) {}

  TemporalVar(TemporalVar<T> &) = default;
  TemporalVar(TemporalVar<T> &&) = default;
  TemporalVar<T> &operator=(TemporalVar<T> &) = default;
  TemporalVar<T> &operator=(TemporalVar<T> &&) = default;

  T earliest(Solver<T> &) const;
  T latest(Solver<T> &) const;

  Literal<T> after(const T t) const;
  Literal<T> before(const T t) const;

  DistanceConstraint<T> after(const TemporalVar<T> &e, const T t = 0) const;
  DistanceConstraint<T> before(const TemporalVar<T> &e, const T t = 0) const;

  T offset() const { return _offset; }

  std::ostream &display(std::ostream &os) const;

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




//{START(job[i]), END(job[j]), -job.getTransitionTime(j, i)},
// s_i -> e_j
// e_j - s_i <= - t


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




template<typename T>
class Job {
public:
//  Job(Solver<T> &s);
  Job(Solver<T> &s, const T mindur, const T maxdur);

  Job(Job<T> &) = default;
  Job(Job<T> &&) = default;
  Job<T> &operator=(Job<T> &) = default;
  Job<T> &operator=(Job<T> &&) = default;

  T getEarliestStart() const;
  T getLatestStart() const;
  T getEarliestEnd() const;
  T getLatestEnd() const;

  bool mustExist() const;
  bool cannotExist() const;

  T minDuration() const;
  T maxDuration() const;

  event getStart() const;
  event getEnd() const;

  int id() const;
  bool operator==(const Job<T> &t) const;

  std::ostream &display(std::ostream &os) const;

private:
//    int _id_{-1};

  Solver<T> &solver;

public:
  TemporalVar<T> start;
  TemporalVar<T> end;

private:
  T min_duration;
  T max_duration;
  BooleanVar<T> optional;
};


template<typename T>
class DisjunctiveResource : public std::vector<Job<T>> {
public:
  template <typename Iter>
  DisjunctiveResource(Solver<T> &s, Iter beg_job, Iter end_job);

  DisjunctiveResource(Solver<T> &s, std::vector<Job<T>> &&container);

  template <typename C> void createOrderVariables(C &container);

private:
  Solver<T> &solver;
};

template<typename T>
template<typename Iter>
DisjunctiveResource<T>::DisjunctiveResource(Solver<T> &s, Iter beg_job, Iter end_job) : std::vector<Job<T>>(beg_job, end_job), solver(s) {}

template<typename T>
DisjunctiveResource<T>::DisjunctiveResource(Solver<T> &s, std::vector<Job<T>>&& container) : std::vector<Job<T>>(container), solver(s) {}

template<typename T>
template<typename C>
void DisjunctiveResource<T>::createOrderVariables(C& container) {
  for (auto a{this->begin()}; a != this->end(); ++a) {
    for (auto b{a + 1}; b != this->end(); ++b) {
      //            container.insert(container.end(),
      //            solver.newDisjunct(a->end.before(b->start),
      //            b->end.before(a->start)));
      container.push_back(
          solver.newDisjunct(a->end.before(b->start), b->end.before(a->start)));
    }
  }
}

template <typename T>
Job<T>::Job(Solver<T> &s, const T mindur, const T maxdur)
    : solver(s), start(solver.newTemporal()),
      end((mindur == maxdur ? TemporalVar(start.id(), mindur)
                            : solver.newTemporal())),
      // end(solver.newTemporalVar()),
      optional(Constant::NoVarx) {
  min_duration = mindur;
  max_duration = maxdur;

  //    std::cout << mindur << ".." << maxdur << std::endl;

  if (start.id() != end.id()) {
    //      std::cout << "\nmindur: " << start.before(end, min_duration) <<
    //      std::endl;
    solver.set(start.before(end, min_duration));

    //      std::cout << "\nmaxdur: " << end.before(start, -max_duration) <<
    //      std::endl;
    solver.set(end.before(start, -max_duration));
  }
}

template <typename T> int Job<T>::id() const { return start.id(); }

template<typename T>
bool Job<T>::operator==(const Job<T>& t) const {
  return id() == t.id();
}

template<typename T>
T Job<T>::getEarliestStart() const {
  return start.earliest(solver);
}

template<typename T>
T Job<T>::getLatestStart() const {
  return start.latest(solver);
}

template<typename T>
T Job<T>::getEarliestEnd() const {
  return end.earliest(solver);
}

template<typename T>
T Job<T>::getLatestEnd() const {
  return end.latest(solver);
}

template <typename T> bool Job<T>::mustExist() const { return true; }

template <typename T> bool Job<T>::cannotExist() const { return false; }

template <typename T> T Job<T>::minDuration() const { return min_duration; }

template <typename T> T Job<T>::maxDuration() const { return max_duration; }

template <typename T> event Job<T>::getStart() const { return start.id(); }

template <typename T> event Job<T>::getEnd() const { return end.id(); }

template<typename T>
std::ostream &Job<T>::display(std::ostream &os) const {
  os << "t" << id() << ": [" << start.earliest(solver) << ".."
     << end.latest(solver) << "]";
  return os;
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const Job<T> &x) {
  return x.display(os);
}



}

#endif // __MODEL_HPP
