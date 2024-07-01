/************************************************
 * Tempo Model.hpp
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

#ifndef __TEMPO_MODEL_HPP
#define __TEMPO_MODEL_HPP


#include "Literal.hpp"
#include "constraints/Cardinality.hpp"

using namespace std;

namespace tempo {

template<typename T> class Solver;

template <typename T> class BooleanExpression;

//! Wrapper/pointer for numeric variables
/*!
Stores the id of the actual numeric variable, and implements various helper methods
 */
template<typename T=int>
class NumericVar {

public:
    NumericVar() {};
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
    var_t _id_{Constant::NoIndex};
};

//! Wrapper/pointer for Boolean variables
/*!
Stores the id of the actual Boolean variable, and implements various helper methods
 */
template<typename T=int>
class BooleanVar {

public:
    BooleanVar() {}
  BooleanVar(const var_t i) : _id_(i) {}

  Literal<T> operator==(const bool t) const;

  BooleanExpression<T> implies(const BooleanVar<T> x) const;
  BooleanExpression<T> implies(const BooleanExpression<T> x) const;

  var_t id() const { return _id_; }

  operator var_t() const { return _id_; }

  std::ostream &display(std::ostream &os) const;

  static bool isNumeric() { return false; }

protected:
    var_t _id_{Constant::NoIndex};
};

template<typename T>
Literal<T> BooleanVar<T>::operator==(const bool t) const {
  return makeBooleanLiteral<T>(t, _id_, Constant::NoSemantic);
}

//! Wrapper/pointer for disjunct variables
/*!
Stores the id of the actual Boolean variable with difference logic semantic, and implements various helper methods
 */
template<typename T=int>
class DisjunctVar : public BooleanVar<T> {

public:
    DisjunctVar() {}
  DisjunctVar(const var_t i, const info_t d) : BooleanVar<T>(i), _edge_id_(d) {}

  Literal<T> operator==(const bool t) const;

  std::ostream &display(std::ostream &os) const;

  static bool isNumeric() { return false; }

private:
    info_t _edge_id_{Constant::NoIndex};
};

template<typename T>
Literal<T> DisjunctVar<T>::operator==(const bool t) const {
  return makeBooleanLiteral<T>(t, BooleanVar<T>::_id_, _edge_id_ + t);
}

//! Wrapper/pointer for temporal variables
/*!
Stores the id of the actual numeric variable, an offset (the variable encode for time t=x+offset where x is the numeric variable) and implements various helper methods
 */
template<typename T=int>
class TemporalVar : public NumericVar<T> {

public:
    TemporalVar() {}
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
    T _offset{0};
};


//! Wrapper for interval variables
/*!
Stores
 - the id of a  temporal variable standing for the start
 - the id of a  temporal variable standing for the end
 - the id of a  Boolean variable standing whether the interval actually is in the schedule
 */
template <typename T = int> class Interval {
public:
  Interval() {}
  Interval(Solver<T> &s, const T mindur = 0,
           const T maxdur = Constant::Infinity<T>,
           const BooleanVar<T> opt = Constant::NoVar);

  //    Interval(const Interval<T>&) = default;

  T getEarliestStart(Solver<T> &s) const;
  T getLatestStart(Solver<T> &s) const;
  T getEarliestEnd(Solver<T> &s) const;
  T getLatestEnd(Solver<T> &s) const;

  bool mustExist(Solver<T> &s) const;
  bool cannotExist(Solver<T> &s) const;

  T minDuration() const;
  T maxDuration() const;

  var_t getStart() const;
  var_t getEnd() const;

  int id() const;
  bool operator==(const Interval<T> &t) const;

  std::ostream &display(std::ostream &os) const;

  TemporalVar<T> start;
  TemporalVar<T> end;

  bool isOptional() const { return exist.id() != Constant::NoVar; }

  BooleanVar<T> exist{Constant::NoVar};

private:
  T min_duration{0};
  T max_duration{Constant::Infinity<T>};
};


//! Wrapper for disjunctive resources
/*!
Contains the list of interval requiring this resource
 */
template <typename T = int>
class DisjunctiveResource : public vector<Interval<T>> {
public:
  using vector<Interval<T>>::vector;
    
    // create and insert in 'disjuncts' the set of disjunctive variables necessary to ensure that this constraint is satisfied
  template <typename C>
  void createOrderVariables(Solver<T> &solver, C &disjuncts);
    // create and insert in 'disjuncts' the set of disjunctive variables necessary to ensure that this constraint is satisfied, according to the options
    template <typename C>
    void createOptionalOrderVariables(Solver<T> &solver, C &disjuncts, C &options);
};


template <typename T = int>
class SchedulingModel : public std::vector<Interval<T>> {
public:
    std::vector<DisjunctiveResource<T>> resources;
    std::vector<Literal<T>> precedences;
};

//! "This expression is not a constraint" exception
class ModelingException : public std::exception {
public:
  ModelingException(std::string msg) : msg(msg) {}

  virtual const char *what() const throw() { return msg.c_str(); }

private:
  std::string msg;
};

////! Variable type exception
// class VarTypeException : public std::exception {
// public:
//
//     VarTypeException()  {}
//
//   virtual const char *what() const throw() { return "Wrong variable type"; }
// };
//
////! "This expression is not a constraint" exception
// class RootException : public std::exception {
// public:
//
//     RootException()  {}
//
//   virtual const char *what() const throw() { return "This expression is not a
//   constraint"; }
// };

template <typename T>
class ExpressionImpl;

template <typename T = int>
class Expression {
    
public:
    Expression() {}
    Expression(ExpressionImpl<T> *i) : impl(i) {}
  
protected:
    ExpressionImpl<T> *impl{NULL};
};


template <typename T = int>
class BooleanExpression : public BooleanVar<T>, public Expression<T> {
    
public:
    BooleanExpression(ExpressionImpl<T> *i) : Expression<T>(i) {}
    BooleanExpression(BooleanVar<T> x) {
        BooleanVar<T>::_id_ = x.id();
    }
    
    void extract(Solver<T>& solver) {
        if(Expression<T>::impl)
            BooleanVar<T>::_id_ = Expression<T>::impl->extract(solver);
    }
    
    void post(Solver<T>& solver) {
        if(Expression<T>::impl)
            Expression<T>::impl->post(solver);
    }

};

template <typename T = int>
class NumericExpression : public NumericVar<T>, public Expression<T> {
    
public:
    NumericExpression(ExpressionImpl<T> *i) : Expression<T>(i) {}
    NumericExpression(NumericVar<T> x) {
        NumericVar<T>::_id_ = x.id();
    }
    
    void extract(Solver<T>& solver) {
        if(Expression<T>::impl)
            NumericVar<T>::_id_ = Expression<T>::impl->extract(solver);
    }
};

template <typename T = int>
class ExpressionImpl {
public:
    virtual var_t extract(Solver<T>& ) = 0;
    virtual void post(Solver<T> &) {
      throw ModelingException("This predicate cannot be a constraint");
    }

//    BooleanVar<T> getBoolean() { throw VarTypeException(); }
//    NumericVar<T> getNumeric() { throw VarTypeException(); }

};


template <typename T = int>
class LogicalAndExpression : public ExpressionImpl<T> {
public:
    
    template <typename Iter>
    LogicalAndExpression(Iter beg_var, Iter end_var) {
        for(auto x{beg_var}; x!=end_var; ++x) {
          boolean_arguments.emplace_back(*x);
        }
    }

    var_t extract(Solver<T>& solver) override {
        std::vector<Literal<T>> L;
        for(auto x : boolean_arguments) {
            x.extract(solver);
            L.push_back(x == false);
        }
        self = solver.newBoolean();
        L.push_back(self == true);
        solver.clauses.add(L.begin(), L.end());
        for(auto x : boolean_arguments) {
            L.clear();
            L.push_back(x == true);
            L.push_back(self == false);
            solver.clauses.add(L.begin(), L.end());
        }
        return self.id();
    }
    
    void post(Solver<T>& solver) override {
        for(auto x : boolean_arguments) {
            x.extract(solver);
            solver.set(x == true);
        }
    }
    
private:
    BooleanVar<T> self;
    std::vector<BooleanExpression<T>> boolean_arguments;
};


template<typename T>
BooleanExpression<T> operator&&(BooleanVar<T>& x, BooleanVar<T>& y) {
  std::vector<BooleanVar<T>> sc{x, y};
  BooleanExpression<T> exp(new LogicalAndExpression<T>(sc.begin(), sc.end()));
  return exp;
}

template <typename T, typename Iterable>
BooleanExpression<T> BigAnd(Iterable &X) {
  BooleanExpression<T> exp(new LogicalAndExpression(X.begin(), X.end()));
  return exp;
}

template <typename T = int>
class LogicalOrExpression : public ExpressionImpl<T> {
public:
  template <typename Iter> LogicalOrExpression(Iter beg_var, Iter end_var) {
    for (auto x{beg_var}; x != end_var; ++x) {
      boolean_arguments.emplace_back(*x);
    }
  }

  var_t extract(Solver<T> &solver) override {
    std::vector<Literal<T>> L;
    for (auto x : boolean_arguments) {
      x.extract(solver);
      L.push_back(x == true);
    }
    self = solver.newBoolean();
    L.push_back(self == false);
    solver.clauses.add(L.begin(), L.end());
    for (auto x : boolean_arguments) {
      L.clear();
      L.push_back(x == false);
      L.push_back(self == true);
      solver.clauses.add(L.begin(), L.end());
    }
    return self.id();
  }

  void post(Solver<T> &solver) override {
    std::vector<Literal<T>> L;
    for (auto x : boolean_arguments) {
      x.extract(solver);
      L.push_back(x == true);
    }
    solver.clauses.add(L.begin(), L.end());
  }

private:
  BooleanVar<T> self;
  std::vector<BooleanExpression<T>> boolean_arguments;
};

template <typename T>
BooleanExpression<T> operator||(BooleanVar<T> &x, BooleanVar<T> &y) {
  std::vector<BooleanVar<T>> sc{x, y};
  BooleanExpression<T> exp(new LogicalOrExpression<T>(sc.begin(), sc.end()));
  return exp;
}

template <typename T, typename Iterable>
BooleanExpression<T> BigOr(Iterable &X) {
  BooleanExpression<T> exp(new LogicalOrExpression(X.begin(), X.end()));
  return exp;
}

template <typename T = int>
class LogicalImplicationExpression : public ExpressionImpl<T> {
public:
  LogicalImplicationExpression(BooleanExpression<T> x, BooleanExpression<T> y)
      : implicant(x), implied(y) {}
  LogicalImplicationExpression(BooleanVar<T> x, BooleanVar<T> y)
      : implicant(x), implied(y) {}

  var_t extract(Solver<T> &solver) override {
    implicant.extract();
    implied.extract();
    self = solver.newBoolean();

    std::vector<Literal<T>> cl{implicant == false, implied == true,
                               self = false};
    solver.clauses.add(cl.begin(), cl.end());

    cl = {self == true, implicant == true};
    solver.clauses.add(cl.begin(), cl.end());

    cl = {self == true, implied == false};
    solver.clauses.add(cl.begin(), cl.end());

    return self.id();
  }

  void post(Solver<T> &solver) override {
    implicant.extract();
    implied.extract();
    std::vector<Literal<T>> cl{implicant == false, implied == true};
    solver.clauses.add(cl.begin(), cl.end());
  }

private:
  BooleanVar<T> self;
  BooleanExpression<T> implicant;
  BooleanExpression<T> implied;
};

template <typename T>
BooleanExpression<T> BooleanVar<T>::implies(const BooleanVar<T> x) const {
  BooleanExpression<T> exp(new LogicalImplicationExpression<T>(*this, x));
  return exp;
}

template <typename T>
BooleanExpression<T>
BooleanVar<T>::implies(const BooleanExpression<T> x) const {
  BooleanExpression<T> exp(
      new LogicalImplicationExpression<T>(BooleanExpression<T>(*this), x));
  return exp;
}

// template <typename T = int>
// class LogicalNotExpression : public ExpressionImpl<T> {
// public:
//
//     LogicalNotExpression(BooleanExpression<T> x) argument(x) {}
//     LogicalNotExpression(BooleanVar<T> x) argument(x) {}
//
//     var_t extract(Solver<T>& solver) override {
//         argument.extract();
//         self = solver.newBoolean();
//
//         std::vector<Literal<T>> L{self==false, argument==true};
//         solver.clauses.add(L.begin(), L.end());
//
//         L = {self==true, argument==false};
//         solver.clauses.add(L.begin(), L.end());
//     }
//
//     void post(Solver<T>& solver) override {
//         argument.extract();
//         solver.set(argument == false);
//     }
//
// private:
//     BooleanVar<T> self;
//     BooleanExpression<T> argument;
// };
//
//
// template<typename T>
// BooleanExpression<T> operator~(BooleanExpression<T>& x) {
//     BooleanExpression<T> exp(new LogicalNotExpression<T>(sc.begin(),
//     sc.end())); return exp;
// }
//
//
// template<typename T, typename Iterable>
// BooleanExpression<T> BigOr(Iterable& X) {
//     BooleanExpression<T> exp(new LogicalOrExpression(X.begin(), X.end()));
//     return exp;
// }

template <typename T = int>
class CardinalityExpression : public ExpressionImpl<T> {
public:
    
    template <typename Iter>
    CardinalityExpression(Iter beg_var, Iter end_var, const T l=0, const T u=Constant::Infinity<T>) : lb(l), ub(u) {
        for(auto x{beg_var}; x!=end_var; ++x) {
          boolean_arguments.emplace_back(*x);
        }
    }
    
    var_t extract(Solver<T>& solver) override {
        std::vector<Literal<T>> L;
        for(auto x : boolean_arguments) {
            x.extract(solver);
            L.push_back(x == true);
        }
        self = solver.newNumeric();
        solver.post( new CardinalityLeqVar<T>>(solver, L.begin(), L.end(), self.id()) );
        solver.post( new CardinalityGeqVar<T>>(solver, L.begin(), L.end(), self.id()) );

        return self.id();
    }
    
    void post(Solver<T>& solver) override {
        std::vector<Literal<T>> L;
        for(auto x : boolean_arguments) {
            x.extract(solver);
            L.push_back(x == true);
        }
        T n{static_cast<T>(L.size())};
        if(ub < n) {
            solver.post( new CardinalityConst<T>>(solver, L.begin(), L.end(), ub) );
        }
        if(lb > 0) {
            for(auto &l : L) l = ~l;
            solver.post( new CardinalityConst<T>>(solver, L.begin(), L.end(), n-lb) );
        }
    }
    
private:
    NumericVar<T> self;
    T lb{0};
    T ub{Constant::Infinity<T>};
    std::vector<BooleanExpression<T>> boolean_arguments;
};

template <typename T, typename Iterable>
NumericExpression<T> Cardinality(Iterable &X) {
  NumericExpression<T> exp(new CardinalityExpression(X.begin(), X.end()));
  return exp;
}

template <typename T = int>
class NoOverlapExpressionImpl : public ExpressionImpl<T>,
                                public std::vector<Interval<T>> {
public:
  //    template <typename Iter>
  //    NoOverlapExpression(Iter beg_int, Iter end_int) :
  //    std::vector<Interval<T>>(beg_int, end_int) {
  ////        for(auto i{beg_int}; i!=end_int; ++i) {
  ////            this->push_back(*i);
  ////        }
  //    }

  NoOverlapExpressionImpl(Interval<T> &sched) : schedule(sched) {}

  var_t extract(Solver<T> &) override {
    throw ModelingException("NoOverlap is not a predicate");
    return Constant::NoVar;
  }

  void post(Solver<T> &solver) override {
    for (auto a{this->begin()}; a != this->end(); ++a) {
      for (auto b{a + 1}; b != this->end(); ++b) {
        if (a->isOptional() and b->isOptional()) {
          BooleanVar<T> x = solver.newBoolean();
          std::vector<Literal<T>> cl{x == false, a->exist == true};
          solver.clauses.add(cl.begin(), cl.end());

          cl = {x == false, b->exist == true};
          solver.clauses.add(cl.begin(), cl.end());

          cl = {x == true, a->exist == false, b->exist == false};
          solver.clauses.add(cl.begin(), cl.end());

          disjunct.push_back(solver.newDisjunct(a->end.before(b->start),
                                                b->end.before(a->start), x));

          relevant.push_back(x);
        } else if (a->isOptional()) {
          disjunct.push_back(solver.newDisjunct(
              a->end.before(b->start), b->end.before(a->start), a->exist));
          //                relevant.push_back(a->exist);
        } else if (b->isOptional()) {
          disjunct.push_back(solver.newDisjunct(
              a->end.before(b->start), b->end.before(a->start), b->exist));
          //                relevant.push_back(b->exist);
        } else {
          disjunct.push_back(solver.newDisjunct(a->end.before(b->start),
                                                b->end.before(a->start)));
        }

        solver.addToSearch(disjunct.back());
      }
    }

    if (solver.getOptions().edge_finding)
      solver.postEdgeFinding(schedule, this->begin(), this->end(),
                             this->begDisjunct());

    if (solver.getOptions().transitivity)
      solver.postTransitivity(schedule, this->begin(), this->end(),
                              this->begDisjunct());
  }

  std::vector<DisjunctVar<T>>::iterator begDisjunct() {
    return disjunct.begin();
  }
  std::vector<DisjunctVar<T>>::iterator endDisjunct() { return disjunct.end(); }

private:
  Interval<T> schedule;
  std::vector<DisjunctVar<T>> disjunct;
  std::vector<BooleanVar<T>> relevant;
};

template <typename T = int>
class NoOverlapExpression : public BooleanExpression<T> {

public:
  NoOverlapExpression(NoOverlapExpressionImpl<T> *i)
      : BooleanExpression<T>(i) {}
    
    std::vector<Interval<T>>::iterator begin() {
      return static_cast<NoOverlapExpressionImpl<T> *>(BooleanExpression<T>::impl)
          ->begin();
    }
    std::vector<Interval<T>>::iterator end() {
      return static_cast<NoOverlapExpressionImpl<T> *>(BooleanExpression<T>::impl)
          ->end();
    }

  std::vector<DisjunctVar<T>>::iterator begDisjunct() {
    return static_cast<NoOverlapExpressionImpl<T> *>(BooleanExpression<T>::impl)
        ->begDisjunct();
  }
  std::vector<DisjunctVar<T>>::iterator endDisjunct() {
    return static_cast<NoOverlapExpressionImpl<T> *>(BooleanExpression<T>::impl)
        ->endDisjunct();
  }
    
    void push_back(const Interval<T>& i) { static_cast<NoOverlapExpressionImpl<T> *>(BooleanExpression<T>::impl)
        ->push_back(i); }
};

template <typename T, typename Iterable>
NoOverlapExpression<T> NoOverlap(Interval<T> &schedule, Iterable &X) {
  auto impl{new NoOverlapExpressionImpl<T>(schedule)};
  for (auto x : X)
    impl->push_back(x);
  NoOverlapExpression<T> exp(impl);
  return exp;
}

template <typename T>
NoOverlapExpression<T> NoOverlap(Interval<T> &schedule) {
  auto impl{new NoOverlapExpressionImpl<T>(schedule)};
  NoOverlapExpression<T> exp(impl);
  return exp;
}

//template <typename T = int>
//class LogicalAndExpression {
//public:
//    void extract(Solver<T>& solver) {
//        std::vector<Literal<T>> L;
//        for(auto x : arguments) {
//            x.extract(solver);
//            L.push_back(x.getBoolean() == false);
//        }
//        self = solver.newBoolean();
//        L.push_back(self == true);
//        solver.clauses.add(L.begin(), L.end());
//        for(auto x : arguments) {
//            L.clear();
//            L.push_back(x.getBoolean() == true);
//            L.push_back(self == false);
//            solver.clauses.add(L.begin(), L.end());
//        }
//    }
//    
//private:
//    BooleanVar<T> self;
//    std::vector<Expression<T>> arguments;
//}





/*!
 NumericVar  implementation
*/
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


/*!
 TemporalVar  implementation
*/
template<typename T>
T TemporalVar<T>::earliest(Solver<T>& s) const {
    auto v{NumericVar<T>::min(s)};
    if(v == -Constant::Infinity<T>)
        return v;
    return v + _offset;
    
//  return NumericVar<T>::min(s) + _offset;
}

template<typename T>
T TemporalVar<T>::latest(Solver<T>& s) const {
//    auto r{NumericVar<T>::max(s)};
//    std::cout << "\nlatest:" << r << " + " << _offset << std::endl;
//    return r + _offset;
    
//  return NumericVar<T>::max(s) + _offset;
    auto v{NumericVar<T>::max(s)};
    if(v == Constant::Infinity<T>)
        return v;
    return v + _offset;
}

template<typename T>
Literal<T> TemporalVar<T>::after(const T t) const {
  return geq<T>(NumericVar<T>::_id_, (t == Constant::Infinity<T> ? t : t - _offset));
}

template<typename T>
Literal<T> TemporalVar<T>::before(const T t) const {
  return leq<T>(NumericVar<T>::_id_, (t == Constant::Infinity<T> ? t : t - _offset));
}

template<typename T>
DistanceConstraint<T> TemporalVar<T>::after(const TemporalVar<T>& e, const T t) const {
  return e.before(*this, t);
}

template<typename T>
DistanceConstraint<T> TemporalVar<T>::before(const TemporalVar<T>& e, const T t) const {
  return {e.id(), NumericVar<T>::_id_, (t == Constant::Infinity<T> ? t : e.offset() - _offset - t)};
}


/*!
 BooleanVar  implementation
*/
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


/*!
 DisjunctiveResource  implementation
*/
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

/*!
 DisjunctiveResource  implementation
*/
template <typename T>
template <typename C>
void DisjunctiveResource<T>::createOptionalOrderVariables(Solver<T> &solver,
                                                          C &disjuncts,
                                                          C &options) {
  for (auto a{this->begin()}; a != this->end(); ++a) {
    for (auto b{a + 1}; b != this->end(); ++b) {
      //                  container.insert(container.end(),
      //                  solver.newDisjunct(a->end.before(b->start),
      //                  b->end.before(a->start)));
      options.push_back(solver.newBoolean());
      disjuncts.push_back(solver.newDisjunct(options.back().id(),
                                            a->end.before(b->start),
                                            b->end.before(a->start)));
    }
  }
}

/*!
 Interval  implementation
*/
template <typename T>
Interval<T>::Interval(Solver<T> &solver, const T mindur, const T maxdur,
                      const BooleanVar<T> opt)
    : start(solver.newTemporal()),
      end((mindur == maxdur ? TemporalVar(start.id(), mindur)
                            : solver.newTemporal())),
      exist(opt) {
  min_duration = mindur;
  max_duration = maxdur;

  if (start.id() != end.id()) {
    solver.set(start.before(end, min_duration));
    solver.set(end.before(start, -max_duration));
  }
}

template <typename T> int Interval<T>::id() const { return start.id(); }

template <typename T> bool Interval<T>::operator==(const Interval<T> &t) const {
  return id() == t.id();
}

template <typename T> T Interval<T>::getEarliestStart(Solver<T> &solver) const {
  return start.earliest(solver);
}

template <typename T> T Interval<T>::getLatestStart(Solver<T> &solver) const {
  return start.latest(solver);
}

template <typename T> T Interval<T>::getEarliestEnd(Solver<T> &solver) const {
  return end.earliest(solver);
}

template <typename T> T Interval<T>::getLatestEnd(Solver<T> &solver) const {
  return end.latest(solver);
}

template <typename T> bool Interval<T>::mustExist(Solver<T> &) const {
  return true;
}

template <typename T> bool Interval<T>::cannotExist(Solver<T> &) const {
  return false;
}

template <typename T> T Interval<T>::minDuration() const {
  return min_duration;
}

template <typename T> T Interval<T>::maxDuration() const {
  return max_duration;
}

template <typename T> var_t Interval<T>::getStart() const { return start.id(); }

template <typename T> var_t Interval<T>::getEnd() const { return end.id(); }

template <typename T> ostream &Interval<T>::display(ostream &os) const {
  os << "t" << id(); //<< ": [" << start.earliest(solver) << ".." <<
                     // end.latest(solver) << "]";
  return os;
}

template <typename T> ostream &operator<<(ostream &os, const Interval<T> &x) {
  return x.display(os);
}
}

#endif // __MODEL_HPP
