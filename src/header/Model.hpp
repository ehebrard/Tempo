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
#include "constraints/PseudoBoolean.hpp"
#include "constraints/SumConstraint.hpp"

using namespace std;

/// @TODO rewritte:
/// ExpressionImpl should be an interface with var_id(); post(solver);
/// extract(solver); BooleanVar and NumericVar should implement Expression
///  BooleanExpression and NumericExpression are pointers to Expression
///  Relational expressions are ExpressionImpl with typically a
///  vector<BooleanExpression> or  a vector<NumericExpression> or both

namespace tempo {

template<typename T> class Solver;

//template <typename T> class BooleanExpression;

//! "This expression is not a constraint" exception
class ModelingException : public std::exception {
public:
  ModelingException(std::string msg) : msg(msg) {}

  virtual const char *what() const throw() { return msg.c_str(); }

private:
  std::string msg;
};

//! Interface for expression implementations
template <typename T = int> class ExpressionImpl {
public:
    
    virtual ~ExpressionImpl() { /*std::cout << "delete expr\n";*/ }
    
  virtual var_t extract(Solver<T> &) = 0;
  virtual void post(Solver<T> &) {
    throw ModelingException("This predicate cannot be a constraint");
  }

  virtual string name() const { return "some expression"; }
  virtual var_t id() const { return Constant::NoVar; }
  virtual T offset() const { return 0; }
  virtual index_t semantic() const { return Constant::NoSemantic; }
};

//! Wrapper/pointer for numeric variables and expressions
/*!
Stores the id of the actual numeric variable, and implements various helper
methods
 */

struct ExpressionFlag {
  ExpressionFlag(const bool t = false) : _is_expression(t) {}
  bool _is_expression;
};

template <typename T = int> struct NumInfo {
  NumInfo() {}
  NumInfo(const var_t i, const T o) : _id_(i), _offset(o) {}

  var_t _id_{Constant::NoIndex};
  T _offset{0};
};

template <typename T = int> class NumericVar : public ExpressionFlag {

public:
    NumericVar() {};
    NumericVar(ExpressionImpl<T> *i) : ExpressionFlag(true), implem(i) {}
    NumericVar(const var_t i, const T o = 0)
        : ExpressionFlag(false), data(i, o) {}

    T min(Solver<T> &sc) const;
    T max(Solver<T> &sc) const;

    T earliest(Solver<T> &) const;
    T latest(Solver<T> &) const;

    Literal<T> after(const T t) const;
    Literal<T> before(const T t) const;

    DistanceConstraint<T> after(const NumericVar<T> &e, const T t = 0) const;
    DistanceConstraint<T> before(const NumericVar<T> &e, const T t = 0) const;

    var_t id() const { return (ExpressionFlag::_is_expression ? implem->id() : data._id_); }

    std::ostream &display(std::ostream &os) const;

    static bool isNumeric() { return true; }

    void setId(const var_t i) { data._id_ = i; }

    T offset() const { return data._offset; }

    void setOffset(const T o) { data._offset = o; }

    void extract(Solver<T> &solver) {
      if (ExpressionFlag::_is_expression) {
        implem->extract(solver);
        data = NumInfo<T>(implem->id(), implem->offset());
        ExpressionFlag::_is_expression = false;
        solver.trash_bin.push_back(implem);
//        std::cout << "extract " << implem->name() << std::endl;
      }
    }

    void post(Solver<T> &solver) {
      if (ExpressionFlag::_is_expression) {
        throw ModelingException("Numeric expression cannot be constraints");
      } else {
        solver.addToSearch(*this);
      }
    }

  protected:
    union {
      struct NumInfo<T> data;
      ExpressionImpl<T> *implem;
    };
};

template <typename T = int> struct BoolInfo {
  BoolInfo() {}
  BoolInfo(const var_t i, const info_t f) : _id_(i), _edge_id_(f) {}

  var_t _id_{Constant::NoIndex};
  index_t _edge_id_{Constant::NoIndex};
};

//! Wrapper/pointer for Boolean variables and expressions
/*!
Stores the id of the actual Boolean variable, and implements various helper
methods
 */
template <typename T = int> class BooleanVar : public ExpressionFlag {

public:
  BooleanVar() { data._id_ = Constant::NoIndex; }
  BooleanVar(ExpressionImpl<T> *i) : ExpressionFlag(true), implem(i) {}
  BooleanVar(const var_t i, const index_t f = 0)
      : ExpressionFlag(false), data(i, f) {}

  Literal<T> operator==(const bool t) const;

  BooleanVar<T> implies(const BooleanVar<T> x) const;

var_t id() const { return (ExpressionFlag::_is_expression ? implem->id() : data._id_); }
    
    index_t semantic() const { return (ExpressionFlag::_is_expression ? implem->semantic() : data._edge_id_); }

  operator var_t() const { return id(); }

  std::ostream &display(std::ostream &os) const;

  static bool isNumeric() { return false; }

  void setId(const var_t i) { data._id_ = i; }

  void extract(Solver<T> &solver) {
    if (ExpressionFlag::_is_expression) {
      implem->extract(solver);
      data = BoolInfo<T>(implem->id(), implem->semantic());
      ExpressionFlag::_is_expression = false;
      solver.trash_bin.push_back(implem);
//      std::cout << "extract (b) " << implem->name() << std::endl;
    }
  }

  void post(Solver<T> &solver) {
    if (ExpressionFlag::_is_expression) {
      implem->post(solver);
      ExpressionFlag::_is_expression = false;
      solver.trash_bin.push_back(implem);
//      std::cout << "post (b) " << implem->name() << std::endl;
    } else {
      solver.addToSearch(*this);
    }
  }

    
//    Literal<T> lit(const bool sign=true) {
//        return makeBooleanLiteral<T>(sign, id(), semantic()+sign);
//    }
    
protected:
  union {
    struct BoolInfo<T> data;
    ExpressionImpl<T> *implem;
  };
};

template <typename T> Literal<T> BooleanVar<T>::operator==(const bool t) const {
  return makeBooleanLiteral<T>(t, id(), semantic());
}

//! Wrapper for interval variables
/*!
Stores
 - the id of a  temporal variable standing for the start
 - the id of a  temporal variable standing for the end
 - the id of a  Boolean variable standing whether the interval actually is in
the schedule
 */
template <typename T = int> class Interval {
public:
  Interval() {}
  Interval(Solver<T> &solver, const T mindur = 0,
           const T maxdur = Constant::Infinity<T>,
           const T earliest_start = -Constant::Infinity<T>,
           const T latest_start = Constant::Infinity<T>,
           const T earliest_end = -Constant::Infinity<T>,
           const T latest_end = Constant::Infinity<T>,
           const BooleanVar<T> opt = Constant::NoVar);
  //
  //
  //           Solver<T> &s, const T mindur = 0,
  //           const T maxdur = Constant::Infinity<T>,
  //           const BooleanVar<T> opt = Constant::NoVar);

  //    Interval(const Interval<T>&) = default;

  T getEarliestStart(Solver<T> &s) const;
  T getLatestStart(Solver<T> &s) const;
  T getEarliestEnd(Solver<T> &s) const;
  T getLatestEnd(Solver<T> &s) const;

  bool mustExist(Solver<T> &s) const;
  bool cannotExist(Solver<T> &s) const;

  T minDuration(Solver<T> &s) const;
  T maxDuration(Solver<T> &s) const;

  var_t getStart() const;
  var_t getEnd() const;

  int id() const;
  bool operator==(const Interval<T> &t) const;

  std::ostream &display(std::ostream &os) const;

  NumericVar<T> start;
  NumericVar<T> end;
  NumericVar<T> duration;

  bool isOptional() const { return exist.id() != Constant::NoVar; }

  BooleanVar<T> exist{Constant::NoVar};
};

template <typename T = int>
class BooleanExpressionImpl : public ExpressionImpl<T> {
public:
    virtual ~BooleanExpressionImpl() { /*std::cout << "del boolexpr\n";*/ }
    
  virtual var_t id() const override { return self.id(); }
  virtual index_t semantic() const override { return self.semantic(); }

protected:
  BooleanVar<T> self;
};


template <typename T = int>
class NumericExpressionImpl : public ExpressionImpl<T> {
public:
    virtual ~NumericExpressionImpl() { /*std::cout << "del numexpr\n";*/ }
    
  virtual var_t id() const override { return self.id(); }
  virtual T offset() const override { return self.offset(); }

protected:
  NumericVar<T> self;
};

template <typename T = int>
class SumExpressionImpl : public NumericExpressionImpl<T> {
public:
  SumExpressionImpl() {}
    virtual ~SumExpressionImpl() { /*std::cout << "del sum\n";*/ }

    virtual string name() const override { return "sum"; }

    var_t extract(Solver<T> &solver) override {

      std::vector<var_t> vars;
      T lb{0};
      T ub{0};
      for (unsigned i{0}; i < arguments.size(); ++i) {
        auto x{arguments[i]};
        auto w{weights[i]};

        if (w < 0) {
          if (lb != -Constant::Infinity<T>) {
            if (x.max(solver) == Constant::Infinity<T>) {
              lb = -Constant::Infinity<T>;
            } else {
              lb += w * x.max(solver);
            }
          }
          if (ub != Constant::Infinity<T>) {
            if (x.min(solver) == -Constant::Infinity<T>) {
              ub = Constant::Infinity<T>;
            } else {
              ub += w * x.min(solver);
            }
          }
        } else {
          if (lb != -Constant::Infinity<T>) {
            if (x.min(solver) == -Constant::Infinity<T>) {
              lb = -Constant::Infinity<T>;
            } else {
              lb += w * x.min(solver);
            }
          }
          if (ub != Constant::Infinity<T>) {
            if (x.max(solver) == Constant::Infinity<T>) {
              ub = Constant::Infinity<T>;
            } else {
              ub += w * x.max(solver);
            }
          }
          //              lb += w * x.min(solver);
          //              ub += w * x.max(solver);
        }

        x.extract(solver);

        NumericExpressionImpl<T>::self.setOffset(
            NumericExpressionImpl<T>::self.offset() + w * x.offset());
        if (x.min(solver) != x.max(solver)) {
          vars.push_back(x.id());
        }
      }

      if (vars.size() == 1 and weights[0] == 1) {
        NumericExpressionImpl<T>::self.setId(*vars.begin());
      } else {
        NumericExpressionImpl<T>::self = solver.newNumeric(lb, ub);

        vars.push_back(NumericExpressionImpl<T>::id());
        weights.push_back(-1);
        T total{NumericExpressionImpl<T>::self.offset()};

        solver.post(new SumConstraint(solver, vars.begin(), vars.end(),
                                      weights.begin(), total));
        for (auto &w : weights)
          w = -w;
        solver.post(new SumConstraint(solver, vars.begin(), vars.end(),
                                      weights.begin(), -total));
      }

      return NumericExpressionImpl<T>::self.id();
  }

  SumExpressionImpl<T> &addTerm(const NumericVar<T> &x, const T w = 1) {
    arguments.push_back(x);
    weights.push_back(w);
    return *this;
  }
  SumExpressionImpl<T> &addTerm(const T k) {
    NumericExpressionImpl<T>::self.setOffset(
        NumericExpressionImpl<T>::self.offset() + k);
    return *this;
  }

private:
  std::vector<NumericVar<T>> arguments;
  std::vector<T> weights;
};

template <typename T>
NumericVar<T> operator+(const NumericVar<T> &x, const T k) {
  auto sum{new SumExpressionImpl<T>()};
  sum->addTerm(x);
  sum->addTerm(k);
  NumericVar<T> exp(sum);
  return exp;
}

template <typename T>
NumericVar<T> operator+(const NumericVar<T> &x, const NumericVar<T> &y) {
  auto sum{new SumExpressionImpl<T>()};
  sum->addTerm(x);
  sum->addTerm(y);
  NumericVar<T> exp(sum);
  return exp;
}

template <typename varit, typename weightit, typename T>
NumericVar<T> Sum(varit beg_var, varit end_var, weightit beg_weight) {
  auto sum{new SumExpressionImpl<T>()};
  auto wp{beg_weight};
  for (auto xp{beg_var}; xp != end_var; ++xp) {
    sum->addTerm(*xp, *wp);
    ++wp;
  }
  NumericVar<T> exp(sum);
  return exp;
}

// x == y+k
template <typename T = int>
class NumEqExpressionImpl : public ExpressionImpl<T> {
public:
  NumEqExpressionImpl(NumericVar<T> x, NumericVar<T> y, const T k)
      : x(x), y(y), k(k) {}

  virtual string name() const override { return "eq"; }

  var_t extract(Solver<T> &solver) override {
    x.extract(solver);
    y.extract(solver);

    auto prec{x.before(y, -k)};
    auto inf = solver.newDisjunct(~prec, prec);
    auto succ{x.after(y, -k)};
    auto sup = solver.newDisjunct(~succ, succ);
    auto conj = (sup and inf);

    conj.extract(solver);
    self = conj;

    return self.id();
  }

  void post(Solver<T> &solver) override {
    x.extract(solver);
    y.extract(solver);
    auto prec{x.before(y, -k)};
    solver.post(prec);
    auto succ{x.after(y, -k)};
    solver.post(succ);
  }

private:
  BooleanVar<T> self;
  NumericVar<T> x;
  NumericVar<T> y;
  T k;
  // self <-> x - y <= k
  // x <= k (y.id() == Constant::NoVar)
  // y <= k (y.id() == Constant::NoVar, this->k == -k)
};

///// <=
template <typename T>
BooleanVar<T> operator==(const NumericVar<T> &x, const NumericVar<T> &y) {
  BooleanVar<T> exp(new NumEqExpressionImpl<T>(x, y, 0));
  return exp;
}

template <typename T>
BooleanVar<T> operator==(const T x, const NumericVar<T> &y) {
  return operator==(NumericVar<T>(Constant::K, x), y);
}

template <typename T>
BooleanVar<T> operator==(const NumericVar<T> &x, const T y) {
  return operator==(x, NumericVar<T>(Constant::K, y));
}

template <typename T = int> class LeqExpressionImpl : public ExpressionImpl<T> {
public:
  LeqExpressionImpl(NumericVar<T> x, NumericVar<T> y, const T k)
      : x(x), y(y), k(k) {
//
//    std::cout << "leq expr " << x << " + " << x.offset() << " + " << k
//              << " <= " << y << " + " << y.offset() << std::endl;
//    std::cout << "leq expr " << x << " + " << x.offset() << " + " << k
//              << " <= " << y << " + " << y.offset() << std::endl;
  }

  virtual string name() const override { return "leq"; }

  var_t extract(Solver<T> &solver) override {
    x.extract(solver);
    y.extract(solver);

    auto prec{x.before(y, -k)};
    self = solver.newDisjunct(~prec, prec);
    return self.id();
  }

  void post(Solver<T> &solver) override {
    x.extract(solver);
    y.extract(solver);
//
//    std::cout << "post prec\n";
//    std::cout << x << " in [" << x.min(solver) << ".." << x.max(solver)
//              << "] + " << x.offset() << std::endl;
//    std::cout << y << " in [" << y.min(solver) << ".." << y.max(solver)
//              << "] + " << y.offset() << std::endl;

    auto prec{x.before(y, -k)};

//    std::cout << " => " << prec << std::endl;

    solver.post(prec);

//    std::cout << "OK\n";
  }

private:
  BooleanVar<T> self;
  NumericVar<T> x;
  NumericVar<T> y;
  T k;
  // self <-> x - y <= k
  // x <= k (y.id() == Constant::NoVar)
  // y <= k (y.id() == Constant::NoVar, this->k == -k)
};

///// <=
template <typename T>
BooleanVar<T> operator<=(const NumericVar<T> &x, const NumericVar<T> &y) {
  BooleanVar<T> exp(new LeqExpressionImpl<T>(x, y, 0));
  return exp;
}

template <typename T>
BooleanVar<T> operator<=(const T x, const NumericVar<T> &y) {
  return operator<=(NumericVar<T>(Constant::K, x), y);
}

template <typename T>
BooleanVar<T> operator<=(const NumericVar<T> &x, const T y) {
  return operator<=(x, NumericVar<T>(Constant::K, y));
}

///// <
template <typename T>
BooleanVar<T> operator<(const NumericVar<T> &x, const NumericVar<T> &y) {
  BooleanVar<T> exp(new LeqExpressionImpl<T>(x, y, -Gap<T>::epsilon()));
  return exp;
}

template <typename T>
BooleanVar<T> operator<(const T x, const NumericVar<T> &y) {
  BooleanVar<T> exp(new LeqExpressionImpl<T>(NumericVar<T>(Constant::K, x), y,
                                             -Gap<T>::epsilon()));
  return exp;
}

template <typename T>
BooleanVar<T> operator<(const NumericVar<T> &x, const T y) {
  BooleanVar<T> exp(new LeqExpressionImpl<T>(x, NumericVar<T>(Constant::K, y),
                                             -Gap<T>::epsilon()));
  return exp;
}

///// >=
template <typename T>
BooleanVar<T> operator>=(const NumericVar<T> &x, const NumericVar<T> &y) {
  BooleanVar<T> exp(new LeqExpressionImpl<T>(y, x, 0));
  return exp;
}

template <typename T>
BooleanVar<T> operator>=(const T x, const NumericVar<T> &y) {
  BooleanVar<T> exp(
      new LeqExpressionImpl<T>(y, NumericVar<T>(Constant::K, x), 0));
  return exp;
}

template <typename T>
BooleanVar<T> operator>=(const NumericVar<T> &x, const T y) {
  BooleanVar<T> exp(
      new LeqExpressionImpl<T>(NumericVar<T>(Constant::K, y), x, 0));
  return exp;
}

///// >
template <typename T>
BooleanVar<T> operator>(const NumericVar<T> &x, const NumericVar<T> &y) {
  BooleanVar<T> exp(new LeqExpressionImpl<T>(y, x, -Gap<T>::epsilon()));
  return exp;
}

template <typename T>
BooleanVar<T> operator>(const NumericVar<T> &x, const T k) {
  BooleanVar<T> exp(new LeqExpressionImpl<T>(
      NumericVar<T>(Constant::K, k + Gap<T>::epsilon()), x, 0));
  return exp;
}

template <typename T>
BooleanVar<T> operator>(const T k, const NumericVar<T> &x) {
  BooleanVar<T> exp(new LeqExpressionImpl<T>(
      x, NumericVar<T>(Constant::K, k + Gap<T>::epsilon()), 0));
  return exp;
}

template <typename T = int>
class LogicalAndExpressionImpl : public BooleanExpressionImpl<T> {
public:
  template <typename Iter>
  LogicalAndExpressionImpl(Iter beg_var, Iter end_var) {
    for (auto x{beg_var}; x != end_var; ++x) {
      boolean_arguments.emplace_back(*x);
    }
  }

  virtual string name() const override { return "and"; }

  var_t extract(Solver<T> &solver) override {
    std::vector<Literal<T>> L;
    for (auto x : boolean_arguments) {
      x.extract(solver);
      L.push_back(x == false);
    }
      
      auto y = solver.newBoolean();
    L.push_back(y == true);
      
    solver.clauses.add(L.begin(), L.end());
    for (auto x : boolean_arguments) {
      L.clear();
      L.push_back(x == true);
      L.push_back(y == false);
      solver.clauses.add(L.begin(), L.end());
    }
      
      BooleanExpressionImpl<T>::self = y;
    return y.id();
  }

  void post(Solver<T> &solver) override {
    for (auto x : boolean_arguments) {
      x.extract(solver);
      solver.post(x == true);
    }
  }

private:
  std::vector<BooleanVar<T>> boolean_arguments;
};

template <typename T>
BooleanVar<T> operator&&(BooleanVar<T> &x, BooleanVar<T> &y) {
  std::vector<BooleanVar<T>> sc{x, y};
  BooleanVar<T> exp(new LogicalAndExpressionImpl<T>(sc.begin(), sc.end()));
  return exp;
}

template <typename T, typename Iterable> BooleanVar<T> BigAnd(Iterable &X) {
  BooleanVar<T> exp(new LogicalAndExpressionImpl(X.begin(), X.end()));
  return exp;
}

template <typename T = int>
class LogicalOrExpressionImpl : public BooleanExpressionImpl<T> {
public:
  template <typename Iter> LogicalOrExpressionImpl(Iter beg_var, Iter end_var) {
    for (auto x{beg_var}; x != end_var; ++x) {
      boolean_arguments.emplace_back(*x);
    }
  }

  virtual string name() const override { return "or"; }

  var_t extract(Solver<T> &solver) override {
    std::vector<Literal<T>> L;
    for (auto x : boolean_arguments) {
      x.extract(solver);
        L.push_back(x == true);
    }
      auto y = solver.newBoolean();
    L.push_back(y == false);
    solver.clauses.add(L.begin(), L.end());
    for (auto x : boolean_arguments) {
      L.clear();
      L.push_back(x == false);
      L.push_back(y == true);
      solver.clauses.add(L.begin(), L.end());
    }
      BooleanExpressionImpl<T>::self = y;
    return y.id();
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
  std::vector<BooleanVar<T>> boolean_arguments;
};

template <typename T>
BooleanVar<T> operator||(BooleanVar<T> &x, BooleanVar<T> &y) {
  std::vector<BooleanVar<T>> sc{x, y};
  BooleanVar<T> exp(new LogicalOrExpressionImpl<T>(sc.begin(), sc.end()));
  return exp;
}

template <typename T, typename Iterable> BooleanVar<T> BigOr(Iterable &X) {
  BooleanVar<T> exp(new LogicalOrExpressionImpl(X.begin(), X.end()));
  return exp;
}

template <typename T = int>
class LogicalImplicationExpression : public BooleanExpressionImpl<T> {
public:
  LogicalImplicationExpression(BooleanVar<T> x, BooleanVar<T> y)
      : implicant(x), implied(y) {}

  virtual string name() const override { return "implication"; }

  var_t extract(Solver<T> &solver) override {
    implicant.extract(solver);
    implied.extract(solver);
      BooleanExpressionImpl<T>::self = solver.newBoolean();

    std::vector<Literal<T>> cl{solver.boolean.getLiteral(false, implicant),
                               solver.boolean.getLiteral(true, implied),
                               solver.boolean.getLiteral(false, BooleanExpressionImpl<T>::self)};
    solver.clauses.add(cl.begin(), cl.end());

      auto x{BooleanExpressionImpl<T>::self};
    cl = {x == true, implicant == true};
    solver.clauses.add(cl.begin(), cl.end());

    cl = {x == true, implied == false};
    solver.clauses.add(cl.begin(), cl.end());

    return BooleanExpressionImpl<T>::self.id();
  }

  void post(Solver<T> &solver) override {
    implicant.extract(solver);
    implied.extract(solver);
    std::vector<Literal<T>> cl{solver.boolean.getLiteral(false, implicant),
                               solver.boolean.getLiteral(true, implied)};
    solver.clauses.add(cl.begin(), cl.end());
  }

private:
  BooleanVar<T> implicant;
  BooleanVar<T> implied;
};

template <typename T>
BooleanVar<T> BooleanVar<T>::implies(const BooleanVar<T> x) const {
  BooleanVar<T> exp(new LogicalImplicationExpression<T>(*this, x));
  return exp;
}

template <typename T = int>
class CardinalityExpressionImpl : public NumericExpressionImpl<T> {
public:
  template <typename Iter>
  CardinalityExpressionImpl(Iter beg_var, Iter end_var, const T l = 0,
                            const T u = Constant::Infinity<T>)
      : lb(l),
        ub(std::min(u, static_cast<T>(std::distance(beg_var, end_var)))) {
    for (auto x{beg_var}; x != end_var; ++x) {
      boolean_arguments.emplace_back(*x);
    }
  }

  virtual string name() const override { return "cardinality"; }

  var_t extract(Solver<T> &solver) override {
    std::vector<Literal<T>> L;
    for (auto x : boolean_arguments) {
      x.extract(solver);
      L.push_back(x == true);
    }
    NumericExpressionImpl<T>::self = solver.newNumeric();
    solver.post(NumericExpressionImpl<T>::self.after(lb));
    solver.post(NumericExpressionImpl<T>::self.before(ub));
    solver.post(new CardinalityLeqVar<T>(solver, L.begin(), L.end(),
                                         NumericExpressionImpl<T>::self.id()));
    solver.post(new CardinalityGeqVar<T>(solver, L.begin(), L.end(),
                                         NumericExpressionImpl<T>::self.id()));

    return NumericExpressionImpl<T>::self.id();
  }

  void post(Solver<T> &solver) override {
    std::vector<Literal<T>> L;
    for (auto x : boolean_arguments) {
      x.extract(solver);
      L.push_back(x == true);
    }
    T n{static_cast<T>(L.size())};
    if (ub < n) {
      solver.post(new CardinalityConst<T>(solver, L.begin(), L.end(), ub));
    }
    if (lb > 0) {
      for (auto &l : L)
        l = ~l;
      solver.post(new CardinalityConst<T>(solver, L.begin(), L.end(), n - lb));
    }
  }

private:
  T lb{0};
  T ub{Constant::Infinity<T>};
  std::vector<BooleanVar<T>> boolean_arguments;
};

template <typename T, typename Iterable>
NumericVar<T> Cardinality(Iterable &X) {
  NumericVar<T> exp(new CardinalityExpressionImpl(X.begin(), X.end()));
  return exp;
}

template <typename T = int>
class NoOverlapExpressionImpl : public ExpressionImpl<T>,
                                public std::vector<Interval<T>> {
public:
  NoOverlapExpressionImpl(Interval<T> &sched) : schedule(sched) {}

  virtual ~NoOverlapExpressionImpl() {
//    std::cout << "del " << name() << std::endl;
  }

  virtual string name() const override { return "no-overlap"; }

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

          disjunct.push_back(solver.newDisjunct(
              a->end.before(b->start), b->end.before(a->start), x.id()));

          relevant.push_back(x);
        } else if (a->isOptional()) {
          disjunct.push_back(solver.newDisjunct(
              a->end.before(b->start), b->end.before(a->start), a->exist.id()));
          //                relevant.push_back(a->exist);
        } else if (b->isOptional()) {
          disjunct.push_back(solver.newDisjunct(
              a->end.before(b->start), b->end.before(a->start), b->exist.id()));
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

  std::vector<BooleanVar<T>>::iterator begDisjunct() {
    return disjunct.begin();
  }
  std::vector<BooleanVar<T>>::iterator endDisjunct() { return disjunct.end(); }

private:
  Interval<T> schedule;
  std::vector<BooleanVar<T>> disjunct;
  std::vector<BooleanVar<T>> relevant;
};

template <typename T = int> class NoOverlapExpression : public BooleanVar<T> {

public:
  NoOverlapExpression(NoOverlapExpressionImpl<T> *i) : BooleanVar<T>(i) {}

  std::vector<Interval<T>>::iterator begin() {
    return static_cast<NoOverlapExpressionImpl<T> *>(BooleanVar<T>::implem)
        ->begin();
  }
  std::vector<Interval<T>>::iterator end() {
    return static_cast<NoOverlapExpressionImpl<T> *>(BooleanVar<T>::implem)
        ->end();
  }

  std::vector<BooleanVar<T>>::iterator begDisjunct() {
    return static_cast<NoOverlapExpressionImpl<T> *>(BooleanVar<T>::implem)
        ->begDisjunct();
  }
  std::vector<BooleanVar<T>>::iterator endDisjunct() {
    return static_cast<NoOverlapExpressionImpl<T> *>(BooleanVar<T>::implem)
        ->endDisjunct();
  }

  void push_back(const Interval<T> &i) {
    static_cast<NoOverlapExpressionImpl<T> *>(BooleanVar<T>::implem)
        ->push_back(i);
  }
};

template <typename T, typename Iterable>
NoOverlapExpression<T> NoOverlap(Interval<T> &schedule, Iterable &X) {
  auto impl{new NoOverlapExpressionImpl<T>(schedule)};
  for (auto x : X)
    impl->push_back(x);
  NoOverlapExpression<T> exp(impl);
  return exp;
}

template <typename T> NoOverlapExpression<T> NoOverlap(Interval<T> &schedule) {
  auto impl{new NoOverlapExpressionImpl<T>(schedule)};
  NoOverlapExpression<T> exp(impl);
  return exp;
}

/*!
 NumericVar  impl
*/
template<typename T>
T NumericVar<T>::min(Solver<T>& s) const {
  auto v{s.numeric.lower(id())};
  if (v == -Constant::Infinity<T>)
    return v;
  return v + offset();
}

template<typename T>
T NumericVar<T>::max(Solver<T>& s) const {
  auto v{s.numeric.upper(id())};
  if (v == Constant::Infinity<T>)
    return v;
  return v + offset();
}

template <typename T> T NumericVar<T>::earliest(Solver<T> &s) const {
  return min(s);
}

template <typename T> T NumericVar<T>::latest(Solver<T> &s) const {
  return max(s);
}

template <typename T> Literal<T> NumericVar<T>::after(const T t) const {
  return geq<T>(id(), (t == Constant::Infinity<T> ? t : t - offset()));
}

template <typename T> Literal<T> NumericVar<T>::before(const T t) const {
  return leq<T>(id(), (t == Constant::Infinity<T> ? t : t - offset()));
}

template <typename T>
DistanceConstraint<T> NumericVar<T>::after(const NumericVar<T> &e,
                                           const T t) const {
  return e.before(*this, t);
}

template <typename T>
DistanceConstraint<T> NumericVar<T>::before(const NumericVar<T> &e,
                                            const T t) const {
  return {e.id(), id(),
          (t == Constant::Infinity<T> ? t : e.offset() - offset() - t)};
}

/*!
 BooleanVar  impl
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
std::ostream &operator<<(std::ostream &os, const BooleanVar<T> &x) {
  return x.display(os);
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const NumericVar<T> &x) {
  return x.display(os);
}

/*!
 Interval  impl
*/
// template <typename T>
// Interval<T>::Interval(Solver<T> &solver, const T mindur, const T maxdur,
//                       const BooleanVar<T> opt) {
template <typename T>
Interval<T>::Interval(Solver<T> &solver, const T mindur, const T maxdur,
                      const T earliest_start, const T latest_start,
                      const T earliest_end, const T latest_end,
                      const BooleanVar<T> opt) {
  //    : start(solver.newNumeric()),
  //      end(mindur == maxdur ? NumericVar(start.id(), mindur)
  //                           : solver.newNumeric()),
  //      duration(mindur == maxdur ? NumericVar(Constant::K, mindur)
  //                                : solver.newNumeric(mindur, maxdur)),
  //      exist(opt) {

    if(earliest_start == latest_start) {
        start = NumericVar(Constant::K, 0);
    } else {
        start = solver.newNumeric(earliest_start, latest_start);
    }
  if (mindur == maxdur) {
    end = NumericVar(start.id(), mindur);
    duration = NumericVar(Constant::K, mindur);
  } else {
      auto s{start.min(solver)};
      if(s != start.max(solver)) {
          end = solver.newNumeric(earliest_end, latest_end);
          duration = solver.newNumeric(mindur, maxdur);
          solver.post((start + duration) == end);
          solver.post(start.before(end, mindur));
          
          if (maxdur != Constant::Infinity<T>)
              solver.post(end.before(start, -maxdur));
      } else {
          auto ect{std::max(earliest_end, s+mindur)};
          auto lct{std::min(latest_end, s+maxdur)};
          end = solver.newNumeric(ect, lct);
          duration = NumericVar(end.id(), -s);
      }

  }
  exist = opt;
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

template <typename T> T Interval<T>::minDuration(Solver<T> &solver) const {
  return duration.min(solver);
}

template <typename T> T Interval<T>::maxDuration(Solver<T> &solver) const {
  return duration.max(solver);
}

template <typename T> var_t Interval<T>::getStart() const { return start.id(); }

template <typename T> var_t Interval<T>::getEnd() const { return end.id(); }

template <typename T> ostream &Interval<T>::display(ostream &os) const {
  os << "t" << id(); //<< ": [" << start.earliest(solver) << ".." <<
                     // end.latest(solver) << "]";

  os << ": " << start << "/" << duration << "/" << end; //<< "/" << exist;
  return os;
}

template <typename T> ostream &operator<<(ostream &os, const Interval<T> &x) {
  return x.display(os);
}
}

#endif // __MODEL_HPP
