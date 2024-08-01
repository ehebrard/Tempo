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

//#define DBG_EXTRACT
//#define DBG_EXTRACT_SUM

#include <utility>

#include "util/traits.hpp"
#include "Literal.hpp"
#include "constraints/Cardinality.hpp"
#include "constraints/PseudoBoolean.hpp"
#include "constraints/SumConstraint.hpp"


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
  ModelingException(std::string msg) : msg(std::move(msg)) {}

  virtual const char *what() const throw() { return msg.c_str(); }

private:
  std::string msg;
};

//! Interface for expression implementations
template <typename T = int> class ExpressionImpl {
public:
    virtual ~ExpressionImpl() {
#ifdef DBG_EXTRACT
        std::cout << "delete expr\n";
#endif
    }

  // if target is not Constant::NoVar, reference to it instead of creating a
  // new var
  virtual var_t extract(Solver<T> &, const var_t target = Constant::NoVar) = 0;
  virtual void post(Solver<T> &) {
    throw ModelingException("This predicate cannot be a constraint");
  }

  virtual std::string name() const { return "some expression"; }
  virtual var_t id() const = 0; //{ return Constant::NoVar; }
  virtual T offset() const { return 0; }
  virtual index_t semantic() const { return Constant::NoSemantic; }

  //    bool extraction_flag{false};
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
  //  NumInfo(const var_t i, const T o) : _id_(2 * i + 1), _offset(o) {}
  NumInfo(const var_t i, const T o) : _id_(i), _offset(o) {}

  var_t _id_{Constant::NoIndex};
  T _offset{0};
};

template <typename T = int> class NumericVar : public ExpressionFlag {

public:
    constexpr NumericVar() noexcept: ExpressionFlag(false), data(Constant::NoVar, 0) {};
    NumericVar(ExpressionImpl<T> *i) : ExpressionFlag(true), implem(i) {}
    NumericVar(const var_t i, const T o = 0)
        : ExpressionFlag(false), data(i, o) {}

    template<concepts::distance_provider S>
    T min(const S &sc) const;

    template<concepts::distance_provider S>
    T max(const S &sc) const;

    template<concepts::distance_provider S>
    T earliest(const S &) const;

    template<concepts::distance_provider S>
    T latest(const S &) const;

    Literal<T> after(const T t) const;
    Literal<T> before(const T t) const;

    DistanceConstraint<T> after(const NumericVar<T> &e, const T t = 0) const;
    DistanceConstraint<T> before(const NumericVar<T> &e, const T t = 0) const;

    var_t id() const {
      return (ExpressionFlag::_is_expression ? implem->id() : data._id_);
    }
    //    var_t id() const {
    //      return (ExpressionFlag::_is_expression ? implem->id() : data._id_ >>
    //      1);
    //    }
    //    bool sign() const {
    //      return (ExpressionFlag::_is_expression ? 1 : data._id_ & 1);
    //    }

    std::ostream &display(std::ostream &os) const;

    static bool isNumeric() { return true; }

    //    void negate() { data._id_ ^= 1; }
    //    void setSign(const bool s) { data._id_ = data._id_ s; }
    //    void setId(const var_t i, const bool s=1) { data._id_ = 2 * i + s; }
    void setId(const var_t i) { data._id_ = i; }

    T offset() const {
      return (ExpressionFlag::_is_expression ? implem->offset() : data._offset);
      //        return data._offset;
    }

    void setOffset(const T o) { data._offset = o; }

    void extract(Solver<T> &solver, const var_t target = Constant::NoVar) {

#ifdef DBG_EXTRACT
      std::cout << "beg extract " << *this << std::endl;
#endif

      if (ExpressionFlag::_is_expression) {

#ifdef DBG_EXTRACT
        std::cout << " (expr)" << std::endl;
#endif

        if (implem->id() == Constant::NoVar) {

#ifdef DBG_EXTRACT
          std::cout << " extract expr" << std::endl;
#endif

          implem->extract(solver, target);

#ifdef DBG_EXTRACT
          std::cout << " ==> " << implem->id() << std::endl;
#endif

          solver.trash_bin.push_back(implem);
        }

#ifdef DBG_EXTRACT
        else
          std::cout << " expr already extracted" << std::endl;
#endif
        data = NumInfo<T>(implem->id(), implem->offset());
        ExpressionFlag::_is_expression = false;
      }

#ifdef DBG_EXTRACT
      std::cout << "end extract " << *this << std::endl;
#endif
    }

    void post(Solver<T> &solver) {
      if (ExpressionFlag::_is_expression) {
        throw ModelingException("Numeric expression cannot be constraints");
      } else {
        solver.addToSearch(*this);
      }
    }

    //  protected:
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
  constexpr BooleanVar() noexcept: ExpressionFlag(false), data(Constant::NoIndex, Constant::NoSemantic) {}
  BooleanVar(ExpressionImpl<T> *i) : ExpressionFlag(true), implem(i) {}
  BooleanVar(const var_t i, const index_t f = 0)
      : ExpressionFlag(false), data(i, f) {}

  Literal<T> operator==(const bool t) const;

  BooleanVar<T> implies(const BooleanVar<T> x) const;

var_t id() const { return (ExpressionFlag::_is_expression ? implem->id() : data._id_); }

    index_t semantic() const { return (ExpressionFlag::_is_expression ? implem->semantic() : data._edge_id_); }
    
    bool isTrue(Solver<T>& s) const { return s.boolean.isTrue(id()); }
    bool isFalse(Solver<T>& s) const { return s.boolean.isFalse(id()); }

  operator var_t() const { return id(); }

  std::ostream &display(std::ostream &os) const;

  static bool isNumeric() { return false; }

  void setId(const var_t i) { data._id_ = i; }

  void extract(Solver<T> &solver, const var_t target = Constant::NoVar) {

#ifdef DBG_EXTRACT
    std::cout << "beg extract " << *this << std::endl;
#endif

    if (ExpressionFlag::_is_expression) {

#ifdef DBG_EXTRACT
      std::cout << " (expr)" << std::endl;
#endif

      if (implem->id() == Constant::NoVar) {

#ifdef DBG_EXTRACT
        std::cout << " extract expr" << std::endl;
#endif

        implem->extract(solver, target);

#ifdef DBG_EXTRACT
          std::cout << " ==> " << implem->id() << std::endl;
#endif

        solver.trash_bin.push_back(implem);
      }

#ifdef DBG_EXTRACT
      else
        std::cout << " expr already extracted" << std::endl;
#endif

      data = BoolInfo<T>(implem->id(), implem->semantic());
      ExpressionFlag::_is_expression = false;
    }

#ifdef DBG_EXTRACT
    std::cout << "end extract " << *this << std::endl;
#endif
  }

  void post(Solver<T> &solver) {

#ifdef DBG_EXTRACT
    std::cout << "beg post" << std::endl;
#endif

    if (ExpressionFlag::_is_expression) {

        if (implem->id() == Constant::NoVar) {
            implem->post(solver);
        }

      ExpressionFlag::_is_expression = false;
    } else {
      solver.addToSearch(*this);
    }

#ifdef DBG_EXTRACT
    std::cout << "end post" << std::endl;
#endif
  }

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
protected:
  Interval(NumericVar<T> start, NumericVar<T> end, NumericVar<T> duration) :
    _id_(start.id()), start(start), end(end), duration(duration) {}
public:
  constexpr Interval() noexcept : start(), end(), duration() {}


  Interval(Solver<T> &solver, const T mindur = 0,
           const T maxdur = Constant::Infinity<T>,
           const T earliest_start = -Constant::Infinity<T>,
           const T latest_start = Constant::Infinity<T>,
           const T earliest_end = -Constant::Infinity<T>,
           const T latest_end = Constant::Infinity<T>,
           const BooleanVar<T> opt = Constant::True);

  Interval(Solver<T> &solver, const NumericVar<T> s, const NumericVar<T> e,
           const NumericVar<T> d, const BooleanVar<T> opt = Constant::True);

  template<concepts::distance_provider S>
  T getEarliestStart(const S &s) const;

  template<concepts::distance_provider S>
  T getLatestStart(const S &s) const;

  template<concepts::distance_provider S>
  T getEarliestEnd(const S &s) const;

  template<concepts::distance_provider S>
  T getLatestEnd(const S &s) const;

  bool mustExist(Solver<T> &s) const;
  bool cannotExist(Solver<T> &s) const;

  template<concepts::distance_provider S>
  T minDuration(const S &s) const;
  template<concepts::distance_provider S>
  T maxDuration(const S &s) const;

  var_t getStart() const;
  var_t getEnd() const;

  int id() const;
  bool operator==(const Interval<T> &t) const;

  std::ostream &display(std::ostream &os) const;

  index_t _id_{Constant::NoIndex};

  NumericVar<T> start;
  NumericVar<T> end;
  NumericVar<T> duration;

  //  bool isOptional() const { return exist.id() != Constant::NoVar; }
  bool isOptional(Solver<T> &s) const { return not exist.isTrue(s); }

    BooleanVar<T> exist{Constant::True};

    static index_t num_intervals;
};

template <typename T> index_t Interval<T>::num_intervals = 0;

template <typename T = int>
class BooleanExpressionImpl : public ExpressionImpl<T> {
public:
    virtual ~BooleanExpressionImpl() {
#ifdef DBG_EXTRACT
        std::cout << "delete boolexpr\n";
#endif
    }

  virtual var_t id() const override { return self.id(); }
  virtual index_t semantic() const override { return self.semantic(); }

protected:
  BooleanVar<T> self{Constant::NoVar};
};

template <typename T = int>
class NumericExpressionImpl : public ExpressionImpl<T> {
public:
    virtual ~NumericExpressionImpl() {
#ifdef DBG_EXTRACT
        std::cout << "delete numexpr\n";
#endif
    }

  virtual var_t id() const override { return self.id(); }
  virtual T offset() const override { return self.offset(); }

protected:
  NumericVar<T> self{Constant::NoVar};
};

template <typename T = int>
class SumExpressionImpl : public NumericExpressionImpl<T> {
public:
  SumExpressionImpl() {}
    virtual ~SumExpressionImpl() {
#ifdef DBG_EXTRACT
        std::cout << "delete sum expression\n";
#endif
    }

  virtual std::string name() const override { return "sum"; }

  void increaseBias(const T v) {
    NumericExpressionImpl<T>::self.setOffset(
        NumericExpressionImpl<T>::self.offset() + v);
  }

  void preprocessNumeric(Solver<T> &solver) {

    if (not numeric_arguments.empty()) {
      // preprocessing to remove duplicates and constants
      // there has been a first preprocessing to move all the variables' offsets
      // into self._offset

      // order the arguments
      std::vector<index_t> order;
      for (index_t i{0}; i < static_cast<index_t>(numeric_arguments.size());
           ++i) {
        order.push_back(i);
      }
      std::sort(order.begin(), order.end(),
                [&](const index_t i, const index_t j) {
                  return numeric_arguments[i].id() < numeric_arguments[j].id();
                });

      // store a copy of the arguments (in index order)
      std::vector<NumericVar<T>> args;
      std::vector<var_t> vars;
      std::vector<T> ws;

      for (index_t i{0}; i < static_cast<index_t>(numeric_arguments.size());
           ++i) {
        args.push_back(numeric_arguments[order[i]]);
        ws.push_back(num_weights[order[i]]);
        args.back().extract(solver);
      }

#ifdef DBG_EXTRACT_SUM
      numeric_arguments = args;
      num_weights = ws;
      this->display(std::cout);
      std::cout << std::endl;
#endif

      // the actual preprocessing
      //            T lb{0};
      //            T ub{0};
      for (index_t i{static_cast<index_t>(args.size())}; i-- > 0;) {
        bool rm{false};
        auto x{args[i]};
        auto w{ws[i]};

#ifdef DBG_EXTRACT_SUM
        std::cout << w << "*[" << x << "]";
#endif

        if (w == 0) {
          // coefficient is 0 -> just ignore
          rm = true;
        } else if (i > 0 and x.id() == args[i - 1].id()) {
          // duplicate -> increase
          ws[i - 1] += w;
          rm = true;
        } else if (x.min(solver) == x.max(solver)) {
          // constant (after extraction) -> move it to the bias
#ifdef DBG_EXTRACT_SUM
          std::cout << ": => increase the bias by " << (w * x.min(solver));
#endif
          increaseBias(w * x.min(solver));
          rm = true;
        }

#ifdef DBG_EXTRACT_SUM
        std::cout << std::endl;
#endif

        if (rm) {

#ifdef DBG_EXTRACT_SUM
          std::cout << " rm term " << w << "*" << x << std::endl;
#endif
          args[i] = args.back();
          ws[i] = ws.back();
          args.pop_back();
          ws.pop_back();
        }
#ifdef DBG_EXTRACT_SUM
        else {
          std::cout << " keep arg " << args[i] << std::endl;
        }
#endif

        //          increaseBias(ws[i] * args[i].offset());
      }

#ifdef DBG_EXTRACT_SUM
      std::cout << args.size() << " arguments remaining\n";
      numeric_arguments = args;
      num_weights = ws;
      this->display(std::cout);
      std::cout << std::endl;
#endif

      // compute the initial bounds
      for (index_t i{0}; i < static_cast<index_t>(args.size()); ++i) {
        auto x{args[i]};
        auto w{ws[i]};

        //                vars.push_back(x.id());
        if (w < 0) {
          if (num_lb != -Constant::Infinity<T>) {
            if (x.max(solver) == Constant::Infinity<T>) {
              num_lb = -Constant::Infinity<T>;
            } else {
              num_lb += w * x.max(solver);
            }
          }
          if (num_ub != Constant::Infinity<T>) {
            if (x.min(solver) == -Constant::Infinity<T>) {
              num_ub = Constant::Infinity<T>;
            } else {
              num_ub += w * x.min(solver);
            }
          }
        } else {
          if (num_lb != -Constant::Infinity<T>) {
            if (x.min(solver) == -Constant::Infinity<T>) {
              num_lb = -Constant::Infinity<T>;
            } else {
              num_lb += w * x.min(solver);
            }
          }
          if (num_ub != Constant::Infinity<T>) {
            if (x.max(solver) == Constant::Infinity<T>) {
              num_ub = Constant::Infinity<T>;
            } else {
              num_ub += w * x.max(solver);
            }
          }
          //              lb += w * x.min(solver);
          //              ub += w * x.max(solver);
        }
      }

      numeric_arguments = args;
      num_weights = ws;
    }
  }

  void preprocessBoolean(Solver<T> &solver) {

    if (not boolean_arguments.empty()) {
      // preprocessing to remove duplicates and constants

      // order the arguments
      std::vector<index_t> order;
      for (index_t i{0}; i < static_cast<index_t>(boolean_arguments.size());
           ++i) {
        order.push_back(i);
      }
      std::sort(order.begin(), order.end(),
                [&](const index_t i, const index_t j) {
                  return boolean_arguments[i].id() < boolean_arguments[j].id();
                });

      // store a copy of the arguments (in index order)
      std::vector<BooleanVar<T>> args;
      std::vector<var_t> vars;
      std::vector<T> ws;

      for (index_t i{0}; i < static_cast<index_t>(boolean_arguments.size());
           ++i) {
        args.push_back(boolean_arguments[order[i]]);
        ws.push_back(bool_weights[order[i]]);
        args.back().extract(solver);
      }

#ifdef DBG_EXTRACT_SUM
      boolean_arguments = args;
      bool_weights = ws;
      this->display(std::cout);
      std::cout << std::endl;
#endif

      // the actual preprocessing
      //            T lb{0};
      //            T ub{0};
      for (index_t i{static_cast<index_t>(args.size())}; i-- > 0;) {
        bool rm{false};
        auto x{args[i]};
        auto w{ws[i]};

#ifdef DBG_EXTRACT_SUM
        std::cout << w << "*[" << x << "]";
#endif

        if (w == 0 or x.isFalse(solver)) {
          // coefficient is 0 -> just ignore
          rm = true;
        } else if (i > 0 and x.id() == args[i - 1].id()) {
          // duplicate -> increase
          ws[i - 1] += w;
          rm = true;
        } else if (x.isTrue(solver)) {
          // constant (after extraction) -> move it to the bias
#ifdef DBG_EXTRACT_SUM
          std::cout << ": => increase the bias by " << w;
#endif
          increaseBias(w);
          rm = true;
        }

#ifdef DBG_EXTRACT_SUM
        std::cout << std::endl;
#endif

        if (rm) {

#ifdef DBG_EXTRACT_SUM
          std::cout << " rm term " << w << "*" << x << std::endl;
#endif
          args[i] = args.back();
          ws[i] = ws.back();
          args.pop_back();
          ws.pop_back();
        }
#ifdef DBG_EXTRACT_SUM
        else {
          std::cout << " keep arg " << args[i] << std::endl;
        }
#endif

        //          increaseBias(ws[i] * args[i].offset());
      }

#ifdef DBG_EXTRACT_SUM
      std::cout << args.size() << " arguments remaining\n";
      boolean_arguments = args;
      bool_weights = ws;
      this->display(std::cout);
      std::cout << std::endl;
#endif

      // compute the initial bounds
      for (index_t i{0}; i < static_cast<index_t>(args.size()); ++i) {
        auto w{ws[i]};
        if (w < 0) {
          bool_lb += w;
        } else {
          bool_ub += w;
        }
      }

      boolean_arguments = args;
      bool_weights = ws;
    }
  }

  var_t extract(Solver<T> &solver,
                const var_t target = Constant::NoVar) override {

#ifdef DBG_EXTRACT_SUM
    this->display(std::cout);
    std::cout << std::endl;
#endif

    //      NumericVar<T> num_part;

    preprocessBoolean(solver);

    if (not boolean_arguments.empty()) {
      //        std::cout << bool_lb << ".." << bool_ub << std::endl;
      std::vector<Literal<T>> L;
      for (auto x : boolean_arguments) {
        L.push_back(x == true);
      }

      auto bool_part{solver.newNumeric(bool_lb, bool_ub)};

      solver.post(new PseudoBooleanLeqVar<T>(
          solver, L.begin(), L.end(), bool_weights.begin(), bool_part.id()));

      solver.post(new PseudoBooleanGeqVar<T>(
          solver, L.begin(), L.end(), bool_weights.begin(), bool_part.id()));

      assert(numeric_arguments.empty());

      NumericExpressionImpl<T>::self = bool_part;

      return NumericExpressionImpl<T>::self.id();
    }

    preprocessNumeric(solver);

    auto bias{NumericExpressionImpl<T>::self.offset()};

    if (not numeric_arguments.empty()) {

      if (not boolean_arguments.empty()) {
        std::cout << "TODO: mixed sum!\n";
        exit(0);
      }

      std::vector<var_t> nvars;
      for (auto x : numeric_arguments)
        nvars.push_back(x.id());

      if (nvars.size() == 1 and num_weights[0] == 1) {
        NumericExpressionImpl<T>::self.setId(*nvars.begin());
        NumericExpressionImpl<T>::self._is_expression = false;
        //                      num_part.setId(*nvars.begin());
        //                      num_part.setOffset(bias);
        //                      num_part._is_expression = false;
      } else {
        if (target == Constant::NoVar) {
          NumericExpressionImpl<T>::self = solver.newNumeric(num_lb, num_ub);
          //                            num_part = solver.newNumeric(num_lb,
          //                            num_ub);
        } else {
          NumericExpressionImpl<T>::self.data = NumInfo<T>(target, -bias);
          //                            num_part.data = NumInfo<T>(target,
          //                            -bias);
          solver.post(NumericExpressionImpl<T>::self >= num_lb);
          solver.post(NumericExpressionImpl<T>::self <= num_ub);
        }

        nvars.push_back(NumericExpressionImpl<T>::id());
        //                      nvars.push_back(num_part.id());
        num_weights.push_back(-1);
        //        T total{offset};

        solver.post(new SumConstraint(solver, nvars.begin(), nvars.end(),
                                      num_weights.begin(), -bias));
        for (auto &w : num_weights)
          w = -w;
        solver.post(new SumConstraint(solver, nvars.begin(), nvars.end(),
                                      num_weights.begin(), bias));
      }

      //
    } else {

      NumericExpressionImpl<T>::self.setId(Constant::K);
      NumericExpressionImpl<T>::self._is_expression = false;
      //            num_part.setId(Constant::K);
      //            num_part.setOffset(bias);
      //            num_part._is_expression = false;
    }

    //      NumericExpressionImpl<T>::self = num_part;

//    std::cout << " => " << NumericExpressionImpl<T>::self << std::endl;

    return NumericExpressionImpl<T>::self.id();
  }

  SumExpressionImpl<T> &addTerm(const NumericVar<T> &x, const T w = 1) {
    increaseBias(w * x.offset());
    if (x.id() != Constant::K) {
      numeric_arguments.push_back(x);
      num_weights.push_back(w);
    }
    return *this;
  }

  SumExpressionImpl<T> &addTerm(const BooleanVar<T> &x, const T w = 1) {
    assert(x.id() != Constant::K);

    boolean_arguments.push_back(x);
    bool_weights.push_back(w);

    return *this;
  }

  SumExpressionImpl<T> &addTerm(const T k) {
    increaseBias(k);
    return *this;
  }

  std::ostream &display(std::ostream &os) const;

private:
  T num_lb{0};
  T num_ub{0};
  T bool_lb{0};
  T bool_ub{0};
  std::vector<NumericVar<T>> numeric_arguments;
  std::vector<T> num_weights;
  std::vector<BooleanVar<T>> boolean_arguments;
  std::vector<T> bool_weights;
  //  T bias{0};
};

template <typename T>
std::ostream &SumExpressionImpl<T>::display(std::ostream &os) const {
  os << "{";
  if (not numeric_arguments.empty()) {
    os << num_weights[0] << "*" << numeric_arguments[0];
    for (unsigned i{1}; i < numeric_arguments.size(); ++i) {
      os << " + " << num_weights[i] << "*" << numeric_arguments[i];
    }
  }
  if (not boolean_arguments.empty()) {
    os << bool_weights[0] << "*" << boolean_arguments[0];
    for (unsigned i{1}; i < boolean_arguments.size(); ++i) {
      os << " + " << bool_weights[i] << "*" << boolean_arguments[i];
    }
  }
  os << "} bias=" << NumericExpressionImpl<T>::self.offset();
  return os;
}

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

template <typename T>
NumericVar<T> operator-(const NumericVar<T> &x, const T k) {
  auto sum{new SumExpressionImpl<T>()};
  sum->addTerm(x);
  sum->addTerm(-k);
  NumericVar<T> exp(sum);
  return exp;
}

template <typename T>
NumericVar<T> operator-(const NumericVar<T> &x, const NumericVar<T> &y) {
  auto sum{new SumExpressionImpl<T>()};
  sum->addTerm(x);
  sum->addTerm(y, -1);
  NumericVar<T> exp(sum);
  return exp;
}

// template <typename varit, typename weightit, typename T>
// NumericVar<T> Sum(varit beg_var, varit end_var, weightit beg_weight) {
//   auto sum{new SumExpressionImpl<T>()};
//   auto wp{beg_weight};
//   for (auto xp{beg_var}; xp != end_var; ++xp) {
//     sum->addTerm(*xp, *wp);
//     ++wp;
//   }
//   NumericVar<T> exp(sum);
//   return exp;
// }

template <typename T = int>
NumericVar<T> Sum(const std::vector<NumericVar<T>> &X,
                  const std::vector<T> &W) {
  auto sum{new SumExpressionImpl<T>()};
  auto wp{W.begin()};
  for (auto xp{X.begin()}; xp != X.end(); ++xp) {
    sum->addTerm(*xp, *wp);
    ++wp;
  }
  NumericVar<T> exp(sum);
  return exp;
}

template <typename T = int>
NumericVar<T> Sum(const std::vector<NumericVar<T>> &X) {
  auto sum{new SumExpressionImpl<T>()};
  for (auto xp{X.begin()}; xp != X.end(); ++xp) {
    sum->addTerm(*xp, 1);
  }
  NumericVar<T> exp(sum);
  return exp;
}

template <typename T = int>
NumericVar<T> Sum(const std::vector<BooleanVar<T>> &X,
                  const std::vector<T> &W) {
  auto sum{new SumExpressionImpl<T>()};
  auto wp{W.begin()};
  for (auto xp{X.begin()}; xp != X.end(); ++xp) {
    sum->addTerm(*xp, *wp);
    ++wp;
  }
  NumericVar<T> exp(sum);
  return exp;
}

template <typename T = int>
NumericVar<T> Sum(const std::vector<NumericVar<T>> &X,
                  const std::vector<BooleanVar<T>> &B,
                  const std::vector<T> &W) {
  auto sum{new SumExpressionImpl<T>()};
  auto wp{W.begin()};
  for (auto xp{X.begin()}; xp != X.end(); ++xp) {
    sum->addTerm(*xp, *wp);
    ++wp;
  }
  for (auto xp{B.begin()}; xp != B.end(); ++xp) {
    sum->addTerm(*xp, *wp);
    ++wp;
  }
  NumericVar<T> exp(sum);
  return exp;
}

// x == y+k
template <typename T = int>
class NumEqExpressionImpl : public BooleanExpressionImpl<T> {
public:
  NumEqExpressionImpl(NumericVar<T> x, NumericVar<T> y, const T k)
      : x(x), y(y), k(k) {}

    virtual ~NumEqExpressionImpl() {
#ifdef DBG_EXTRACT
        std::cout << "delete eq expression\n";
#endif
    }

  virtual std::string name() const override { return "eq"; }

  var_t extract(Solver<T> &solver, const var_t = Constant::NoVar) override {
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
    if (x._is_expression) {
      y.extract(solver);
      x.extract(solver, y.id());
    } else if (y._is_expression) {
      x.extract(solver);
      y.extract(solver, x.id());
    } else {
      x.extract(solver);
      y.extract(solver);
      auto prec{x.before(y, -k)};
      solver.post(prec);
      auto succ{x.after(y, -k)};
      solver.post(succ);
    }
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

template <typename T = int> class LeqExpressionImpl : public BooleanExpressionImpl<T> {
public:
  LeqExpressionImpl(NumericVar<T> x, NumericVar<T> y, const T k)
      : x(x), y(y), k(k) {
//
//    std::cout << "leq expr " << x << " + " << x.offset() << " + " << k
//              << " <= " << y << " + " << y.offset() << std::endl;
//    std::cout << "leq expr " << x << " + " << x.offset() << " + " << k
//              << " <= " << y << " + " << y.offset() << std::endl;
  }

  virtual std::string name() const override { return "leq"; }

  var_t extract(Solver<T> &solver, const var_t = Constant::NoVar) override {
    x.extract(solver);
    y.extract(solver);
    auto prec{x.before(y, -k)};
    self = solver.newDisjunct(~prec, prec);
    return self.id();
  }

  void post(Solver<T> &solver) override {
    x.extract(solver);
    y.extract(solver);

    auto prec{x.before(y, -k)};
    solver.post(prec);
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

  virtual std::string name() const override { return "and"; }

  var_t extract(Solver<T> &solver, const var_t = Constant::NoVar) override {
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

  virtual std::string name() const override { return "or"; }

  var_t extract(Solver<T> &solver, const var_t = Constant::NoVar) override {
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

  virtual std::string name() const override { return "implication"; }

  var_t extract(Solver<T> &solver, const var_t = Constant::NoVar) override {
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

  virtual std::string name() const override { return "cardinality"; }

  var_t extract(Solver<T> &solver, const var_t = Constant::NoVar) override {
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

template <typename T>
NumericVar<T> Cardinality(const std::vector<BooleanVar<T>> &X) {
  NumericVar<T> exp(new CardinalityExpressionImpl(X.begin(), X.end()));
  return exp;
}

template <typename T = int>
class NoOverlapExpressionImpl : public BooleanExpressionImpl<T>,
                                public std::vector<Interval<T>> {
public:
  NoOverlapExpressionImpl(Interval<T> &sched) : schedule(sched) {}

  virtual ~NoOverlapExpressionImpl() {
#ifdef DBG_EXTRACT
    std::cout << "delete " << name() << std::endl;
#endif
  }

  virtual std::string name() const override { return "no-overlap"; }

  var_t extract(Solver<T> &, const var_t) override {
    throw ModelingException("NoOverlap is not a predicate");
    return Constant::NoVar;
  }

  void post(Solver<T> &solver) override {
    size_t k{0};

    bool opt_flag{false};
    for (auto a{this->begin()}; a != this->end(); ++a) {
      for (auto b{a + 1}; b != this->end(); ++b) {
        auto t_ab{0};
        auto t_ba{0};
        auto ai{static_cast<size_t>(a - this->begin())};
        auto bi{static_cast<size_t>(b - this->begin())};
        if (ai < transitions.size() and bi < transitions[ai].size()) {
          t_ab = transitions[ai][bi];
        }
        if (bi < transitions.size() and ai < transitions[bi].size()) {
          t_ba = transitions[bi][ai];
        }

        auto a_before_b{a->end.before(b->start, t_ab)};
        auto b_before_a{b->end.before(a->start, t_ba)};

        if (a->isOptional(solver) or b->isOptional(solver)) {

          auto x_ab{solver.newDisjunct(Constant::NoEdge<T>, a_before_b)};
          auto x_ba{solver.newDisjunct(Constant::NoEdge<T>, b_before_a)};

          // (a->exist and b->exist) -> (x_ab or x_ba)
          // ~a->exist or ~b->exist or x_ab or x_ba
          std::vector<Literal<T>> cl{x_ab == true, x_ba == true};
          if (a->isOptional(solver))
            cl.push_back(a->exist == false);
          if (b->isOptional(solver))
            cl.push_back(b->exist == false);
          solver.clauses.add(cl.begin(), cl.end());

          disjunct.push_back(x_ab);
          disjunct.push_back(x_ba);

          opt_flag = true;
        }

        else {
          disjunct.push_back(solver.newDisjunct(a_before_b, b_before_a));
        }

        while (k < disjunct.size())
          solver.addToSearch(disjunct[k++]);
      }
    }

      if(std::distance(this->begin(),this->end()) > 1) {
        if (solver.getOptions().edge_finding) {
          if (opt_flag) {
            std::cout
                << "edge-finding for optional intervals is not implemented, "
                   "please use '--no-edge-finding'. Aborthing\n";
            exit(0);
          }

          solver.postEdgeFinding(schedule, this->begin(), this->end(),
                                 this->begDisjunct());
        }

        if (solver.getOptions().transitivity) {

          if (opt_flag) {
            std::cout
                << "transitivity for optional intervals is not implemented, "
                   "please use '--no-transitivity'. Aborthing\n";
            exit(0);
          }

          solver.postTransitivity(schedule, this->begin(), this->end(),
                                  this->begDisjunct());
        }
      }
  }

  template <typename Matrix> void setTransisions(const Matrix &D) {
    transitions.resize(D.size());
    index_t i{0};
    for (auto &row : D) {
      for (auto t : row) {
        transitions[i].push_back(t);
      }
      ++i;
    }
  }

  std::vector<BooleanVar<T>>::iterator begDisjunct() {
    return disjunct.begin();
  }
  std::vector<BooleanVar<T>>::iterator endDisjunct() { return disjunct.end(); }

private:
  Interval<T> schedule;
  std::vector<BooleanVar<T>> disjunct;
//  std::vector<BooleanVar<T>> relevant;
  std::vector<std::vector<T>> transitions;
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

  template <typename Matrix> void setTransisions(const Matrix &D) {
    static_cast<NoOverlapExpressionImpl<T> *>(BooleanVar<T>::implem)
        ->setTransisions(D);
  }
};

template <typename T, typename Iterable, typename Matrix>
NoOverlapExpression<T> NoOverlap(Interval<T> &schedule, const Iterable &X,
                                 const Matrix &D) {
  auto impl{new NoOverlapExpressionImpl<T>(schedule)};
  //  for (auto x : X)
  //    impl->push_back(x);
  NoOverlapExpression<T> exp(impl);

  for (auto x : X)
    exp.push_back(x);

  exp.setTransisions(D);

  return exp;
}

template <typename T> NoOverlapExpression<T> NoOverlap(Interval<T> &schedule) {
  auto impl{new NoOverlapExpressionImpl<T>(schedule)};
  NoOverlapExpression<T> exp(impl);
  return exp;
}

/*!
 Cumulative Resource
 */

template <typename T = int>
class CumulativeExpressionImpl : public BooleanExpressionImpl<T>,
                                 public std::vector<Interval<T>> {
public:
  CumulativeExpressionImpl(const NumericVar<T> c) : capacity(c) {}

  virtual ~CumulativeExpressionImpl() {
#ifdef DBG_EXTRACT
    std::cout << "delete " << name() << std::endl;
#endif
  }

  virtual std::string name() const override { return "cumulative"; }

  var_t extract(Solver<T> &, const var_t) override {
    throw ModelingException("Cumulative is not a predicate");
    return Constant::NoVar;
  }

  void provide(const NumericVar<T> d, const Interval<T> i) {
    this->push_back(i);
    demand.push_back(d);
  }

  void post(Solver<T> &solver) override {
    size_t k{0};

    for (auto a{this->begin()}; a != this->end(); ++a) {
      for (auto b{a + 1}; b != this->end(); ++b) {

        auto ea_before_sb{a->end.before(b->start)};
        auto eb_before_sa{b->end.before(a->start)};

        disjunct.push_back(solver.newDisjunct(~ea_before_sb, ea_before_sb));
        disjunct.push_back(solver.newDisjunct(~eb_before_sa, eb_before_sa));

        while (k < disjunct.size())
          solver.addToSearch(disjunct[k++]);
      }
    }

    if (std::distance(this->begin(), this->end()) > 1) {
      solver.postCumulative(capacity, this->begin(), this->end(),
                            this->begDemand(), this->begDisjunct());
    }
  }

  std::vector<BooleanVar<T>>::iterator begDisjunct() {
    return disjunct.begin();
  }
  std::vector<BooleanVar<T>>::iterator endDisjunct() { return disjunct.end(); }

  std::vector<NumericVar<T>>::iterator begDemand() { return demand.begin(); }
  std::vector<NumericVar<T>>::iterator endDemand() { return demand.end(); }

private:
  NumericVar<T> capacity;
  std::vector<BooleanVar<T>> disjunct;
  std::vector<NumericVar<T>> demand;
};

template <typename T = int> class CumulativeExpression : public BooleanVar<T> {

public:
  CumulativeExpression(CumulativeExpressionImpl<T> *i) : BooleanVar<T>(i) {}

  std::vector<Interval<T>>::iterator begin() {
    return static_cast<CumulativeExpressionImpl<T> *>(BooleanVar<T>::implem)
        ->begin();
  }
  std::vector<Interval<T>>::iterator end() {
    return static_cast<CumulativeExpressionImpl<T> *>(BooleanVar<T>::implem)
        ->end();
  }

  std::vector<BooleanVar<T>>::iterator begDisjunct() {
    return static_cast<CumulativeExpressionImpl<T> *>(BooleanVar<T>::implem)
        ->begDisjunct();
  }
  std::vector<BooleanVar<T>>::iterator endDisjunct() {
    return static_cast<CumulativeExpressionImpl<T> *>(BooleanVar<T>::implem)
        ->endDisjunct();
  }

  std::vector<BooleanVar<T>>::iterator begDemand() {
    return static_cast<CumulativeExpressionImpl<T> *>(BooleanVar<T>::implem)
        ->begDemand();
  }
  std::vector<BooleanVar<T>>::iterator endDemand() {
    return static_cast<CumulativeExpressionImpl<T> *>(BooleanVar<T>::implem)
        ->endDemand();
  }

  //    void provide(const NumericVar<T> d, const Interval<T> i) {
  //        static_cast<CumulativeExpressionImpl<T> *>(BooleanVar<T>::implem)
  //            ->provide(d,i);
  //    }
};

template <typename T>
CumulativeExpression<T> Cumulative(const NumericVar<T> c,
                                   const std::vector<Interval<T>> &I,
                                   const std::vector<NumericVar<T>> &D) {
  auto impl{new CumulativeExpressionImpl<T>(c)};

  auto i{I.begin()};
  auto d{D.begin()};
  while (i != I.end()) {
    impl->provide(*d, *i);
    ++i;
    ++d;
  }

  CumulativeExpression<T> exp(impl);

  return exp;
}

/*!
 NumericVar  impl
*/
template<typename T>
template<concepts::distance_provider S>
T NumericVar<T>::min(const S& s) const {
  T v = s.numeric.lower(id());
  if (v == -Constant::Infinity<T>)
    return v;
  return v + offset();
}

template<typename T>
template<concepts::distance_provider S>
T NumericVar<T>::max(const S& s) const {
  T v = s.numeric.upper(id());
  if (v == Constant::Infinity<T>)
    return v;
  return v + offset();
}

template <typename T>
template<concepts::distance_provider S>
T NumericVar<T>::earliest(const S &s) const {
  return min(s);
}

template <typename T>
template<concepts::distance_provider S>
T NumericVar<T>::latest(const S &s) const {
  return max(s);
}

template <typename T> Literal<T> NumericVar<T>::after(const T t) const {
  assert(t != Constant::Infinity<T>);
  assert(t != -Constant::Infinity<T>);

  //    if(sign())
          return geq<T>(id(), (t == Constant::Infinity<T> ? t : t - offset()));
  //    else
  //        return leq<T>(id(), (t == Constant::Infinity<T> ? t : t -
  //        offset()));

  //    if(sign())
//  return geq<T>(id(), t - offset());
  //    else
  //        return leq<T>(id(), offset() - t);
}

template <typename T> Literal<T> NumericVar<T>::before(const T t) const {
//  assert(t != Constant::Infinity<T>);
//  assert(t != -Constant::Infinity<T>);

    return leq<T>(id(), (t == Constant::Infinity<T> ? t : t - offset()));

  //    if(sign())
//  return leq<T>(id(), t - offset());
  //    else
  //        return geq<T>(id(), offset() - t);
}

template <typename T>
DistanceConstraint<T> NumericVar<T>::after(const NumericVar<T> &e,
                                           const T t) const {
  return e.before(*this, t);
}

template <typename T>
DistanceConstraint<T> NumericVar<T>::before(const NumericVar<T> &e,
                                            const T t) const {
  //  return {e.id(), id(),
  //          (t == Constant::Infinity<T> ? t : e.offset() - offset() - t)};

  //    if(e.id() == Constant::K) {
  //        return this->before(e.offset() - t);
  //    }
  //    if(id() == Constant::K) {
  //        e.after(offset() + t);
  //    }

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
  if (id() == Constant::NoVar) {
    os << "undef";
    if (offset() > 0)
      os << " (+ " << offset() << ")";
    else if (offset() < 0)
      os << " (- " << -offset() << ")";
  } else if (id() == Constant::K) {
    //    assert(sign() == 1);
    os << offset();
  } else {
    //    if (not sign())
    //      os << "-";
    os << "x" << id();
    if (offset() > 0)
      os << " + " << offset();
    else if (offset() < 0)
      os << " - " << -offset();
  }
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

template <typename T>
Interval<T>::Interval(Solver<T> &solver, const NumericVar<T> s,
                      const NumericVar<T> e, const NumericVar<T> d,
                      const BooleanVar<T> o)
    : _id_(num_intervals++), start(s), end(e), duration(d), exist(o) {

  //  std::cout << "\ncreate start\n";
  start.extract(solver);

  //  std::cout << "\ncreate end\n";
  end.extract(solver);

  //  std::cout << "\ncreate duration\n";
  duration.extract(solver);

  exist.extract(solver);

  //
  //  solver.post((start + duration) == end);
  if (start.id() != end.id()) {
    solver.post(start.before(end, duration.min(solver)));
    //
    if (duration.max(solver) != Constant::Infinity<T>)
      solver.post(end.before(start, -duration.max(solver)));
  }
}

template <typename T>
Interval<T>::Interval(Solver<T> &solver, const T mindur, const T maxdur,
                      const T earliest_start, const T latest_start,
                      const T earliest_end, const T latest_end,
                      const BooleanVar<T> opt)
    : _id_(num_intervals++), exist(opt) {

  //                          exist.extract(solver);
  //                          std::cout << "here " << opt << " / " <<
  //                          exist._is_expression << " / " << exist.id() << " /
  //                          " << Constant::NoVar << std::endl;

  //                          if(exist.id() == Constant::NoVar) {
  //                              exist = BooleanVar<T>(Constant::True);
  //                          }

  if (earliest_start == latest_start) {
    start = NumericVar(Constant::K, earliest_start);
  } else {
    start = solver.newNumeric(earliest_start, latest_start);
  }

  if (mindur == maxdur) {
    end = NumericVar(start.id(), mindur);
    duration = NumericVar(Constant::K, mindur);
  } else {
    auto s{start.min(solver)};
    if (s != start.max(solver)) {
      end = solver.newNumeric(earliest_end, latest_end);
      duration = solver.newNumeric(mindur, maxdur);
      solver.post((start + duration) == end);
      solver.post(start.before(end, mindur));

      if (maxdur != Constant::Infinity<T>)
        solver.post(end.before(start, -maxdur));
    } else {
      auto ect{std::max(earliest_end, s + mindur)};
      auto lct{std::min(latest_end, s + maxdur)};
      end = solver.newNumeric(ect, lct);
      duration = NumericVar(end.id(), -s);
    }
  }
}

template <typename T> int Interval<T>::id() const {
  return /*start.id()*/ _id_;
}

template <typename T> bool Interval<T>::operator==(const Interval<T> &t) const {
  return id() == t.id();
}

template <typename T>
template<concepts::distance_provider S>
T Interval<T>::getEarliestStart(const S &solver) const {
  return start.earliest(solver);
}

template <typename T>
template<concepts::distance_provider S>
T Interval<T>::getLatestStart(const S &solver) const {
  return start.latest(solver);
}

template <typename T>
template<concepts::distance_provider S>
T Interval<T>::getEarliestEnd(const S &solver) const {
  return end.earliest(solver);
}

template <typename T>
template<concepts::distance_provider S>
T Interval<T>::getLatestEnd(const S &solver) const {
  return end.latest(solver);
}

template <typename T> bool Interval<T>::mustExist(Solver<T> &) const {
  return true;
}

template <typename T> bool Interval<T>::cannotExist(Solver<T> &) const {
  return false;
}

template <typename T>
template<concepts::distance_provider S>
T Interval<T>::minDuration(const S &solver) const {
  return duration.min(solver);
}

template <typename T>
template<concepts::distance_provider S>
T Interval<T>::maxDuration(const S &solver) const {
  return duration.max(solver);
}

template <typename T> var_t Interval<T>::getStart() const { return start.id(); }

template <typename T> var_t Interval<T>::getEnd() const { return end.id(); }

template <typename T> std::ostream &Interval<T>::display(std::ostream &os) const {
  os << "t" << id(); //<< ": [" << start.earliest(solver) << ".." <<
                     // end.latest(solver) << "]";

  //    auto est{start.earliest(solver)};
  //    auto lst{start.latest(solver)};
  //    auto ect{end.earliest(solver)};
  //    auto lct{end.latest(solver)};
  //    auto pmin{duration.min(solver)};
  //    auto pmax{duration.max(solver)};
  //
  //   os << ": [" << est << "-" << lst << ".." << ect << "-" << lct << "] (" <<
  //   pmin << "-" << pmax << ")";

  os << ": from " << start << " to " << end << " for "
     << duration; //<< "/" << exist;
  return os;
}

template <typename T> std::ostream &operator<<(std::ostream &os, const Interval<T> &x) {
  return x.display(os);
}
}

#endif // __MODEL_HPP
