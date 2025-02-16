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
#include <ranges>
#include <algorithm>
#include <memory>
#include <variant>
#include <Iterators.hpp>

#include "Literal.hpp"
#include "constraints/Cardinality.hpp"
#include "constraints/Incrementality.hpp"
#include "constraints/PseudoBoolean.hpp"
#include "constraints/SumConstraint.hpp"
#include "util/Matrix.hpp"
#include "util/edge_distance.hpp"
#include "util/traits.hpp"

/// @TODO rewritte:
/// ExpressionImpl should be an interface with var_id(); post(solver);
/// extract(solver); BooleanVar and NumericVar should implement Expression
///  BooleanExpression and NumericExpression are pointers to Expression
///  Relational expressions are ExpressionImpl with typically a
///  vector<BooleanExpression> or  a vector<NumericExpression> or both

namespace tempo {
template <typename T> class Solver;

// template <typename T> class BooleanExpression;

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
  ExpressionImpl() = default;

  virtual ~ExpressionImpl() = default;

  ExpressionImpl(const ExpressionImpl &) = default;

  ExpressionImpl(ExpressionImpl &&) = default;

  ExpressionImpl &operator=(const ExpressionImpl &) = default;

  ExpressionImpl &operator=(ExpressionImpl &&) = default;

  // if target is not Constant::NoVar, reference to it instead of creating a
  // new var
  //        virtual var_t extract(Solver<T> &, const var_t target =
  //        Constant::NoVar) = 0;
  virtual var_t extract(Solver<T> &, const var_t target = Constant::NoVar) = 0;

  virtual void post(Solver<T> &) {
    throw ModelingException("This predicate cannot be a constraint");
  }

  virtual std::string name() const { return "some expression"; }

  virtual var_t id() const = 0; //{ return Constant::NoVar; }
  virtual T offset() const { return 0; }
  virtual index_t semantic() const { return Constant::NoSemantic; }

  virtual std::ostream &display(std::ostream &os) const {
    os << "unknown expression";
    return os;
  };
  //    bool extraction_flag{false};
};

template <concepts::scalar T>
using ExpressionPtr = std::shared_ptr<ExpressionImpl<T>>;

template <concepts::scalar T>
using cExpressionPtr = std::shared_ptr<const ExpressionImpl<T>>;

/**
 * @brief Manages either an expression implementation or a variable information
 * struct
 * @tparam T timing type
 * @tparam InfoType info struct type
 */
template <concepts::scalar T, template <concepts::scalar> typename InfoType>
class InfoStorage {
  std::variant<InfoType<T>, ExpressionPtr<T>> content;

public:
  template <typename... Args>
  explicit constexpr InfoStorage(Args &&...args)
      : content(std::forward<Args>(args)...) {}

  /**
   * Whether contains an expression
   * @return
   */
  [[nodiscard]] bool isExpression() const noexcept {
    return std::holds_alternative<ExpressionPtr<T>>(content);
  }

  /**
   * Clears expression implementation (if any) and sets data to information
   * struct
   * @tparam Val information type
   * @param val new information struct
   */
  template <typename Val> void setInfo(Val &&val) {
    content = std::forward<Val>(val);
  }

  /**
   * Marks instance as non-expression with empty info
   */
  void clearExpressionStatus() {
    if (isExpression()) {
      setInfo(InfoType<T>{});
    }
  }

  /**
   * Gets the contained implementation pointer
   * @throws bad_variant_accesss if data is not an expression
   */
  decltype(auto) implem() const noexcept {
    return std::get<ExpressionPtr<T>>(content);
  }

  /**
   * Gets the contained implementation pointer
   * @throws bad_variant_accesss if data is not an expression
   */
  decltype(auto) implem() noexcept {
    return std::get<ExpressionPtr<T>>(content);
  }

  /**
   * Gets the contained info data
   * @throws bad_variant_accesss if data is an expression
   */
  decltype(auto) data() const noexcept {
    return std::get<InfoType<T>>(content);
  }

  /**
   * Gets the contained info data
   * @throws bad_variant_accesss if data is an expression
   */
  decltype(auto) data() noexcept { return std::get<InfoType<T>>(content); }
};

template <typename T = int> struct NumInfo {
  NumInfo() {}
  //  NumInfo(const var_t i, const T o) : _id_(2 * i + 1), _offset(o) {}
  NumInfo(const var_t i, const T o) : _id_(i), _offset(o) {}

  var_t _id_{Constant::NoIndex};
  T _offset{0};
};

/**
 * @brief Wrapper/pointer for numeric variables and expressions
 * @details @copybrief
 * Stores the id of the actual numeric variable, and implements various helper
 * methods
 * @tparam T timing type
 */
template <typename T = int> class NumericVar : public InfoStorage<T, NumInfo> {
public:
  using TimingType = T;
  constexpr NumericVar() noexcept
      : InfoStorage<T, NumInfo>(std::in_place_type<NumInfo<T>>, Constant::NoVar,
                                0) {}

  NumericVar(ExpressionPtr<T> i) noexcept
      : InfoStorage<T, NumInfo>(std::in_place_type<ExpressionPtr<T>>,
                                std::move(i)) {}

  NumericVar(const var_t i, const T o = 0) noexcept
      : InfoStorage<T, NumInfo>(std::in_place_type<NumInfo<T>>, i, o) {}

  template <distance_provider S> T min(const S &sc) const;

  template <distance_provider S> T max(const S &sc) const;

  template <distance_provider S> T earliest(const S &) const;

  template <distance_provider S> T latest(const S &) const;

  Literal<T> after(const T t) const;

  Literal<T> before(const T t) const;

  Literal<T> greaterThanOrEqual(const T t) const;

  Literal<T> lessThanOrEqual(const T t) const;

  DistanceConstraint<T> after(const NumericVar<T> &e, const T t = 0) const;

  DistanceConstraint<T> before(const NumericVar<T> &e, const T t = 0) const;

  DistanceConstraint<T> greaterThanOrEqual(const NumericVar<T> &e,
                                           const T t = 0) const;

  DistanceConstraint<T> lessThanOrEqual(const NumericVar<T> &e,
                                        const T t = 0) const;

  var_t id() const {
    return (this->isExpression() ? this->implem()->id() : this->data()._id_);
  }

  std::ostream &display(std::ostream &os) const;

  static bool isNumeric() { return true; }

  //    void negate() { data._id_ ^= 1; }
  //    void setSign(const bool s) { data._id_ = data._id_ s; }
  void setId(const var_t i) { this->data()._id_ = i; }

  T offset() const {
    return (this->isExpression() ? this->implem()->offset()
                                 : this->data()._offset);
  }

  void setOffset(const T o) { this->data()._offset = o; }

  template <distance_provider S>
  void extract(S &solver, const var_t target = Constant::NoVar) {
#ifdef DBG_EXTRACT
    std::cout << "beg extract " << *this << std::endl;
#endif

    if (this->isExpression()) {
#ifdef DBG_EXTRACT
      std::cout << " (expr)" << std::endl;
#endif

      if (this->implem()->id() == Constant::NoVar) {
#ifdef DBG_EXTRACT
        std::cout << " extract expr" << std::endl;
#endif

        this->implem()->extract(solver, target);

#ifdef DBG_EXTRACT
        std::cout << " ==> " << this->implem()->id() << std::endl;
#endif
      }

#ifdef DBG_EXTRACT
      else
        std::cout << " expr already extracted" << std::endl;
#endif
      this->setInfo(NumInfo<T>(this->implem()->id(), this->implem()->offset()));
    }

#ifdef DBG_EXTRACT
    std::cout << "end extract " << *this << std::endl;
#endif
  }

  template <distance_provider S> void post(S &solver) {
    if (this->isExpression()) {
      throw ModelingException("Numeric expression cannot be constraints");
    } else {
      solver.addToSearch(*this);
    }
  }
};

template <typename T = int> struct BoolInfo {
  constexpr BoolInfo() noexcept = default;

  constexpr BoolInfo(const var_t i, const info_t f) noexcept
      : _id_(i), _edge_id_(f) {}

  constexpr bool operator==(const BoolInfo &) const noexcept = default;

  var_t _id_{Constant::NoIndex};
  index_t _edge_id_{Constant::NoIndex};
};

/**
 * @brief Wrapper/pointer for Boolean variables and expressions
 * @details @copybrief
 * Stores the id of the actual Boolean variable, and implements various helper
 * methods
 * @tparam T timing type
 */
template <typename T = int> class BooleanVar : public InfoStorage<T, BoolInfo> {
public:
  constexpr BooleanVar() noexcept
      : InfoStorage<T, BoolInfo>(std::in_place_type<BoolInfo<T>>,
                                 Constant::NoIndex, Constant::NoSemantic) {}

  BooleanVar(ExpressionPtr<T> i)
      : InfoStorage<T, BoolInfo>(std::in_place_type<ExpressionPtr<T>>,
                                 std::move(i)) {}

  constexpr BooleanVar(const var_t i, const index_t f = 0) noexcept
      : InfoStorage<T, BoolInfo>(std::in_place_type<BoolInfo<T>>, i, f) {}

  constexpr explicit BooleanVar(Literal<T> l)
      : InfoStorage<T, BoolInfo>(std::in_place_type<BoolInfo<T>>, l.variable(),
                                 l.semantic()) {}

  Literal<T> operator==(const bool t) const;

//  constexpr bool operator==(const BooleanVar &other) const {
//      
//      std::cout << "should not do that!";
//      exit(1);
//      
////    if (this->isExpression() != other.isExpression()) {
////      return false;
////    }
////
////    if (this->isExpression()) {
////      return this->implem()->id() == other.implem()->id() and
////             this->implem()->semantic() == other.implem()->semantic();
////    }
////
////    return this->data() == other.data();
//  }

  BooleanVar<T> implies(const BooleanVar<T> x) const;

  var_t id() const {
    return (this->isExpression() ? this->implem()->id() : this->data()._id_);
  }

  index_t semantic() const {
    return (this->isExpression() ? this->implem()->semantic()
                                 : this->data()._edge_id_);
  }

  bool isTrue(Solver<T> &s) const { return s.boolean.isTrue(id()); }
  bool isFalse(Solver<T> &s) const { return s.boolean.isFalse(id()); }

  operator var_t() const { return id(); }

  std::ostream &display(std::ostream &os) const;

  static bool isNumeric() { return false; }

  void setId(const var_t i) { this->data()._id_ = i; }

  template <distance_provider S>
  void extract(S &solver, const var_t target = Constant::NoVar) {
#ifdef DBG_EXTRACT
    std::cout << "beg extract " << *this << std::endl;
#endif

    if (!this->isExpression()) {
      return;
    }
#ifdef DBG_EXTRACT
    std::cout << " (expr)" << std::endl;
#endif

    if (this->implem()->id() == Constant::NoVar) {
#ifdef DBG_EXTRACT
      std::cout << " extract expr" << std::endl;
#endif

      this->implem()->extract(solver, target);

#ifdef DBG_EXTRACT
      std::cout << " ==> " << this->implem()->id() << std::endl;
#endif
    }

#ifdef DBG_EXTRACT
    else
      std::cout << " expr already extracted" << std::endl;
#endif

    this->setInfo(
        BoolInfo<T>(this->implem()->id(), this->implem()->semantic()));

#ifdef DBG_EXTRACT
    std::cout << "end extract " << *this << std::endl;
#endif
  }

  template <distance_provider S> void post(S &solver) {
#ifdef DBG_EXTRACT
    std::cout << "beg post" << std::endl;
#endif

    if (this->isExpression()) {
      if (this->implem()->id() == Constant::NoVar) {
        this->implem()->post(solver);
      }

      this->clearExpressionStatus();
    } else {
      solver.addToSearch(*this);
    }

#ifdef DBG_EXTRACT
    std::cout << "end post" << std::endl;
#endif
  }
};

template <typename T> Literal<T> BooleanVar<T>::operator==(const bool t) const {
  return makeBooleanLiteral<T>(t, id(), semantic());
}

template <typename T = int> class CumulativeExpression;
template <typename T = int> class NoOverlapExpression;

/**
* @brief Wrapper for interval variables
* @details stores
* - the id of a  temporal variable standing for the start
* - the id of a  temporal variable standing for the end
* - the id of a  Boolean variable standing whether the interval actually is in
the schedule
* @tparam T timing type
*/
template <typename T = int> class Interval {
  //    protected:
  //        Interval(NumericVar<T> start, NumericVar<T> end, NumericVar<T>
  //        duration) : _id_(start.id()), start(start), end(end),
  //        duration(duration) {}

public:
  constexpr Interval() noexcept : start(), end(), duration() {}

  //        template<distance_provider S>
  //        Interval(S &solver, const T mindur = 0,
  //                 const T maxdur = Constant::Infinity<T>,
  //                 const T earliest_start = -Constant::Infinity<T>,
  //                 const T latest_start = Constant::Infinity<T>,
  //                 const T earliest_end = -Constant::Infinity<T>,
  //                 const T latest_end = Constant::Infinity<T>,
  //                 const BooleanVar<T> opt = Solver<T>::truism());

  //        template<distance_provider S>
  Interval(/*S &solver,*/ const NumericVar<T> s, const NumericVar<T> e,
           const NumericVar<T> d, const BooleanVar<T> opt = Constant::True);

  template <distance_provider S> T getEarliestStart(const S &s) const;

  template <distance_provider S> T getLatestStart(const S &s) const;

  template <distance_provider S> T getEarliestEnd(const S &s) const;

  template <distance_provider S> T getLatestEnd(const S &s) const;

  bool mustExist(Solver<T> &s) const;

  bool cannotExist(Solver<T> &s) const;

  template <distance_provider S> T minDuration(const S &s) const;

  template <distance_provider S> T maxDuration(const S &s) const;

  var_t getStart() const;

  var_t getEnd() const;

  int id() const;

  bool operator==(const Interval<T> &t) const;

  std::ostream &display(std::ostream &os) const;

  index_t _id_{Constant::NoIndex};

  NumericVar<T> start;
  NumericVar<T> end;
  NumericVar<T> duration;

  void require(NoOverlapExpression<T> &R);

  //    void require(const T d, CumulativeExpression<T>& R);
  void require(const NumericVar<T> d, CumulativeExpression<T> &R);

  //  bool isOptional() const { return exist.id() != Constant::NoVar; }
  bool isOptional(Solver<T> &s) const { return not exist.isTrue(s); }

  void extract(Solver<T> &solver);

  BooleanVar<T> exist{Constant::True};

  static index_t num_intervals;
};

template <typename T> index_t Interval<T>::num_intervals = 0;

template <typename T = int>
class BooleanExpressionImpl : public ExpressionImpl<T> {
public:
  var_t id() const override { return self.id(); }
  index_t semantic() const override { return self.semantic(); }

protected:
  BooleanVar<T> self{Constant::NoVar};
};

template <typename T = int>
class NumericExpressionImpl : public ExpressionImpl<T> {
public:
  var_t id() const override { return self.id(); }
  T offset() const override { return self.offset(); }

protected:
  NumericVar<T> self{Constant::NoVar};
};

template <typename T = int>
class SumExpressionImpl : public NumericExpressionImpl<T> {
public:
  std::string name() const override { return "sum"; }

  void increaseBias(const T v) {
    NumericExpressionImpl<T>::self.setOffset(
        NumericExpressionImpl<T>::self.offset() + v);
  }

  template <distance_provider S> void preprocessNumeric(S &solver) {
#ifdef DBG_EXTRACT_SUM
    std::cout << "preprocess\n";
#endif

    if (numeric_arguments.empty()) {
      return;
    }
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
#ifdef DBG_EXTRACT_SUM
        std::cout << ": null weight => remove\n";
#endif

        // coefficient is 0 -> just ignore
        rm = true;
      } else if (i > 0 and x.id() == args[i - 1].id()) {
#ifdef DBG_EXTRACT_SUM
        std::cout << ": duplicate => change weight from " << ws[i - 1] << " to "
                  << (ws[i - 1] + w) << "  and remove\n";
#endif

        // duplicate -> increase
        ws[i - 1] += w;
        rm = true;
      } else if (x.min(solver) == x.max(solver)) {
        // constant (after extraction) -> move it to the bias
#ifdef DBG_EXTRACT_SUM
        std::cout << ": constant => increase the bias by "
                  << (w * x.min(solver));
#endif
        if (x.id() != Constant::K) // otherwise the bias was already counted
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
        //          increaseBias(-w * x.min(solver));
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

  template <distance_provider S> void preprocessBoolean(S &solver) {
    if (boolean_arguments.empty()) {
      return;
    }
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

  var_t extract(Solver<T> &solver,
                const var_t target = Constant::NoVar) override {
#ifdef DBG_EXTRACT_SUM
    std::cout << "extract ";
    this->display(std::cout);
    std::cout << " with target " << target << std::endl;
#endif

    //      NumericVar<T> num_part;

    preprocessBoolean(solver);

    if (not boolean_arguments.empty()) {
      //        std::cout << bool_lb << ".." << bool_ub << std::endl;
      std::vector<Literal<T>> L;
      for (auto x : boolean_arguments) {
        L.push_back(x == true);
      }

      auto bool_part{target};
      if (target == Constant::NoVar) {
        //            auto bool_part{solver.newNumeric(bool_lb, bool_ub)};
        NumericExpressionImpl<T>::self = solver.newNumeric(bool_lb, bool_ub);
        bool_part = NumericExpressionImpl<T>::self.id();
      } else {
        NumericExpressionImpl<T>::self.data() = NumInfo<T>(bool_part, 0);
        solver.set(leq<T>(target, bool_ub));
        solver.set(geq<T>(target, bool_lb));
      }

      solver.post(new PseudoBooleanLeqVar<T>(solver, L.begin(), L.end(),
                                             bool_weights.begin(), bool_part));

      solver.post(new PseudoBooleanGeqVar<T>(solver, L.begin(), L.end(),
                                             bool_weights.begin(), bool_part));

      assert(numeric_arguments.empty());

      return bool_part;
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
        NumericExpressionImpl<T>::self.clearExpressionStatus();
        //                      num_part.setId(*nvars.begin());
        //                      num_part.setOffset(bias);
        //                      num_part._is_expression = false;
      } else {
        if (target == Constant::NoVar) {
          NumericExpressionImpl<T>::self = solver.newNumeric(num_lb, num_ub);
          //                            num_part = solver.newNumeric(num_lb,
          //                            num_ub);
        } else {
          NumericExpressionImpl<T>::self.data() = NumInfo<T>(target, -bias);
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
        for (auto &w : num_weights) {
          w = -w;
        }

        solver.post(new SumConstraint(solver, nvars.begin(), nvars.end(),
                                      num_weights.begin(), bias));
      }

      //
    } else {
      NumericExpressionImpl<T>::self.setId(Constant::K);
      NumericExpressionImpl<T>::self.clearExpressionStatus();
      //            num_part.setId(Constant::K);
      //            num_part.setOffset(bias);
      //            num_part._is_expression = false;

      //        std::cout << " extract constant sum = " <<
      //        NumericExpressionImpl<T>::self.offset() << std::endl;
    }

    //      NumericExpressionImpl<T>::self = num_part;

    //    std::cout << " => " << NumericExpressionImpl<T>::self << std::endl;

    return NumericExpressionImpl<T>::self.id();
  }

  SumExpressionImpl<T> &addTerm(const NumericVar<T> &x, const T w = 1) {
    //      std::cout << "increase bias by " << w << "*" << x.offset() << ", was
    //      " << NumericExpressionImpl<T>::self.offset();

    increaseBias(w * x.offset());

    //      std::cout << " now: " << NumericExpressionImpl<T>::self.offset() <<
    //      std::endl;

    if (x.id() != Constant::K) {
      numeric_arguments.push_back(x);
      num_weights.push_back(w);
    }

#ifdef DBG_EXTRACT_SUM
    std::cout << "add term " << w << "*" << x
              << " => offset = " << NumericExpressionImpl<T>::self.offset()
              << std::endl;
#endif

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

#ifdef DBG_EXTRACT_SUM
    std::cout << "add term " << k
              << " => offset = " << NumericExpressionImpl<T>::self.offset()
              << std::endl;
#endif

    return *this;
  }

    std::ostream &displayNumArg(std::ostream &os, const unsigned i) const;
  std::ostream &display(std::ostream &os) const override;

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
std::ostream &SumExpressionImpl<T>::displayNumArg(std::ostream &os, const unsigned i) const {
    if(numeric_arguments[i].isExpression()) {
        if (num_weights[i] == 1) {
          os << "+ " << numeric_arguments[i];
        } else if (num_weights[i] == -1) {
          os << "- " << numeric_arguments[i];
        } else if (num_weights[i] > 0) {
          os << "+ " << num_weights[i] << "*" << numeric_arguments[i];
        } else if (num_weights[i] < 0) {
          os << "- " << -num_weights[i] << "*" << numeric_arguments[i];
        }
    } else {
        if (num_weights[i] == 1) {
          os << "+ x" << numeric_arguments[i].id();
        } else if (num_weights[i] == -1) {
          os << "- x" << numeric_arguments[i].id();
        } else if (num_weights[i] > 0) {
          os << "+ " << num_weights[i] << "*x" << numeric_arguments[i].id();
        } else if (num_weights[i] < 0) {
          os << "- " << -num_weights[i] << "*x" << numeric_arguments[i].id();
        }
    }
    
    return os;
}


template <typename T>
std::ostream &SumExpressionImpl<T>::display(std::ostream &os) const {
  //        os << "{";
  //        if (not numeric_arguments.empty()) {
  //            os << num_weights[0] << "*(" << numeric_arguments[0] << ")";
  //            for (unsigned i{1}; i < numeric_arguments.size(); ++i) {
  //                os << " + " << num_weights[i] << "*(" <<
  //                numeric_arguments[i] << ")";
  //            }
  //        }
  //        if (not boolean_arguments.empty()) {
  //            os << bool_weights[0] << "*" << boolean_arguments[0];
  //            for (unsigned i{1}; i < boolean_arguments.size(); ++i) {
  //                os << " + " << bool_weights[i] << "*" <<
  //                boolean_arguments[i];
  //            }
  //        }
  //        os << "} bias=" << NumericExpressionImpl<T>::self.offset();

  os << "(";
  if (not numeric_arguments.empty()) {
      displayNumArg(os, 0);
      
//      if (num_weights[0] == 1) {
//          os << numeric_arguments[0];
//      } else if (num_weights[0] == -1) {
//          os << "-" << numeric_arguments[0];
//      } else {
//          os << num_weights[0] << "*" << numeric_arguments[0];
//      }
    for (unsigned i{1}; i < numeric_arguments.size(); ++i) {
        os << " ";
        displayNumArg(os, i);
//      if (num_weights[i] == 1) {
//        os << " + " << numeric_arguments[i];
//      } else if (num_weights[i] == -1) {
//        os << " - " << numeric_arguments[i];
//      } else if (num_weights[i] > 0) {
//        os << " + " << num_weights[i] << "*" << numeric_arguments[i];
//      } else if (num_weights[i] < 0) {
//        os << " - " << -num_weights[i] << "*" << numeric_arguments[i];
//      }
    }
  }
  if (not boolean_arguments.empty()) {
    os << bool_weights[0] << "*" << boolean_arguments[0];
    for (unsigned i{1}; i < boolean_arguments.size(); ++i) {
      if (bool_weights[i] == 1) {
        os << " + " << boolean_arguments[i];
      } else if (bool_weights[i] == -1) {
        os << " - " << boolean_arguments[i];
      } else if (bool_weights[i] > 0) {
        os << " + " << bool_weights[i] << "*" << boolean_arguments[i];
      } else if (bool_weights[i] < 0) {
        os << " - " << bool_weights[i] << "*" << boolean_arguments[i];
      }
    }
  }
  os << ")";
  if (NumericExpressionImpl<T>::self.offset() > 0)
    os << " + " << NumericExpressionImpl<T>::self.offset();
  else if (NumericExpressionImpl<T>::self.offset() < 0)
    os << " - " << -NumericExpressionImpl<T>::self.offset();

  return os;
}

template <typename T>
NumericVar<T> operator+(const NumericVar<T> &x, const T k) {
#ifdef DBG_EXTRACT_SUM
  std::cout << " create expression " << x << " + " << k << std::endl;
#endif

  auto sum = std::make_shared<SumExpressionImpl<T>>();
  sum->addTerm(x);
  sum->addTerm(k);
  NumericVar<T> exp(std::move(sum));
  return exp;
}

template <typename T>
NumericVar<T> operator+(const NumericVar<T> &x, const NumericVar<T> &y) {
#ifdef DBG_EXTRACT_SUM
  std::cout << " create expression " << x << " + " << y << std::endl;
#endif

  auto sum = std::make_shared<SumExpressionImpl<T>>();
  sum->addTerm(x);
  sum->addTerm(y);
  NumericVar<T> exp(std::move(sum));
  return exp;
}

template <typename T>
NumericVar<T> operator-(const NumericVar<T> &x, const T k) {
#ifdef DBG_EXTRACT_SUM
  std::cout << " create expression " << x << " - " << k << std::endl;
#endif

//  auto sum = std::make_shared<SumExpressionImpl<T>>();
//  sum->addTerm(x);
//  sum->addTerm(-k);
//  NumericVar<T> exp(std::move(sum));
    
    
  return NumericVar<T>(x.id(), x.offset+k);
}

template <typename T>
NumericVar<T> operator-(const NumericVar<T> &x, const NumericVar<T> &y) {
#ifdef DBG_EXTRACT_SUM
  std::cout << " create expression " << x << " - " << y << std::endl;
#endif

  auto sum = std::make_shared<SumExpressionImpl<T>>();
  sum->addTerm(x);
  sum->addTerm(y, -1);
  NumericVar<T> exp(std::move(sum));
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
  auto sum = std::make_shared<SumExpressionImpl<T>>();
  auto wp{W.begin()};
  for (auto xp{X.begin()}; xp != X.end(); ++xp) {
    sum->addTerm(*xp, *wp);
    ++wp;
  }
  NumericVar<T> exp(std::move(sum));
  return exp;
}

template <typename T = int>
NumericVar<T> Sum(const std::vector<NumericVar<T>> &X) {
  auto sum = std::make_shared<SumExpressionImpl<T>>();
  for (auto xp{X.begin()}; xp != X.end(); ++xp) {
    sum->addTerm(*xp, 1);
  }
  NumericVar<T> exp(std::move(sum));
  return exp;
}

template <typename T = int>
NumericVar<T> Sum(const std::vector<BooleanVar<T>> &X,
                  const std::vector<T> &W) {
  auto sum = std::make_shared<SumExpressionImpl<T>>();
  auto wp{W.begin()};
  for (auto xp{X.begin()}; xp != X.end(); ++xp) {
    sum->addTerm(*xp, *wp);
    ++wp;
  }
  NumericVar<T> exp(std::move(sum));
  return exp;
}

template <typename T = int>
NumericVar<T> Sum(const std::vector<NumericVar<T>> &X,
                  const std::vector<BooleanVar<T>> &B,
                  const std::vector<T> &W) {
  auto sum = std::make_shared<SumExpressionImpl<T>>();
  auto wp{W.begin()};
  for (auto xp{X.begin()}; xp != X.end(); ++xp) {
    sum->addTerm(*xp, *wp);
    ++wp;
  }
  for (auto xp{B.begin()}; xp != B.end(); ++xp) {
    sum->addTerm(*xp, *wp);
    ++wp;
  }
  NumericVar<T> exp(std::move(sum));
  return exp;
}

// x == y
template <typename T = int>
class BoolEqExpressionImpl : public BooleanExpressionImpl<T> {
public:
  BoolEqExpressionImpl(BooleanVar<T> x, BooleanVar<T> y)
      : x(x), y(y) {}

  virtual std::ostream &display(std::ostream &os) const override {
    os << x << " <=> " << y;
    return os;
  };

  std::string name() const override { return "eq"; }

  var_t extract(Solver<T> &solver, const var_t = Constant::NoVar) override {
    x.extract(solver);
    y.extract(solver);
    BooleanExpressionImpl<T>::self = solver.newBoolean();
      
      auto& e{BooleanExpressionImpl<T>::self};
      
      std::vector<Literal<T>> cl1{x == true, x == false, e == false};
      solver.clauses.add(cl1);
      std::vector<Literal<T>> cl2{x == false, x == true, e == false};
      solver.clauses.add(cl2);
      std::vector<Literal<T>> cl3{x == true, x == true, e == true};
      solver.clauses.add(cl3);
      std::vector<Literal<T>> cl4{x == false, x == false, e == true};
      solver.clauses.add(cl4);
      
    return BooleanExpressionImpl<T>::self.id();
  }

  void post(Solver<T> &solver) override {
    if (x.isExpression()) {
      y.extract(solver);
      x.extract(solver, y.id());
    } else if (y.isExpression()) {
      x.extract(solver);
      y.extract(solver, x.id());
    } else {
      x.extract(solver);
      y.extract(solver);
        std::vector<Literal<T>> cl1{x == true, x == false};
        solver.clauses.add(cl1);
        std::vector<Literal<T>> cl2{x == false, x == true};
        solver.clauses.add(cl2);
    }
  }

private:
    BooleanVar<T> x;
    BooleanVar<T> y;
};

// x == y+k
template <typename T = int>
class NumEqExpressionImpl : public BooleanExpressionImpl<T> {
public:
  NumEqExpressionImpl(NumericVar<T> x, NumericVar<T> y, const T k)
      : x(x), y(y), k(k) {}

  virtual std::ostream &display(std::ostream &os) const override {
    os << x << " = " << y;
    if (k > 0)
      os << " + " << k;
    else if (k < 0)
      os << " - " << -k;
    return os;
  };

  std::string name() const override { return "eq"; }

  var_t extract(Solver<T> &solver, const var_t = Constant::NoVar) override {
    x.extract(solver);
    y.extract(solver);

    auto prec{x.before(y, -k)};
    auto inf = solver.newDisjunct(~prec, prec);
    auto succ{x.after(y, -k)};
    auto sup = solver.newDisjunct(~succ, succ);
    auto conj = (sup and inf);

    conj.extract(solver);
    BooleanExpressionImpl<T>::self = conj;

    return BooleanExpressionImpl<T>::self.id();
  }

  void post(Solver<T> &solver) override {
    if (x.isExpression()) {
      y.extract(solver);
      x.extract(solver, y.id());
    } else if (y.isExpression()) {
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
  //  BooleanVar<T> self;
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
  BooleanVar<T> exp(std::make_shared<NumEqExpressionImpl<T>>(x, y, 0));
  return exp;
}

template <typename T>
BooleanVar<T> operator==(const BooleanVar<T> &x, const BooleanVar<T> &y) {
  BooleanVar<T> exp(std::make_shared<BoolEqExpressionImpl<T>>(x, y));
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

template <typename T = int>
class LeqExpressionImpl : public BooleanExpressionImpl<T> {
public:
  LeqExpressionImpl(NumericVar<T> x, NumericVar<T> y, const T k)
      : x(x), y(y), k(k) {
    //
    //    std::cout << "leq expr " << x << " + " << x.offset() << " + " << k
    //              << " <= " << y << " + " << y.offset() << std::endl;
    //    std::cout << "leq expr " << x << " + " << x.offset() << " + " << k
    //              << " <= " << y << " + " << y.offset() << std::endl;
  }

  virtual std::ostream &display(std::ostream &os) const override {
    os << x;
    if (k > 0)
      os << " + " << k;
    else if (k < 0)
      os << " - " << -k;
    os << " <= " << y;
    return os;
  };

  std::string name() const override { return "leq"; }

  var_t extract(Solver<T> &solver, const var_t = Constant::NoVar) override {
    //      std::cout << "extract x <= y expression\n - extract x\n";

    x.extract(solver);

    //      std::cout << " -> " << x << "\n - extract y\n";

    y.extract(solver);

    //      std::cout << " -> " << y << "\n create precedence constraint\n";

    auto prec{x.before(y, -k)};

    //      std::cout << " -> " << prec << "\n create disjunct\n";

    BooleanExpressionImpl<T>::self = solver.newDisjunct(~prec, prec);

    //      std::cout << " -> " << solver.pretty(BooleanExpressionImpl<T>::self
    //      == true) << " <> " << solver.pretty(BooleanExpressionImpl<T>::self
    //      == false) << " (" << BooleanExpressionImpl<T>::self.id() << ")\n";

    return BooleanExpressionImpl<T>::self.id();
  }

  void post(Solver<T> &solver) override {
    x.extract(solver);
    y.extract(solver);

    auto prec{x.before(y, -k)};
    solver.post(prec);
  }

private:
  //  BooleanVar<T> self;
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
  BooleanVar<T> exp(std::make_shared<LeqExpressionImpl<T>>(x, y, 0));
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
  BooleanVar<T> exp(
      std::make_shared<LeqExpressionImpl<T>>(x, y, -Gap<T>::epsilon()));
  return exp;
}

template <typename T>
BooleanVar<T> operator<(const T x, const NumericVar<T> &y) {
  BooleanVar<T> exp(std::make_shared<LeqExpressionImpl<T>>(
      NumericVar<T>(Constant::K, x), y, -Gap<T>::epsilon()));
  return exp;
}

template <typename T>
BooleanVar<T> operator<(const NumericVar<T> &x, const T y) {
  BooleanVar<T> exp(std::make_shared<LeqExpressionImpl<T>>(
      x, NumericVar<T>(Constant::K, y), -Gap<T>::epsilon()));
  return exp;
}

///// >=
template <typename T>
BooleanVar<T> operator>=(const NumericVar<T> &x, const NumericVar<T> &y) {
  BooleanVar<T> exp(std::make_shared<LeqExpressionImpl<T>>(y, x, 0));
  return exp;
}

template <typename T>
BooleanVar<T> operator>=(const T x, const NumericVar<T> &y) {
  BooleanVar<T> exp(std::make_shared<LeqExpressionImpl<T>>(
      y, NumericVar<T>(Constant::K, x), 0));
  return exp;
}

template <typename T>
BooleanVar<T> operator>=(const NumericVar<T> &x, const T y) {
  BooleanVar<T> exp(std::make_shared<LeqExpressionImpl<T>>(
      NumericVar<T>(Constant::K, y), x, 0));
  return exp;
}

///// >
template <typename T>
BooleanVar<T> operator>(const NumericVar<T> &x, const NumericVar<T> &y) {
  BooleanVar<T> exp(
      std::make_shared<LeqExpressionImpl<T>>(y, x, -Gap<T>::epsilon()));
  return exp;
}

template <typename T>
BooleanVar<T> operator>(const NumericVar<T> &x, const T k) {
  BooleanVar<T> exp(std::make_shared<LeqExpressionImpl<T>>(
      NumericVar<T>(Constant::K, k + Gap<T>::epsilon()), x, 0));
  return exp;
}

template <typename T>
BooleanVar<T> operator>(const T k, const NumericVar<T> &x) {
  BooleanVar<T> exp(std::make_shared<LeqExpressionImpl<T>>(
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

  virtual std::ostream &display(std::ostream &os) const override {
    os << "(" << *(boolean_arguments.begin());
    for (auto x{boolean_arguments.begin() + 1}; x != boolean_arguments.end();
         ++x)
      os << " AND " << *x;
    os << ")";
    return os;
  };

  std::string name() const override { return "and"; }

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
  BooleanVar<T> exp(
      std::make_shared<LogicalAndExpressionImpl<T>>(sc.begin(), sc.end()));
  return exp;
}

template <concepts::ttyped_range<BooleanVar> Iterable>
auto BigAnd(Iterable &X) {
  using T = typename std::ranges::range_value_t<Iterable>::TimingType;
  BooleanVar<T> exp(
      std::make_shared<LogicalAndExpressionImpl<T>>(X.begin(), X.end()));
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

  virtual std::ostream &display(std::ostream &os) const override {
    os << "(" << *(boolean_arguments.begin());
    for (auto x{boolean_arguments.begin() + 1}; x != boolean_arguments.end();
         ++x)
      os << " OR " << *x;
    os << ")";
    return os;
  };

  std::string name() const override { return "or"; }

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
BooleanVar<T> operator||(const BooleanVar<T> &x, const BooleanVar<T> &y) {
  std::vector<BooleanVar<T>> sc{x, y};
  BooleanVar<T> exp(
      std::make_shared<LogicalOrExpressionImpl<T>>(sc.begin(), sc.end()));
  return exp;
}

template <typename T> BooleanVar<T> BigOr(const std::vector<BooleanVar<T>> &X) {
  BooleanVar<T> exp(
      std::make_shared<LogicalOrExpressionImpl<T>>(X.begin(), X.end()));
  return exp;
}

template <typename T = int>
class LogicalImplicationExpression : public BooleanExpressionImpl<T> {
public:
  LogicalImplicationExpression(BooleanVar<T> x, BooleanVar<T> y)
      : implicant(x), implicate(y) {}

  virtual std::ostream &display(std::ostream &os) const override {
    os << "(" << implicant << " implies " << implicate << ")";
    return os;
  };

  std::string name() const override { return "implication"; }

  var_t extract(Solver<T> &solver, const var_t = Constant::NoVar) override {
    implicant.extract(solver);
    implicate.extract(solver);
    BooleanExpressionImpl<T>::self = solver.newBoolean();

    std::vector<Literal<T>> cl{
        solver.boolean.getLiteral(false, implicant),
        solver.boolean.getLiteral(true, implicate),
        solver.boolean.getLiteral(false, BooleanExpressionImpl<T>::self)};
    solver.clauses.add(cl.begin(), cl.end());

    auto x{BooleanExpressionImpl<T>::self};
    cl = {x == true, implicant == true};
    solver.clauses.add(cl.begin(), cl.end());

    cl = {x == true, implicate == false};
    solver.clauses.add(cl.begin(), cl.end());

    return BooleanExpressionImpl<T>::self.id();
  }

  void post(Solver<T> &solver) override {
    //      std::cout << "post implication\n - extract implicant\n";

    implicant.extract(solver);

    //      std::cout << " -> " << implicant << "\n - extract implicate\n";

    implicate.extract(solver);

    //      std::cout << " -> " << implicate << std::endl;

    std::vector<Literal<T>> cl{solver.boolean.getLiteral(false, implicant),
                               solver.boolean.getLiteral(true, implicate)};

    //      std::cout << "add clause\n";

    solver.clauses.add(cl.begin(), cl.end());

    //      std::cout << "done\n";
  }

private:
  BooleanVar<T> implicant;
  BooleanVar<T> implicate;
};

template <typename T>
BooleanVar<T> BooleanVar<T>::implies(const BooleanVar<T> x) const {
  BooleanVar<T> exp(
      std::make_shared<LogicalImplicationExpression<T>>(*this, x));
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

  virtual std::ostream &display(std::ostream &os) const override {
    os << "(" << *(boolean_arguments.begin());
    for (auto x{boolean_arguments.begin() + 1}; x != boolean_arguments.end();
         ++x)
      os << " + " << *x;
    os << ") in [" << lb << ", " << ub << "]";
    return os;
  };

  std::string name() const override { return "cardinality"; }

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
  NumericVar<T> exp(
      std::make_shared<CardinalityExpressionImpl<T>>(X.begin(), X.end()));
  return exp;
}

template <typename T>
BooleanVar<T> Cardinality(const std::vector<BooleanVar<T>> &X, const T l, const T u) {
  BooleanVar<T> exp(
      std::make_shared<CardinalityExpressionImpl<T>>(X.begin(), X.end(), l, u));
  return exp;
}

template <typename T = int>
class AtMostExpressionImpl : public BooleanExpressionImpl<T> {
public:
  template <typename Iter>
  AtMostExpressionImpl(Iter beg_var, Iter end_var, const T k, const bool sign)
      : k(k), sign(sign) {
    for (auto x{beg_var}; x != end_var; ++x) {
      boolean_arguments.emplace_back(*x);
    }
  }

  virtual std::ostream &display(std::ostream &os) const override {
    os << "(" << *(boolean_arguments.begin());
    for (auto x{boolean_arguments.begin() + 1}; x != boolean_arguments.end();
         ++x)
      os << " + " << *x;
    os << ") <= " << k;
    return os;
  };

  std::string name() const override { return "at-most"; }

  var_t extract(Solver<T> &, const var_t) override {
    throw ModelingException("AtMost is not a predicate");
    return Constant::NoVar;
  }

  void post(Solver<T> &solver) override {
    solver.postCardinality(boolean_arguments.begin(), boolean_arguments.end(),
                           sign, k);
  }

private:
  T k;
  bool sign;
  std::vector<BooleanVar<T>> boolean_arguments;
};

template <typename T>
BooleanVar<T> AtMost(const T k, const std::vector<BooleanVar<T>> &X) {
  BooleanVar<T> exp(
      std::make_shared<AtMostExpressionImpl<T>>(X.begin(), X.end(), k, true));
  return exp;
}

template <typename T>
BooleanVar<T> AtLeast(const T k, const std::vector<BooleanVar<T>> &X) {
  BooleanVar<T> exp(std::make_shared<AtMostExpressionImpl<T>>(
      X.begin(), X.end(), X.size() - k, false));
  return exp;
}

template <typename T = int>
class NoOverlapExpressionImpl : public BooleanExpressionImpl<T>,
                                public std::vector<Interval<T>> {
public:
  NoOverlapExpressionImpl(const Interval<T> &sched) : schedule(sched) {}

  std::string name() const override { return "no-overlap"; }

  virtual std::ostream &display(std::ostream &os) const override {
    os << "NoOverlap(" << schedule;
    for (auto x{this->begin()}; x != this->end(); ++x)
      os << ", " << *x;
    os << ")";
    return os;
  };

  var_t extract(Solver<T> &, const var_t) override {
    throw ModelingException("NoOverlap is not a predicate");
    return Constant::NoVar;
  }

  void post(Solver<T> &solver) override {
    using iterators::const_enumerate;
    using namespace std::views;

    for (auto &inter : *this) {
      inter.extract(solver);
    }

    bool opt_flag{false};
    disjunctiveLiterals.fill(this->size(), this->size(), Contradiction<T>);
    for (auto [intervalIdxA, intervalA] : const_enumerate(*this)) {
      //      if(a->isOptional(solver))
      //          std::cout << " ===> " << a->id() << ": " << a->exist <<
      //          std::endl;

      for (auto [intervalIdxB, intervalB] :
           const_enumerate(*this | drop(intervalIdxA + 1), intervalIdxA + 1)) {
        auto t_ab{0};
        auto t_ba{0};
        if (intervalIdxA < transitions.size() and
            intervalIdxB < transitions[intervalIdxA].size()) {
          t_ab = transitions[intervalIdxA][intervalIdxB];
        }
        if (intervalIdxB < transitions.size() and
            intervalIdxA < transitions[intervalIdxB].size()) {
          t_ba = transitions[intervalIdxB][intervalIdxA];
        }

        auto a_before_b{intervalA.end.before(intervalB.start, t_ab)};
        auto b_before_a{intervalB.end.before(intervalA.start, t_ba)};

        if (intervalA.isOptional(solver) or intervalB.isOptional(solver)) {
          auto x_ab{solver.newDisjunct(Constant::NoEdge<T>, a_before_b)};
          auto x_ba{solver.newDisjunct(Constant::NoEdge<T>, b_before_a)};

          // (intervalA.exist and intervalB.exist) -> (x_ab or x_ba)
          // ~intervalA.exist or ~intervalB.exist or x_ab or x_ba
          std::vector<Literal<T>> cl{x_ab == true, x_ba == true};
          if (intervalA.isOptional(solver))
            cl.push_back(intervalA.exist == false);
          if (intervalB.isOptional(solver))
            cl.push_back(intervalB.exist == false);
          solver.clauses.add(cl.begin(), cl.end());

          // to avoid setting useless constraints
          // (~intervalA.exist -> ~x_ab) and (~intervalA.exist -> ~x_ba)
          // intervalA.exist or ~x_ab
          cl.clear();
          if (intervalA.isOptional(solver)) {
            cl.push_back(intervalA.exist == true);
            cl.push_back(x_ab == false);
            solver.clauses.add(cl.begin(), cl.end());

            cl.pop_back();
            cl.push_back(x_ba == false);
            solver.clauses.add(cl.begin(), cl.end());
          }

          cl.clear();
          if (intervalB.isOptional(solver)) {
            cl.push_back(intervalB.exist == true);
            cl.push_back(x_ab == false);
            solver.clauses.add(cl.begin(), cl.end());

            cl.pop_back();
            cl.push_back(x_ba == false);
            solver.clauses.add(cl.begin(), cl.end());
          }

          // to enforce the disjunction
          cl.clear();
          cl.push_back(x_ab == false);
          cl.push_back(x_ba == false);
          solver.clauses.add(cl.begin(), cl.end());

          disjunctiveLiterals(intervalIdxA, intervalIdxB) =
              solver.boolean.getLiteral(true, x_ab);
          disjunctiveLiterals(intervalIdxB, intervalIdxA) =
              solver.boolean.getLiteral(true, x_ba);
          solver.addToSearch(x_ab);
          solver.addToSearch(x_ba);
          opt_flag = true;
        } else {
          auto var = solver.newDisjunct(a_before_b, b_before_a);
          disjunctiveLiterals(intervalIdxA, intervalIdxB) =
              solver.boolean.getLiteral(true, var);
          disjunctiveLiterals(intervalIdxB, intervalIdxA) =
              solver.boolean.getLiteral(false, var);
          solver.addToSearch(var);
        }
      }
    }

    if (std::distance(this->begin(), this->end()) > 1) {
      if (solver.getOptions().edge_finding) {
        if (opt_flag) {
          std::cout
              << "edge-finding for optional intervals is not implemented, "
                 "please use '--no-edge-finding'. Aborting\n";
          exit(0);
        }

        solver.postEdgeFinding(
            schedule, std::ranges::subrange(this->begin(), this->end()),
            disjunctiveLiterals);
      }

      if (solver.getOptions().transitivity) {
        if (opt_flag) {
          std::cout
              << "transitivity for optional intervals is not implemented, "
                 "please use '--no-transitivity'. Aborting\n";
          exit(0);
        }

        solver.postTransitivity(
            schedule, std::ranges::subrange(this->begin(), this->end()),
            disjunctiveLiterals);
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

  const auto &getDisjunctiveLiterals() const noexcept {
    return disjunctiveLiterals;
  }

private:
  Interval<T> schedule;
  Matrix<Literal<T>> disjunctiveLiterals;
  //  std::vector<BooleanVar<T>> relevant;
  std::vector<std::vector<T>> transitions;
};

template <concepts::scalar T>
using NoOverlapExpressionPtr = std::shared_ptr<NoOverlapExpressionImpl<T>>;

template <concepts::scalar T>
using cNoOverlapExpressionPtr =
    std::shared_ptr<const NoOverlapExpressionImpl<T>>;

template <typename T> class NoOverlapExpression : public BooleanVar<T> {
public:
  NoOverlapExpression(NoOverlapExpressionPtr<T> i)
      : BooleanVar<T>(std::move(i)) {}

  [[nodiscard]] size_t size() const {
    return static_cast<NoOverlapExpressionImpl<T> *>(
               BooleanVar<T>::implem().get())
        ->size();
  }

  std::vector<Interval<T>>::iterator begin() {
    return static_cast<NoOverlapExpressionImpl<T> *>(
               BooleanVar<T>::implem().get())
        ->begin();
  }

  std::vector<Interval<T>>::const_iterator begin() const {
    return const_cast<NoOverlapExpression *>(this)->begin();
  }

  std::vector<Interval<T>>::iterator end() {
    return static_cast<NoOverlapExpressionImpl<T> *>(
               BooleanVar<T>::implem().get())
        ->end();
  }

  std::vector<Interval<T>>::const_iterator end() const {
    return const_cast<NoOverlapExpression *>(this)->end();
  }

  const auto &getDisjunctiveLiterals() const noexcept {
    return static_cast<NoOverlapExpressionImpl<T> *>(
               BooleanVar<T>::implem().get())
        ->getDisjunctiveLiterals();
  }

  void provide(const Interval<T> &i) {
    static_cast<NoOverlapExpressionImpl<T> *>(BooleanVar<T>::implem().get())
        ->push_back(i);
  }

  template <typename Matrix> void setTransisions(const Matrix &D) {
    static_cast<NoOverlapExpressionImpl<T> *>(BooleanVar<T>::implem().get())
        ->setTransisions(D);
  }
};

template <typename T, std::ranges::range Iterable, typename Matrix>
NoOverlapExpression<T> NoOverlap(const Interval<T> &schedule, const Iterable &X,
                                 const Matrix &D) {
  NoOverlapExpression<T> exp(
      std::make_shared<NoOverlapExpressionImpl<T>>(schedule));
  for (auto x : X) {
    exp.provide(x);
  }
  exp.setTransisions(D);
  return exp;
}

template <typename T>
NoOverlapExpression<T> NoOverlap(Interval<T> &schedule,
                                 const std::vector<Interval<T>> &X) {
  NoOverlapExpression<T> exp(
      std::make_shared<NoOverlapExpressionImpl<T>>(schedule));
  for (auto x : X) {
    exp.provide(x);
  }
  return exp;
}

template <typename T> NoOverlapExpression<T> NoOverlap(Interval<T> &schedule) {
  NoOverlapExpression<T> exp(
      std::make_shared<NoOverlapExpressionImpl<T>>(schedule));
  return exp;
}

/*!
Cumulative Resource
*/

template <typename T = int>
class CumulativeExpressionImpl : public BooleanExpressionImpl<T>,
                                 public std::vector<Interval<T>> {
public:
  CumulativeExpressionImpl(const Interval<T> &s, const NumericVar<T> c)
      : schedule(s), capacity(c) {}

  std::string name() const override { return "cumulative"; }

  var_t extract(Solver<T> &, const var_t) override {
    throw ModelingException("Cumulative is not a predicate");
    return Constant::NoVar;
  }

  void provide(const NumericVar<T> d, const Interval<T> i) {
    this->push_back(i);
    demand.push_back(d);
  }

  void post(Solver<T> &solver) override {
    using iterators::const_zip_enumerate;
    using std::views::drop;

    bool opt_flag{false};
    disjunctiveLiterals.fill(this->size(), this->size(), Contradiction<T>);
    for (auto [taskIdxA, intervalA, demandA] :
         const_zip_enumerate(*this, demand)) {
      for (auto [taskIdxB, intervalB, demandB] :
           const_zip_enumerate(*this | drop(taskIdxA + 1),
                               demand | drop(taskIdxA + 1), taskIdxA + 1)) {
        auto ea_before_sb{intervalA.end.before(intervalB.start)};
        auto eb_before_sa{intervalB.end.before(intervalA.start)};

        if (intervalA.isOptional(solver) or intervalB.isOptional(solver)) {
          opt_flag = true;

          auto x_ab{solver.newDisjunct(Constant::NoEdge<T>, ea_before_sb)};
          auto x_ba{solver.newDisjunct(Constant::NoEdge<T>, eb_before_sa)};
          auto x_anb{solver.newDisjunct(Constant::NoEdge<T>, ~ea_before_sb)};
          auto x_bna{solver.newDisjunct(Constant::NoEdge<T>, ~eb_before_sa)};

          disjunctiveLiterals(taskIdxA, taskIdxB) =
              solver.boolean.getLiteral(true, x_anb);
          // a is not before b
          disjunctiveLiterals(taskIdxB, taskIdxA) =
              solver.boolean.getLiteral(true, x_bna);
          // b is not before a
          solver.addToSearch(x_ab);
          solver.addToSearch(x_ba);
          solver.addToSearch(x_anb);
          solver.addToSearch(x_bna);

          // if a and b exist, then either a is not before b or b is not before
          // a
          std::vector<Literal<T>> cl{x_anb == true, x_bna == false};
          if (intervalA.isOptional(solver))
            cl.push_back(intervalA.exist == false);
          if (intervalB.isOptional(solver))
            cl.push_back(intervalB.exist == false);
          solver.clauses.add(cl.begin(), cl.end());

          // to avoid setting useless constraints
          // (~intervalA.exist -> ( ~x_ab and ~x_ba and ~x_anb and ~x_bna)
          cl.clear();
          if (intervalA.isOptional(solver)) {
            cl.push_back(intervalA.exist == true);
            cl.push_back(x_ab == false);
            solver.clauses.add(cl.begin(), cl.end());

            cl.pop_back();
            cl.push_back(x_ba == false);
            solver.clauses.add(cl.begin(), cl.end());

            cl.pop_back();
            cl.push_back(x_anb == false);
            solver.clauses.add(cl.begin(), cl.end());

            cl.pop_back();
            cl.push_back(x_bna == false);
            solver.clauses.add(cl.begin(), cl.end());
          }

          cl.clear();
          if (intervalB.isOptional(solver)) {
            cl.push_back(intervalB.exist == true);
            cl.push_back(x_ab == false);
            solver.clauses.add(cl.begin(), cl.end());

            cl.pop_back();
            cl.push_back(x_ba == false);
            solver.clauses.add(cl.begin(), cl.end());

            cl.pop_back();
            cl.push_back(x_anb == false);
            solver.clauses.add(cl.begin(), cl.end());

            cl.pop_back();
            cl.push_back(x_bna == false);
            solver.clauses.add(cl.begin(), cl.end());
          }

          // to enforce the disjunction // is it necessary ?
          cl.clear();
          cl.push_back(x_ab == false);
          cl.push_back(x_ba == false);
          solver.clauses.add(cl.begin(), cl.end());

          // if a is before b, then a cannot not be before b
          cl.clear();
          cl.push_back(x_ab == false);
          cl.push_back(x_bna == false);
          solver.clauses.add(cl.begin(), cl.end());
          cl.clear();
          cl.push_back(x_ba == false);
          cl.push_back(x_anb == false);
          solver.clauses.add(cl.begin(), cl.end());
        } else {
          if (demandA.min(solver) + demandB.min(solver) >
              capacity.max(solver)) {
            auto d{solver.newDisjunct(eb_before_sa, ea_before_sb)};
            disjunctiveLiterals(taskIdxA, taskIdxB) =
                solver.boolean.getLiteral(false, d);
            // a is not before b (<-> b is before a)
            disjunctiveLiterals(taskIdxB, taskIdxA) =
                solver.boolean.getLiteral(true, d);
            // b is not before a (<-> a is before b)
            solver.addToSearch(d);
          } else {
            auto ba{solver.newDisjunct(~ea_before_sb, ea_before_sb)};
            auto bb{solver.newDisjunct(eb_before_sa, ~eb_before_sa)};
            disjunctiveLiterals(taskIdxA, taskIdxB) =
                solver.boolean.getLiteral(false, ba);
            // a is not before b
            disjunctiveLiterals(taskIdxB, taskIdxA) =
                solver.boolean.getLiteral(true, bb);
            // b is not before a
            solver.addToSearch(ba);
            solver.addToSearch(bb);

            // if a ends before b starts then b cannot end before a
            // starts
            std::vector<Literal<T>> cl{disjunctiveLiterals(taskIdxA, taskIdxB),
                                       disjunctiveLiterals(taskIdxB, taskIdxA)};
            solver.clauses.add(cl.begin(), cl.end());
          }
        }
      }
    }

    if (std::distance(this->begin(), this->end()) > 1) {
      if (opt_flag) {
        std::cout << "cumulative for optional intervals is not implemented. "
                     "Aborting\n";
        exit(0);
      }

      solver.postCumulative(capacity, *this, demand, disjunctiveLiterals);

      //        solver.postCumulativeIncrementality(*this,
      //        disjunctiveLiterals);

      // #ifdef DBG_SEF
      if (solver.getOptions().edge_finding) {
        auto incr = std::make_unique<Incrementality<T>>(solver, *this,
                                                        disjunctiveLiterals);

        solver.postStrongEdgeFinding(
            schedule, capacity, this->begin(), this->end(), this->begDemand(),
            solver.getOptions().tt_edge_finding, incr.get(),
            solver.getOptions().incomplete_edge_finding);

        solver.post(incr.release());
      }
      // #endif

      if (solver.getOptions().overlap_finding) {
        solver.postOverlapFinding(schedule, capacity, this->begin(),
                                  this->end(), this->begDemand(),
                                  disjunctiveLiterals);
      }

      if (solver.getOptions().time_tabling) {
        solver.postTimetabling(capacity, this->begin(), this->end(),
                               this->begDemand());

        //          solver.postTimetablingFixedDemand(capacity, this->begin(),
        //          this->end(),
        //                                 this->begDemand());
      }
    }
  }

  const auto &getDisjunctiveLiterals() const noexcept {
    return disjunctiveLiterals;
  }

  std::vector<NumericVar<T>>::iterator begDemand() { return demand.begin(); }
  std::vector<NumericVar<T>>::const_iterator begDemand() const {
    return demand.begin();
  }
  std::vector<NumericVar<T>>::iterator endDemand() { return demand.end(); }
  std::vector<NumericVar<T>>::const_iterator endDemand() const {
    return demand.end();
  }

private:
  Interval<T> schedule;
  NumericVar<T> capacity;
  Matrix<Literal<T>> disjunctiveLiterals;
  std::vector<NumericVar<T>> demand;
};

template <concepts::scalar T>
using CumulativeExpressionPtr = std::shared_ptr<CumulativeExpressionImpl<T>>;

template <concepts::scalar T>
using cCumulativeExpressionPtr =
    std::shared_ptr<const CumulativeExpressionImpl<T>>;

template <typename T> class CumulativeExpression : public BooleanVar<T> {
public:
  CumulativeExpression(CumulativeExpressionPtr<T> i)
      : BooleanVar<T>(std::move(i)) {}

  std::vector<Interval<T>>::iterator begin() {
    return static_cast<CumulativeExpressionImpl<T> *>(
               BooleanVar<T>::implem().get())
        ->begin();
  }

  std::vector<Interval<T>>::const_iterator begin() const {
    return const_cast<CumulativeExpression *>(this)->begin();
  }

  std::vector<Interval<T>>::iterator end() {
    return static_cast<CumulativeExpressionImpl<T> *>(
               BooleanVar<T>::implem().get())
        ->end();
  }

  std::vector<Interval<T>>::const_iterator end() const {
    return const_cast<CumulativeExpression *>(this)->end();
  }

  const auto &getDisjunctiveLiterals() const noexcept {
    return static_cast<const CumulativeExpressionImpl<T> *>(
               BooleanVar<T>::implem().get())
        ->getDisjunctiveLiterals();
  }

  std::vector<BooleanVar<T>>::iterator begDemand() {
    return static_cast<CumulativeExpressionImpl<T> *>(
               BooleanVar<T>::implem().get())
        ->begDemand();
  }

  std::vector<BooleanVar<T>>::const_iterator begDemand() const {
    return const_cast<CumulativeExpression *>(this)->begDemand();
  }

  std::vector<BooleanVar<T>>::iterator endDemand() {
    return static_cast<CumulativeExpressionImpl<T> *>(
               BooleanVar<T>::implem().get())
        ->endDemand();
  }

  std::vector<BooleanVar<T>>::const_iterator endDemand() const {
    const_cast<CumulativeExpression *>(this)->endDemand();
  }

  void provide(const NumericVar<T> d, const Interval<T> i) {
    static_cast<CumulativeExpressionImpl<T> *>(BooleanVar<T>::implem().get())
        ->provide(d, i);
  }
};

template <typename T = int>
CumulativeExpression<T> Cumulative(const Interval<T> &s, const NumericVar<T> c,
                                   const std::vector<Interval<T>> &I,
                                   const std::vector<NumericVar<T>> &D) {
  auto impl = std::make_shared<CumulativeExpressionImpl<T>>(s, c);
  auto i{I.begin()};
  auto d{D.begin()};
  while (i != I.end()) {
    impl->provide(*d, *i);
    ++i;
    ++d;
  }

  CumulativeExpression<T> exp(std::move(impl));
  return exp;
}

namespace detail {
template <typename T>
concept literal_matrix =
    concepts::same_template<T, Matrix> and
    concepts::same_template<typename std::remove_cvref_t<T>::value_type,
                            Literal>;

template <typename T>
class single_forward_view
    : public std::ranges::view_interface<single_forward_view<T>> {
  T data;

public:
  template <typename U>
  explicit constexpr single_forward_view(U &&data)
      : data(std::forward<U>(data)) {}

  constexpr auto begin() const noexcept { return &data; }

  constexpr auto begin() noexcept { return &data; }

  constexpr auto end() const noexcept { return &data + 1; }

  constexpr auto end() noexcept { return &data + 1; }
};

template <typename T> constexpr auto single_forward(T &&data) {
  return single_forward_view<T>(std::forward<T>(data));
}
} // namespace detail

/**
 * @brief concept for different resource expressions (NoOverlapExpression,
 * CumulativeExpression, ...)
 * @tparam E expression type
 */
template <typename E>
concept resource_expression =
    tempo::concepts::ttyped_range<E, tempo::Interval> and
    requires(E expression) {
      { expression.getDisjunctiveLiterals() } -> detail::literal_matrix;
    };

namespace detail {
template <resource_expression R> class timing_type_from_resource {
  using LitT = typename std::remove_cvref_t<
      decltype(std::declval<R>().getDisjunctiveLiterals())>::value_type;

public:
  using type = decltype(std::declval<LitT>().value());
};

template <resource_expression R>
using timing_type_from_resource_t = typename timing_type_from_resource<R>::type;
} // namespace detail

/**
 * @brief concept modelling a range of resource expressions
 * @tparam T range type
 */
template <typename T>
concept resource_range = std::ranges::range<T> and
                         resource_expression<std::ranges::range_value_t<T>>;

/**
 * Extracts the unique boolean varaibles from a range of resource expressions
 * @tparam Resources resource expression range type
 * @param resources range of resources
 * @return vector of unique boolean variables
 */
template <resource_range Resources>
auto booleanVarsFromResources(Resources &&resources) {
  using T = detail::timing_type_from_resource_t<
      std::ranges::range_value_t<Resources>>;
  using namespace std::views;
  std::vector<BooleanVar<T>> variables;
  for (const auto &resource : std::forward<Resources>(resources)) {
    auto resVars =
        resource.getDisjunctiveLiterals().rawData() |
        filter([](auto lit) { return lit != Contradiction<T>; }) |
        transform([](auto lit) { return tempo::BooleanVar<T>(lit); });
    std::ranges::copy(resVars, std::back_inserter(variables));
  }

  std::ranges::sort(variables, {}, [](const auto &var) { return var.id(); });
    auto res = std::ranges::unique(variables, [](const auto &var1, const auto &var2) {return var1.id() == var2.id();});
  variables.erase(res.begin(), res.end());
  variables.shrink_to_fit();
  return variables;
}

/**
 * overload of booleanVarsFromResources accepting a single resource
 * @tparam R resource type
 * @param resource resource
 * @return vector of unique boolean variables
 */
template <resource_expression R> auto booleanVarsFromResources(R &&resource) {
  auto range = detail::single_forward(std::forward<R>(resource));
  return booleanVarsFromResources(range);
}

/*!
NumericVar  impl
*/
template <typename T>
template <distance_provider S>
T NumericVar<T>::min(const S &s) const {
  T v = s.numeric.lower(id());
  if (v == -Constant::Infinity<T>)
    return v;
  return v + offset();
}

template <typename T>
template <distance_provider S>
T NumericVar<T>::max(const S &s) const {
  T v = s.numeric.upper(id());
  if (v == Constant::Infinity<T>)
    return v;
  return v + offset();
}

template <typename T>
template <distance_provider S>
T NumericVar<T>::earliest(const S &s) const {
  return min(s);
}

template <typename T>
template <distance_provider S>
T NumericVar<T>::latest(const S &s) const {
  return max(s);
}

template <typename T> Literal<T> NumericVar<T>::after(const T t) const {
  return geq<T>(id(), (t == Constant::Infinity<T> ? t : t - offset()));
}

template <typename T> Literal<T> NumericVar<T>::before(const T t) const {
  return leq<T>(id(), (t == Constant::Infinity<T> ? t : t - offset()));
}

template <typename T>
Literal<T> NumericVar<T>::greaterThanOrEqual(const T t) const {
  return after(t);
}

template <typename T>
Literal<T> NumericVar<T>::lessThanOrEqual(const T t) const {
  return before(t);
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

template <typename T>
DistanceConstraint<T> NumericVar<T>::lessThanOrEqual(const NumericVar<T> &e,
                                                     const T t) const {
  return before(e, t);
}

template <typename T>
DistanceConstraint<T> NumericVar<T>::greaterThanOrEqual(const NumericVar<T> &e,
                                                        const T t) const {
  return e.before(*this, t);
}

/*!
BooleanVar  impl
*/
template <typename T>
std::ostream &BooleanVar<T>::display(std::ostream &os) const {
  if (this->isExpression()) {
//               os << "expr ";
    os << "(";
    this->implem()->display(os);
    os << ")";
  } else {
    os << "b" << id();
  }
  return os;
}

template <typename T>
std::ostream &NumericVar<T>::display(std::ostream &os) const {
  if (this->isExpression()) {
    //            os << *(this->implem());
    this->implem()->display(os);
  } else if (id() == Constant::NoVar) {
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

template <typename T>
std::ostream &operator<<(std::ostream &os, const BooleanVar<T> &x) {
  return x.display(os);
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const NumericVar<T> &x) {
  return x.display(os);
}

/*!
Interval  impl
*/

template <typename T>
//    template<distance_provider S>
Interval<T>::Interval(/*S &solver,*/ const NumericVar<T> s,
                      const NumericVar<T> e, const NumericVar<T> d,
                      const BooleanVar<T> o)
    : _id_(num_intervals++), start(s), end(e), duration(d), exist(o) {

  //            std::cout << "this interval\n";
  //    extract(solver);
}

template<typename T>
void Interval<T>::extract(Solver<T> &solver) {

  start.extract(solver);
  //        std::cout << "\ncreate start " << start << "\n";

  end.extract(solver);
  //        std::cout << "\ncreate end " << end << "\n";

  duration.extract(solver);
  //        std::cout << "\ncreate duration " << duration << "\n";

  exist.extract(solver);
  //        std::cout << "\ncreate optionality " << exist << "\n";

  //
  //  solver.post((start + duration) == end);
  if (start.id() != end.id()) {
    solver.post(start.before(end, duration.min(solver)));
    //
    if (duration.max(solver) != Constant::Infinity<T>)
      solver.post(end.before(start, -duration.max(solver)));
  }
}

//    template<typename T>
//    template<distance_provider S>
//    Interval<T>::Interval(S &solver, const T mindur, const T maxdur,
//                          const T earliest_start, const T latest_start,
//                          const T earliest_end, const T latest_end,
//                          const BooleanVar<T> opt)
//        : _id_(num_intervals++), exist(opt) {
//  
//        if (earliest_start == latest_start) {
//            start = NumericVar(Constant::K, earliest_start);
//        } else {
//            start = solver.newNumeric(earliest_start, latest_start);
//        }
//
//        if (mindur == maxdur) {
//            end = NumericVar(start.id(), mindur);
//            duration = NumericVar(Constant::K, mindur);
//        } else {
//            auto s{start.min(solver)};
//            if (s != start.max(solver)) {
//                end = solver.newNumeric(earliest_end, latest_end);
//                duration = solver.newNumeric(mindur, maxdur);
//                solver.post((start + duration) == end);
//                solver.post(start.before(end, mindur));
//
//                if (maxdur != Constant::Infinity<T>)
//                    solver.post(end.before(start, -maxdur));
//            } else {
//                auto ect{std::max(earliest_end, s + mindur)};
//                auto lct{std::min(latest_end, s + maxdur)};
//                end = solver.newNumeric(ect, lct);
//                duration = NumericVar(end.id(), -s);
//            }
//        }
//            
//            
////            std::cout << *this << std::endl;
//    }

template <typename T> int Interval<T>::id() const {
  return /*start.id()*/ _id_;
}

// @TODO: remove that, should be a constraint
template <typename T> bool Interval<T>::operator==(const Interval<T> &t) const {
  return id() == t.id();
}

template <typename T>
template <distance_provider S>
T Interval<T>::getEarliestStart(const S &solver) const {
  return start.earliest(solver);
}

template <typename T>
template <distance_provider S>
T Interval<T>::getLatestStart(const S &solver) const {
  return start.latest(solver);
}

template <typename T>
template <distance_provider S>
T Interval<T>::getEarliestEnd(const S &solver) const {
  return end.earliest(solver);
}

template <typename T>
template <distance_provider S>
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
template <distance_provider S>
T Interval<T>::minDuration(const S &solver) const {
  return duration.min(solver);
}

template <typename T>
template <distance_provider S>
T Interval<T>::maxDuration(const S &solver) const {
  return duration.max(solver);
}

template <typename T> var_t Interval<T>::getStart() const { return start.id(); }

template <typename T> var_t Interval<T>::getEnd() const { return end.id(); }

template <typename T> void Interval<T>::require(NoOverlapExpression<T> &R) {
  R.provide(*this);
}

// template <typename T>
// void Interval<T>::require(const T d, CumulativeExpression<T>& R) {
//     R.provide(d, *this);
// }

template <typename T>
void Interval<T>::require(const NumericVar<T> d, CumulativeExpression<T> &R) {
  R.provide(d, *this);
}

template <typename T>
std::ostream &Interval<T>::display(std::ostream &os) const {
  //        os << "t" << id(); //<< ": [" << start.earliest(solver) << ".."
  //        <<
  //        // end.latest(solver) << "]";
  //
  //        //    auto est{start.earliest(solver)};
  //        //    auto lst{start.latest(solver)};
  //        //    auto ect{end.earliest(solver)};
  //        //    auto lct{end.latest(solver)};
  //        //    auto pmin{duration.min(solver)};
  //        //    auto pmax{duration.max(solver)};
  //        //
  //        //   os << ": [" << est << "-" << lst << ".." << ect << "-" <<
  //        lct << "] (" <<
  //        //   pmin << "-" << pmax << ")";
  //
  //        os << ": from " << start << " to " << end << " for "
  //                << duration; //<< "/" << exist;
  //
  //        os << " if " << exist;

  os << "[" << start << ", " << end << "]";
  if (exist.id() != 0) {
    os << ":" << exist;
  }
//    os << " dur=" << duration;
  return os;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const Interval<T> &x) {
  return x.display(os);
}

template<typename T>
class StaticBooleanStore;

template<typename T>
class StaticNumericStore;

template<typename T>
class Model {

public:
  Model() {
    newNumeric();
    numeric.set(geq<T>(0, 0));
    numeric.set(leq<T>(0, 0));
  }

  // create an internal Boolean variable and return a model object pointing to
  // it
  BooleanVar<T> newBoolean();
  // create an internal Boolean variable with a difference logic semantic,
  // post the channelling constraints, and return a model object pointing to
  // it
  BooleanVar<T> newDisjunct(const DistanceConstraint<T> &,
                            const DistanceConstraint<T> &);
  // returns the constant 0
  NumericVar<T> zero() { return NumericVar<T>(Constant::K, 0); }
  // returns the constant true
  BooleanVar<T> truism() { return BooleanVar<T>(0); }
  // create an internal numeric variable and return a model object pointing to
  // it
  NumericVar<T> newNumeric(const T lb = -Constant::Infinity<T>,
                           const T ub = Constant::Infinity<T>);
  // create a modeling object representing a constant
  NumericVar<T> newConstant(const T k);
  NumericVar<T> newOffset(NumericVar<T> &x, const T k);
  // create an internal temporal variable and return a model object pointing
  // to it
  //    TemporalVar<T> newTemporal(const T offset = 0);
  //    NumericVar<T> newTemporal(const T offset = 0);
  // create the internal variables (depending on the type of Interval) and
  // return a model object pointing to them
  //    Interval<T> newInterval(const T mindur = 0,
  //                            const T maxdur = Constant::Infinity<T>,
  //                            const T earliest_start =
  //                            -Constant::Infinity<T>, const T latest_start =
  //                            Constant::Infinity<T>, const T earliest_end =
  //                            -Constant::Infinity<T>, const T latest_end =
  //                            Constant::Infinity<T>, const BooleanVar<T> opt
  //                            = Constant::True
  //                            );

  Interval<T> maybe_between(const NumericVar<T> s, const NumericVar<T> e);
  Interval<T> maybe_continuefor(const NumericVar<T> s, const NumericVar<T> d);

  Interval<T> between(const NumericVar<T> s, const NumericVar<T> e);
  Interval<T> continuefor(const NumericVar<T> s, const NumericVar<T> d);

  Interval<T> between(const NumericVar<T> s, const NumericVar<T> e,
                      const BooleanVar<T> opt);
  Interval<T> continuefor(const NumericVar<T> s, const NumericVar<T> d,
                          const BooleanVar<T> opt);

  void post(BooleanVar<T> expr) { constraint.push_back(expr); }

  void post(DistanceConstraint<T> c) { precedence.push_back(c); }

  void post(Literal<T> l) { literal.push_back(l); }

  std::ostream &display(std::ostream &os) const;

  std::vector<BooleanVar<T>> boolean_var;
  std::vector<NumericVar<T>> numeric_var;
  std::vector<Interval<T>> interval;

  std::vector<DistanceConstraint<T>> precedence;
  std::vector<Literal<T>> literal;
  //        std::vector<std::vector<Literal<T>>> clause;

  std::vector<BooleanVar<T>> constraint;
  std::vector<NoOverlapExpression<T>> disjunctive_resource;
  std::vector<CumulativeExpression<T>> cumulative_resource;

  StaticBooleanStore<T> boolean;
  StaticNumericStore<T> numeric;
};


template <typename T> BooleanVar<T> Model<T>::newBoolean() {
  boolean_var.push_back(boolean.newVar());
  return boolean_var.back();
}

template <typename T>
BooleanVar<T> Model<T>::newDisjunct(const DistanceConstraint<T> &d1,
                                    const DistanceConstraint<T> &d2) {
  boolean_var.push_back(boolean.newVnewDisjunctar(d1, d2));
  return boolean_var.back();
}

template <typename T> NumericVar<T> Model<T>::newConstant(const T k) {
  return NumericVar(Constant::K, k);
}

template <typename T>
NumericVar<T> Model<T>::newOffset(NumericVar<T> &x, const T k) {
  return NumericVar<T>(x.id(), k);
}

template <typename T>
NumericVar<T> Model<T>::newNumeric(const T lb, const T ub) {
  if (lb > ub) {
    throw Failure<T>(Constant::NoReason<T>);
  } else if (lb == ub) {
    return NumericVar(Constant::K, lb);
  } else {
    auto x{numeric.newVar(Constant::Infinity<T>)};
    numeric_var.push_back(x);
    //        constraint.push_back(x >= lb);
    //        constraint.push_back(x <= ub);
    numeric.set(geq<T>(x.id(), lb));
    numeric.set(leq<T>(x.id(), ub));
    return numeric_var.back();
  }
}

// template <typename T>
// Interval<T> Model<T>::newInterval(const T mindur, const T maxdur,
//                                    const T earliest_start, const T
//                                    latest_start, const T earliest_end, const
//                                    T latest_end, const BooleanVar<T> opt) {
//     interval.emplace_back(*this, mindur, maxdur, earliest_start,
//     latest_start, earliest_end, latest_end, opt); return interval.back();
// }

template <typename T>
Interval<T> Model<T>::between(const NumericVar<T> s, const NumericVar<T> e) {
  interval.emplace_back(s, e, e - s, BooleanVar<T>(Constant::True));
  //    numeric.set(interval.back().duration >= 0);
  //    numeric.set(geq<T>(interval.back().duration.id(), 0));
  post(interval.back().duration >= 0);
  return interval.back();
}

template <typename T>
Interval<T> Model<T>::continuefor(const NumericVar<T> s, const NumericVar<T> d) {
  interval.emplace_back(s, s + d, d, BooleanVar<T>(Constant::True));
  return interval.back();
}

template <typename T>
Interval<T> Model<T>::maybe_between(const NumericVar<T> s, const NumericVar<T> e) {
  interval.emplace_back(s, e, e - s, newBoolean());
  post(interval.back().duration >= 0);
  //    numeric.set(geq<T>(interval.back().duration.id(), 0));
  return interval.back();
}

template <typename T>
Interval<T> Model<T>::maybe_continuefor(const NumericVar<T> s, const NumericVar<T> d) {
  interval.emplace_back(s, s + d, d, newBoolean());
  return interval.back();
}

template <typename T>
Interval<T> Model<T>::between(const NumericVar<T> s, const NumericVar<T> e, const BooleanVar<T> optional) {
  interval.emplace_back(s, e, e - s, optional);
  post(interval.back().duration >= 0);
  //    numeric.set(geq<T>(interval.back().duration.id(), 0));
  return interval.back();
}

template <typename T>
Interval<T> Model<T>::continuefor(const NumericVar<T> s, const NumericVar<T> d, const BooleanVar<T> optional) {
  interval.emplace_back(s, s + d, d, optional);
  return interval.back();
}


template <typename T>
std::ostream &Model<T>::display(std::ostream &os) const {
  if (not boolean_var.empty())
    os << "b1,...,b" << boolean_var.back().id() << std::endl;
    os << "numeric variables:\n";
  for (auto x : numeric_var) {
    os << "x" << x.id() << " in [" << x.min(*this) << ".." << x.max(*this)
       << "]\n";
  }
    os << "intervals:\n";
    for (auto i : interval) {
        os << i << std::endl;
    }
    os << "expressions:\n";
    for (auto c : constraint) {
      os << c << std::endl;
    }
//    os << "precedences:\n";
//    for (auto p : precedence) {
//      os << p << std::endl;
//    }
return os;
}



////template<typename T>
////class SchedulingModel;
//
//
//
//template<typename T>
//class SchedulingInstance {
//    
//public:
//    
//    struct ModNum {
//        int id;
//        T offset;
//        
//        bool isConstant() {return id == 0;}
//        
//        std::ostream &display(std::ostream &os) const {
//            if(id > 0) {
//                os << "x" << id;
//                if(offset > 0) {
//                    os << "+" << offset;
//                } else if(offset < 0) {
//                    os << offset;
//                }
//            } else {
//                os << offset;
//            }
//            return os;
//        }
//    };
//    
//private:
//    
//    struct ModPrecedence {
//        ModNum start{0};
//        ModNum end{0};
//        T lag{0};
//    };
//    
//    struct ModInterval {
//        ModNum start{0};
//        ModNum end{0};
//        ModNum duration{0};
//        int optional{-2};
//    };
//    
//    struct CardConstraint {
//        std::vector<int> vars;
//        T lb;
//        T ub;
//    };
//    
//public:
//    SchedulingInstance() {
//        LBs.push_back(0);
//        UBs.push_back(0);
//        variables.push_back()
//    }
//    
////    size_t numConstant() { return constant.size(); }
//    size_t numNumeric() { return LBs.size(); }
//    size_t numBoolean() { return static_cast<size_t>(num_booleans); }
//    
//    ModNum getStart(const int i) const {
//        return interval_rep[i].start;
//    }
//    
//    ModNum getEnd(const int i) const {
//        return interval_rep[i].end;
//    }
//    
//    ModNum getDuration(const int i) const {
//        return interval_rep[i].duration;
//    }
//    
//    int getOptional(const int i) const {
//        return interval_rep[i].optional;
//    }
//    
//    T getLb(const ModNum x) {
//        return LBs[x.id] + x.offset;
//    }
//    
//    T getUb(const ModNum x) {
//        return UBs[x.id] + x.offset;
//    }
//    
////    void encode(Solver<T>& solver, Interval<T>& schedule, std::vector<BooleanVar<T>>& boolean, std::vector<NumericVar<T>>& numeric, std::vector<Interval<T>>& interval, std::vector<DistanceConstraint<T>>& constraint);
////    void encode(Solver<T>& solver, SchedulingModel<T>& model);
//    
//    int addBoolean() {
//        return num_booleans++;
//        //        return (val == ? num_booleans++ : val);
//    }
//    ModNum addNumeric(const T lb=-Constant::Infinity<T>, const T ub=Constant::Infinity<T>) {
//        if(lb != ub) {
//            auto x{static_cast<int>(LBs.size())};
//            LBs.push_back(lb);
//            UBs.push_back(ub);
//            variables.emplace_back(x,0);
////            return ModNum(x,0);
//            variables.back();
//        } else {
////            variables.emplace_back(x,0);
//            return ModNum(0,lb);
//        }
//    }
//    ModNum addView(const ModNum x, const T offset) {
//        variables.emplace_back(x.id,x.offset + offset);
//        return variables.back();
////        return ModNum(x.id,x.offset + offset);
//    }
//    ModNum addConstant(const T v) {
//        return ModNum(0,v);
//    }
//    void declareCardinalityConstraints(const int n) {
//        cardinalities.resize(n);
//    }
//    void addCardinalityArgument(const int i, const int c) {
//        cardinalities[c].vars.push_back(i);
//    }
//    void setCardinalityBounds(const int c, const T l, const T u) {
//        cardinalities[c].lb = l;
//        cardinalities[c].ub = u;
//    }
//    void declareDisjunctiveResources(const int n) {
//        disjunctive_resources.resize(disjunctive_resources.size() + n);
//        disjunctive_resource_transitions.resize(disjunctive_resources.size());
//    }
//    void declareCumulativeResources(const int n) {
//        cumulative_resources.resize(cumulative_resources.size() + n);
//    }
//    void addDisjunctiveResourceUsage(const int i, const int r) {
//        disjunctive_resources[r].push_back(i);
//    }
//    void addCumulativeResourceUsage(const int i, const int r, const int x) {
//        cumulative_resources[r].emplace_back(i,x);
//    }
//    int addPrecedence(const ModNum t1, const ModNum t2, const T lag=0) {
//        auto p{static_cast<int>(precedence_rep.size())};
//        precedence_rep.emplace_back(t1, t2, lag);
//        return p;
//    }
//    int addInterval(const ModNum t1, const ModNum t2, const ModNum p, const int b=-1) {
////        auto l{addNumeric(min_duration, max_duration)};
//        auto i{static_cast<int>(interval_rep.size())};
//        interval_rep.emplace_back(t1, t2, p, b);
//        return i;
//    }
////    int addFixedDurationIntervalFrom(const int t, const T duration) {
////        auto l{addNumeric(duration, duration)};
////        auto i{static_cast<int>(interval_rep.size())};
////        interval_rep.emplace_back(t, -1, l, -1);
////        return i;
////    }
////    int addFixedDurationOptionalIntervalFrom(const int t, const T duration) {
////        auto l{addNumeric(duration, duration)};
////        auto b{addBoolean()};
////        auto i{static_cast<int>(interval_rep.size())};
////        interval_rep.emplace_back(t, -1, l, b);
////        return i;
////    }
//    
//    void addSameAllocation(const int i, const int j) {
////        equivalences.resize(interval_rep.size());
//        equivalences.emplace_back(i,j);
//    }
//    
//    //    int addInterval(const int t1, const int t2, const T min_duration=0, const T max_duration=Constant::Infinity<T>, const int b=0) {
//    //        auto l{addNumeric(min_duration, max_duration)};
//    //        interval_rep.emplace_back(t1, t2, l, b);
//    //    }
//    
////    int origin; // var id for origin
//    ModNum makespan; // var id for makespan
//    
//    std::vector<ModNum> variables;
//    std::vector<ModInterval> interval_rep;
//    std::vector<ModPrecedence> precedence_rep;
//    std::vector<std::vector<int>> disjunctive_resources;
//    std::vector<std::vector<std::pair<int,int>>> cumulative_resources;
//    std::vector<std::vector<std::vector<T>>> disjunctive_resource_transitions;
//    std::vector<std::pair<int,int>> equivalences;
//    std::vector<CardConstraint> cardinalities;
//    
//private:
//    
//    // solver agnostic model
//    int num_booleans{0};
//    int num_timepoints{0};
//    std::vector<T> LBs;
//    std::vector<T> UBs;
////    std::vector<T> constant;
//    
//    
//};
//  
//template<typename T>
//class SchedulingModel {
//    
//public:
//    SchedulingModel(SchedulingInstance<T>& data, Solver<T>& solver);
////    {
////        data.encode(solver, schedule, boolean, numeric, interval, constraint);
////    }
//    
//    Interval<T> getScheduleInterval() { return schedule; }
//    
//private:
//    
//    // model in
//    Interval<T> schedule;
//    
//    std::vector<BooleanVar<T>> boolean;
//    std::vector<NumericVar<T>> numeric;
//    std::vector<Interval<T>> interval;
//    
//    std::vector<DistanceConstraint<T>> precedence;
//    std::vector<std::vector<Literal<T>>> clause;
//    
//    std::vector<BooleanVar<T>> constraint;
//    
//    std::vector<NoOverlapExpression<T>> disjunctive_resource;
//    std::vector<CumulativeExpression<T>> cumulative_resource;
//};
//
//template<typename T>
////void SchedulingInstance<T>::encode(Solver<T>& solver, Interval<T>& schedule, std::vector<BooleanVar<T>>& boolean, std::vector<NumericVar<T>>& numeric, std::vector<Interval<T>>& interval, std::vector<DistanceConstraint<T>>& precedence) {
////void SchedulingInstance<T>::encode(Solver<T>& solver, SchedulingModel<T>& model) {
//SchedulingModel<T>::SchedulingModel(SchedulingInstance<T>& data, Solver<T>& solver) {
////
////    model.schedule = solver.newInterval(0, Constant::Infinity<int>, 0, 0, 0,
////                  Constant::Infinity<int>);
////    
////    
////    for(size_t i{0}; i<numBoolean(); ++i) {
////        model.boolean.push_back(solver.newBoolean());
////    }
////    
////    for(size_t i{0}; i<numNumeric(); ++i) {
////        model.numeric.push_back(solver.newNumeric(LBs[i], UBs[i]));
////    }
////    
////    
////    for(auto I : interval_rep) {
////        if(I.optional >= 0) {
////            if(LBs[I.duration] == UBs[I.duration]) {
////                model.interval.push_back(solver.between(model.numeric[I.start], model.numeric[I.start]+LBs[I.duration], model.boolean[I.optional]));
////            } else {
////                model.interval.push_back(solver.between(model.numeric[I.start], model.numeric[I.end], model.boolean[I.optional]));
////            }
////        } else {
////            if(LBs[I.duration] == UBs[I.duration]) {
////                model.interval.push_back(solver.between(model.numeric[I.start], model.numeric[I.start]+LBs[I.duration]));
////            } else {
////                model.interval.push_back(solver.between(model.numeric[I.start], model.numeric[I.end]));
////            }
////        }
////        
////        solver.set(model.interval.back().end.before(model.schedule.end));
////    }
////    
////    
////    for(auto P : precedence_rep) {
////        auto prec{model.numeric[P.end].after(model.numeric[P.start], P.lag)};
////        model.precedence.push_back(prec);
////
//////                std::cout << " add " << constraint.back() << " (" << P.lag << ")"<< std::endl;
////        
////        solver.set(prec);
////    }
////    
////    
////    std::vector<Interval<T>> scope;
////    int i{0};
////    for (auto &tasks : disjunctive_resources) {
////        for (auto j : tasks) {
////            scope.push_back(model.interval[j]);
////        }
////        solver.disjunctive_resource.push_back(NoOverlap(solver.schedule, scope, disjunctive_resource_transitions[i++]));
////        solver.post(solver.disjunctive_resource.back());
////        scope.clear();
////    }
////    
////    for(auto p : equivalences) {
////        solver.addClause({model.boolean[p.first] == false, model.boolean[p.second] == true});
////        solver.addClause({model.boolean[p.first] == true, model.boolean[p.second] == false});
////    }
////    
//////    for(auto )
////    
////    
//////    solver.set({0, 1, solver.numeric.upper(1)});
////    
////
//    
//    
////    schedule = solver.newInterval(0, Constant::Infinity<int>, 0, 0, 0,
////                  Constant::Infinity<int>);
//    
//    
//    for(size_t i{0}; i<data.numBoolean(); ++i) {
//        boolean.push_back(solver.newBoolean());
//    }
//    
////    for(size_t i{0}; i<data.numNumeric(); ++i) {
////        numeric.push_back(solver.newNumeric(data.getLb(i), data.getUb(i)));
////    }
//    for( auto x : data.variables) {
//        if(x.)
//        numeric.push_back(solver.newNumeric(data.getLb(i), data.getUb(i)));
//    }
//    
//    
//    schedule = solver.continuefor(data.origin, data.makespan);
//    
//    
//    std::cout << "add variables\n";
//    std::cout << solver << std::endl;
//    
//    for(auto x : numeric) {
//        std::cout << x << std::endl;
//    }
//    
//    
//    for(auto& I : data.interval_rep) {
//        if(I.optional >= 0) {
//            if(data.getLb(I.duration) == data.getUb(I.duration)) {
//                interval.push_back(solver.between(numeric[I.start], numeric[I.start]+data.getLb(I.duration), boolean[I.optional]));
//            } else {
//                interval.push_back(solver.between(numeric[I.start], numeric[I.end], boolean[I.optional]));
//            }
//        } else {
//            if(data.getLb(I.duration) == data.getUb(I.duration)) {
//                interval.push_back(solver.between(numeric[I.start], numeric[I.start]+data.getLb(I.duration)));
//            } else {
//                interval.push_back(solver.between(numeric[I.start], numeric[I.end]));
//            }
//        }
//        
////        solver.set(interval.back().end.before(schedule.end));
//    }
//    
//    
//    std::cout << "add intervals\n";
//    std::cout << solver << std::endl;
//    
//    
//    for(auto P : data.precedence_rep) {
//        auto prec{numeric[P.end].after(numeric[P.start], P.lag)};
//        precedence.push_back(prec);
//
//                std::cout << " add " << precedence.back() << " (" << P.lag << ")"<< std::endl;
//        
//        solver.set(prec);
//    }
//    
//    
//    std::vector<Interval<T>> scope;
//    int i{0};
//    for (auto &tasks : data.disjunctive_resources) {
//        for (auto j : tasks) {
//            scope.push_back(interval[j]);
//        }
//        disjunctive_resource.push_back(NoOverlap(schedule, scope, data.disjunctive_resource_transitions[i++]));
//        solver.post(disjunctive_resource.back());
//        scope.clear();
//    }
//    
////    for(auto p : data.equivalences) {
////        solver.clauses.add({boolean[p.first] == false, boolean[p.second] == true});
////        solver.clauses.add({boolean[p.first] == true, boolean[p.second] == false});
////    }
// 
//}
//
//
//template<typename T>
//std::ostream &operator<<(std::ostream &os, const typename SchedulingInstance<T>::ModNum &x) {
//    return x.display(os);
//}

template<typename T>
std::ostream &operator<<(std::ostream &os, const Model<T> &x) {
  return x.display(os);
}
} // namespace tempo

#endif // __MODEL_HPP
