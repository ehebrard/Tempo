
#ifndef __TEMPO_LITERAL_HPP
#define __TEMPO_LITERAL_HPP

#include <stdint.h>
#include <variant>

#include "Constant.hpp"

namespace tempo {


/**********************************************
* Literal
**********************************************/

enum bound { lower = 0, upper = 1 };



template<typename T>
struct Literal {
    Literal() = default;
//    Literal(const var_t i, const T d);
    Literal(const var_t i, const std::variant<info_t,T> d);
    Literal(const bool sign, const var_t x, const T v);
    Literal(const bool sign, const var_t x, const info_t v);

    bool isBoolean() const;
    bool hasSemantic() const;
    bool isNumeric() const;

    //    operator int() const;
    operator info_t() const;

    bool sign() const;
    var_t variable() const;

    bool sameVariable(Literal<T> l) const;

    T value() const;
    info_t constraint() const;
    
    std::ostream& display(std::ostream &os) const;
    
    static var_t index(const bool sign, const var_t x) ;
    
    var_t _id_{0};
    std::variant<info_t,T> _data_{0};
};

template <typename T> Literal<T>::Literal(const var_t i, const std::variant<info_t,T> d) : _id_(i), _data_(d) {}

//template <typename T> Literal<T>::Literal(const var_t i, const T d) : _id_(i), _data_(d) {}

template <typename T>
Literal<T>::Literal(const bool sign, const var_t x, const T v)
    : _id_(index(sign,x)), _data_(v) {
}

template <typename T>
Literal<T>::Literal(const bool sign, const var_t x, const info_t v)
    : _id_(index(sign,x)), _data_(v) {
}

template <typename T> bool Literal<T>::sameVariable(Literal<T> l) const {
  return isNumeric() == l.isNumeric() and variable() == l.variable();
}

template <typename T> var_t Literal<T>::index(const bool sign, const var_t x) {
    return 2 * x + sign;
}

template <typename T> Literal<T>::operator info_t() const { return _id_; }
// template <typename T> Literal<T>::operator int() const { return
// static_cast<int>(_id_); }

template <typename T> bool Literal<T>::sign() const { return (_id_ & 1); }

template<typename T>
var_t Literal<T>::variable() const {
    return _id_ / 2;
}

template<typename T>
T Literal<T>::value() const {
    return std::get<T>(_data_);
}

template<typename T>
info_t Literal<T>::constraint() const {
    return std::get<info_t>(_data_) + sign();
}

template<typename T>
bool Literal<T>::isBoolean() const {
    return std::holds_alternative<info_t>(_data_);
}

template<typename T>
bool Literal<T>::isNumeric() const {
    return not isBoolean();
}

template<typename T>
bool Literal<T>::hasSemantic() const {
    return isNumeric() or std::get<info_t>(_data_) != Constant::NoSemantic;
}

template <typename T> Literal<T> operator~(const Literal<T> l) {
  if (l.isNumeric()) {
    return Literal<T>(not l.sign(), l.variable(),
                      -std::get<T>(l._data_) - Gap<T>::epsilon());
  } else {
//    auto d{std::get<info_t>(l._data_)};
//    if (d != Constant::NoSemantic)
//      return Literal<T>(l._id_ ^ 1, d ^ 1);
    return Literal<T>(l._id_ ^ 1, std::get<info_t>(l._data_));
  }
}

template <typename T> std::ostream& Literal<T>::display(std::ostream &os) const {
    if(_id_ == Constant::NoVarx) {
        os << "constant: " << value();
    } else if(isNumeric()) {
        os << (sign() ? "x" : "-x") << variable() << " <= " << value();
    } else {
      os << (sign() ? "b" : "¬b") << variable() << (hasSemantic() ? "*" : "");
    }
    return os;
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const Literal<T> &x) {
  return x.display(os);
}

template <typename T> Literal<T> ub(const var_t x) {
  return Literal<T>(true, x, Constant::Infinity<T>);
}

template <typename T> Literal<T> lb(const var_t x) {
  return Literal<T>(false, x, Constant::Infinity<T>);
}

template<typename T>
Literal<T> leq(const var_t x, const T v) {
    return Literal<T>(true, x, v);
}

template<typename T>
Literal<T> gt(const var_t x, const T v) {
    return Literal<T>(false, x, -v-Gap<T>::epsilon());
}

template<typename T>
Literal<T> lt(const var_t x, const T v) {
    return Literal<T>(true, x, v-Gap<T>::epsilon());
}

template<typename T>
Literal<T> geq(const var_t x, const T v) {
    return Literal<T>(false, x, -v);
}

template<typename T>
Literal<T> pos(const var_t x, const info_t d=0) {
    return Literal<T>(true, x, d);
}

template<typename T>
Literal<T> neg(const var_t x, const info_t d=0) {
    return Literal<T>(false, x, d);
}


}

#endif // __Literal_HPP
