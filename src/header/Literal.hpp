
#ifndef __TEMPO_LITERAL_HPP
#define __TEMPO_LITERAL_HPP

#include <variant>

#include "Constant.hpp"
#include "Global.hpp"

namespace tempo {


/**********************************************
* Literal
**********************************************/

enum bound { lower = 0, upper = 1 };

namespace detail {
    /**
     * Tag class to resolve ctor ambiguity
     * @tparam value held
     */
    template<typename T>
    struct NumericValue {
        explicit constexpr NumericValue(T val) noexcept: val(val) {}

        T val;
    };

    /**
     * Tagged union class that holds a variable ID and either numeric data or constraint info
     * @tparam T data ype
     */
    template<typename T>
    class LiteralStorage {
        static constexpr std::uint32_t SignBit = 2;
        static constexpr std::uint32_t TagBit = 1;
        std::uint32_t info;
        union {
            T numericData;
            info_t semanticInfo;

        };

    public:
        /**
         * Ctor
         * constructs a constraint literal
         * @param sign sign of the literal
         * @param literalId identifier of the literal (not of the variable!)
         * @param semanticInfo constraint information
         * @note the id of the corresponding variable should be litId / 2. See also makeSemanticLit
         */
        constexpr LiteralStorage(info_t literalId, info_t semanticInfo) noexcept:
                info(literalId << 1), semanticInfo(semanticInfo) {}

        /**
         * Ctor
         * constructs a numeric literal
         * @param sign sign of the literal
         * @param literalId identifier of the literal (not of the variable!)
         * @param numericVal numeric value stored in the literal
         * @note see also makeNumericLit
         */
        constexpr LiteralStorage(info_t literalId, NumericValue<T> numericVal) noexcept: info(
                (literalId << 1) | TagBit), numericData(std::move(numericVal.val)) {}

        /**
         * Checks whether the storage holds a numeric value
         * @return true if storage holds a numeric value, false otherwise
         */
        [[nodiscard]] constexpr bool isNumeric() const noexcept {
            return info & TagBit;
        }

        /**
         * Gets the sign of the literal
         * @return returns the sign of the literal
         */
        [[nodiscard]] constexpr bool sign() const noexcept {
            return info & SignBit;
        }

        /**
         * Gets the id of the literal
         * @return the literal's id
         */
        [[nodiscard]] constexpr std::uint32_t id() const noexcept {
            return (info & ~TagBit) >> 1;
        }

        /**
         * Access to the stored value
         * @return the stored value
         * @throws std::bad_variant_access if storage does not contain numeric value
         */
        [[nodiscard]] constexpr T value() const {
            if (not isNumeric()) {
                throw std::bad_variant_access();
            }

            return numericData;
        }

        /**
         * Access to the stored constraint information
         * @return the stored constraint information
         * @throws std::bad_variant_access if storage does not contain constraint information
         */
        [[nodiscard]] constexpr info_t constraint() const {
            if (isNumeric()) {
                throw std::bad_variant_access();
            }

            return semanticInfo;
        }

        void setValue(T value) noexcept {
            numericData = value;
            info |= TagBit;
        }
    };

    /**
     * Helper function for constructing numeric literal storage
     * @tparam T data type
     * @param id identifier of the literal
     * @param value numeric value
     * @return Numeric LiteralStorage
     */
    template<typename T>
    constexpr auto makeNumericLit(var_t id, T value) noexcept(std::is_nothrow_move_constructible_v<T>) {
        return LiteralStorage<T>(id, NumericValue<T>(std::move(value)));
    }

    /**
     * Helper function for constructing constraint literal storage
     * @tparam T data type
     * @param id identifier of the literal
     * @param semantic constraint information
     * @return Semantic LiteralStorage
     */
    template<typename T>
    constexpr auto makeSemanticLit(var_t id, info_t semantic) noexcept {
        return LiteralStorage<T>(id, semantic);
    }
}



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
    bool operator==(const Literal<T> &l) const;

    T value() const;
    info_t constraint() const;
    
    std::ostream& display(std::ostream &os) const;
    
    static var_t index(const bool sign, const var_t x) ;
    static var_t var(const var_t s);
    static bool sign(const var_t s);

    void setValue(T v);

    info_t _id_{0};
    std::variant<info_t,T> _data_;
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

template <typename T> void Literal<T>::setValue(T v) { _data_ = v; }

template <typename T> bool Literal<T>::sameVariable(Literal<T> l) const {
  return isNumeric() == l.isNumeric() and variable() == l.variable();
}

template <typename T> bool Literal<T>::operator==(const Literal<T> &l) const {
  if (isNumeric()) {
    if (not l.isNumeric())
      return false;
    return _id_ == l._id_ and value() == l.value();
  } else if (l.isNumeric())
    return false;
  return _id_ == l._id_;
}

template <typename T> info_t Literal<T>::index(const bool sign, const var_t x) {
    return 2 * x + sign;
}

template <typename T> var_t Literal<T>::var(const info_t x) { return x / 2; }

template <typename T> bool Literal<T>::sign(const info_t x) { return x & 1; }

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
      if (value() == Constant::Infinity<T>)
        os << "infinity";
      else if (value() == -Constant::Infinity<T>)
        os << "-infinity";
      else
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
