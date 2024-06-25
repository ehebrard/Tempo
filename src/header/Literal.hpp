/************************************************
 * Tempo Solver.hpp
 *
 * Copyright 2024 Emmanuel Hebrard and Tim Luchterhand
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

#ifndef __TEMPO_LITERAL_HPP
#define __TEMPO_LITERAL_HPP

#include <variant>

#include "Constant.hpp"
#include "Global.hpp"
#include "util/traits.hpp"

namespace tempo {


/**********************************************
* Literal
**********************************************/

enum bound { lower = 0, upper = 1 };

namespace detail {
    /**
     * Tag class to resolve ctor ambiguity
     */
    struct Numeric {};

    /**
     * Tag class to resolve ctor ambiguity
     */
    struct Boolean {};

    /**
     * Tagged union class that holds a variable ID and either numeric data or constraint info
     * @tparam T data ype
     */
    template<concepts::scalar T>
    class LiteralStorage {
        static constexpr info_t SignBit = 2;
        static constexpr info_t TagBit = 1;
        info_t info;
        union {
            T numericData;
            info_t semanticInfo;

        };

    public:
        /**
         * Ctor
         * constructs a boolean literal storage
         * @param sign sign of the literal
         * @param literalId identifier of the literal (not of the variable!)
         * @param semanticInfo semantic information
         * @note the id of the corresponding variable should be litId / 2. See also makeSemanticLit
         */
        constexpr LiteralStorage(info_t literalId, info_t semanticInfo, Boolean) noexcept:
                info(literalId << 1), semanticInfo(semanticInfo) {}

        /**
         * Ctor
         * constructs a numeric literal storage
         * @param sign sign of the literal
         * @param literalId identifier of the literal (not of the variable!)
         * @param numericVal numeric value stored in the literal
         * @note see also makeNumericLit
         */
        constexpr LiteralStorage(info_t literalId, T numericVal, Numeric) noexcept: info(
                (literalId << 1) | TagBit), numericData(std::move(numericVal)) {}

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
        [[nodiscard]] constexpr info_t id() const noexcept {
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
         * Access to the stored value without safe-guard
         * @return the stored value
         * @note calling this function on a boolean literal is undefined behavior
         */
        [[nodiscard]] constexpr T value_unsafe() const noexcept {
            return numericData;
        }

        /**
         * Access to the stored semantic information
         * @return the stored semantic information
         * @throws std::bad_variant_access if storage does not contain semantic information
         */
        [[nodiscard]] constexpr info_t semantic() const {
            if (isNumeric()) {
                throw std::bad_variant_access();
            }

            return semanticInfo;
        }

        /**
         * Access to the stored semantic information without safe-guard
         * @return the stored semantic information
         * @note calling this function on a non-boolean literal is undefined behavior
         */
        [[nodiscard]] constexpr info_t semantic_unsafe() const noexcept {
            return semanticInfo;
        }

        /**
         * set the stored value, implicitly transforming the literal into a numeric literal
         */
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
    constexpr auto makeNumericLit(info_t id, T value) noexcept(std::is_nothrow_move_constructible_v<T>) {
        return LiteralStorage<T>(id, std::move(value), Numeric{});
    }

    /**
     * Helper function for constructing semantic literal storage
     * @tparam T data type
     * @param id identifier of the literal
     * @param semantic semantic information
     * @return Semantic LiteralStorage
     */
    template<typename T>
    constexpr auto makeSemanticLit(info_t id, info_t semantic) noexcept {
        return LiteralStorage<T>(id, semantic, Boolean{});
    }
}



template<typename T>
struct Literal : public detail::LiteralStorage<T> {
  constexpr Literal() noexcept
      : detail::LiteralStorage<T>(Constant::NoVar, Constant::NoSemantic,
                                  detail::Boolean{}) {}

  /**
   * Ctor
   * constructs a NUMERIC literal
   * @param sign sign of the literal
   * @param x variable identifier
   * @param v numeric value
   * @note see also makeNumericLiteral()
   */
  constexpr Literal(bool sign, var_t x, T v, detail::Numeric) noexcept;

  /**
   * Ctor
   * constructs a BOOLEAN literal
   * @param sign sign of the literal
   * @param x variable identifier
   * @param v semantic information
   * @note see also makeBooleanLiteral()
   */
  constexpr Literal(bool sign, var_t x, info_t info, detail::Boolean) noexcept;
  using detail::LiteralStorage<T>::LiteralStorage;

  [[nodiscard]] constexpr bool isBoolean() const noexcept;
  [[nodiscard]] constexpr bool hasSemantic() const noexcept;

  constexpr operator info_t() const noexcept;

  [[nodiscard]] constexpr var_t variable() const noexcept;

  constexpr bool sameVariable(const Literal<T> &l) const noexcept;
  constexpr bool operator==(const Literal<T> &l) const noexcept;

  [[nodiscard]] constexpr info_t constraint() const;
  std::ostream &display(std::ostream &os) const;

  static constexpr info_t index(bool sign, var_t x) noexcept;
  static constexpr var_t var(info_t l) noexcept;
  static constexpr bool sgn(info_t l) noexcept;
};

template <typename T>
constexpr Literal<T>::Literal(bool sign, var_t x, T v, detail::Numeric) noexcept:
        detail::LiteralStorage<T>(index(sign, x), std::move(v), detail::Numeric{}) {}

template <typename T>
constexpr Literal<T>::Literal(bool sign, var_t x, info_t info, detail::Boolean) noexcept:
        detail::LiteralStorage<T>(index(sign, x), info, detail::Boolean{}) {}


/**
 * Helper function to generate numeric Literals
 * @param sign sign of the literal
 * @param variableId id of the variable (not the index of the literal)
 * @param value numeric value to store
 * @note prefer using this function when constructing numeric literals over the constructor
 */
template<concepts::scalar T>
constexpr auto makeNumericLiteral(bool sign, var_t variableId, T value) noexcept {
    return Literal<T>(sign, variableId, std::move(value), detail::Numeric{});
}

/**
 * Helper function to generate boolean Literals
 * @param sign sign of the literal
 * @param variableId id of the variable (not the index of the literal)
 * @param info semantic info
 * @note prefer using this function when constructing boolean literals over the constructor
 */
template<concepts::scalar T>
constexpr auto makeBooleanLiteral(bool sign, var_t variableId, info_t info) noexcept {
    return Literal<T>(sign, variableId, info, detail::Boolean{});
}

template <typename T>
constexpr bool Literal<T>::sameVariable(const Literal<T> &l) const noexcept{
  return this->isNumeric() == l.isNumeric() and variable() == l.variable();
}

template <typename T>
constexpr bool Literal<T>::operator==(const Literal<T> &l) const noexcept {
    if (this->isNumeric()) {
        return l.isNumeric() and this->id() == l.id() and this->value_unsafe() == l.value_unsafe();
    } else {
        return not l.isNumeric() and this->id() == l.id();
    }
}

template <typename T>
constexpr info_t Literal<T>::index(const bool sign, const var_t x) noexcept {
    return 2 * x + sign;
}

template <typename T>
constexpr var_t Literal<T>::var(info_t l) noexcept { return l / 2; }

template<typename T>
constexpr bool Literal<T>::sgn(info_t l) noexcept {
    return l & 1;
}

template <typename T>
constexpr Literal<T>::operator info_t() const noexcept { return this->id(); }

template<typename T>
constexpr var_t Literal<T>::variable() const noexcept{
    return this->id() / 2;
}

template<typename T>
constexpr info_t Literal<T>::constraint() const {
    return this->semantic() + this->sign();
}

template<typename T>
constexpr bool Literal<T>::isBoolean() const noexcept {
    return not this->isNumeric();
}

template<typename T>
constexpr bool Literal<T>::hasSemantic() const noexcept {
    return this->isNumeric() or this->semantic_unsafe() != Constant::NoSemantic;
}

template<typename T>
Literal<T> operator~(const Literal<T> &l) noexcept {
    if (l.isNumeric()) {
        return makeNumericLiteral(not l.sign(), l.variable(), -l.value_unsafe() - Gap<T>::epsilon());
    } else {
        return Literal<T>(l.id() ^ 1, l.semantic_unsafe(), detail::Boolean{});
    }
}

template <typename T> std::ostream& Literal<T>::display(std::ostream &os) const {
  if (this->id() == Constant::NoVar) {
    if (this->value() == Constant::Infinity<T>)
      os << "infinity";
    else if (this->value() == -Constant::Infinity<T>)
      os << "-infinity";
    else
      os << "constant: " << this->value();
  } else if (this->isNumeric()) {
    os << (this->sign() ? "x" : "-x") << variable() << " <= " << this->value();
  } else {
    os << (this->sign() ? "b" : "¬b") << variable()
       << (hasSemantic() ? "*" : "");
  }
    return os;
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const Literal<T> &x) {
  return x.display(os);
}

template <typename T>
constexpr Literal<T> ub(const var_t x) noexcept {
  return makeNumericLiteral(true, x, Constant::Infinity<T>);
}

template <typename T>
constexpr Literal<T> lb(const var_t x) noexcept {
  return makeNumericLiteral(false, x, Constant::Infinity<T>);
}

template<typename T>
constexpr Literal<T> leq(const var_t x, const T v) noexcept {
    return makeNumericLiteral(true, x, v);
}

template<typename T>
constexpr Literal<T> gt(const var_t x, const T v) noexcept {
    return makeNumericLiteral(false, x, -v-Gap<T>::epsilon());
}

template<typename T>
constexpr Literal<T> lt(const var_t x, const T v) noexcept {
    return makeNumericLiteral(true, x, v-Gap<T>::epsilon());
}

template<typename T>
constexpr Literal<T> geq(const var_t x, const T v) noexcept {
    return makeNumericLiteral(false, x, -v);
}

template<typename T>
constexpr Literal<T> pos(const var_t x, const info_t d=0) noexcept {
    return makeBooleanLiteral<T>(true, x, d);
}

template<typename T>
constexpr Literal<T> neg(const var_t x, const info_t d=0) noexcept {
    return makeBooleanLiteral<T>(false, x, d);
}


}

#endif // __Literal_HPP
