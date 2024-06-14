
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
         * @param semanticInfo semantic information
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
     * Helper function for constructing semantic literal storage
     * @tparam T data type
     * @param id identifier of the literal
     * @param semantic semantic information
     * @return Semantic LiteralStorage
     */
    template<typename T>
    constexpr auto makeSemanticLit(var_t id, info_t semantic) noexcept {
        return LiteralStorage<T>(id, semantic);
    }
}



template<typename T>
struct Literal : public detail::LiteralStorage<T> {
    constexpr Literal() noexcept: detail::LiteralStorage<T>(Constant::NoVarx, Constant::NoSemantic) {}
    constexpr Literal(bool sign, var_t x, T v) noexcept;
    constexpr Literal(bool sign, var_t x, info_t v) noexcept;
    using detail::LiteralStorage<T>::LiteralStorage;

    [[nodiscard]] constexpr bool isBoolean() const noexcept;
    [[nodiscard]] constexpr bool hasSemantic() const noexcept;

    //    operator int() const;
    constexpr operator std::uint32_t() const noexcept;

    [[nodiscard]] constexpr var_t variable() const noexcept;

    constexpr bool sameVariable(const Literal<T> &l) const noexcept;
    constexpr bool operator==(const Literal<T> &l) const noexcept;

    [[nodiscard]] constexpr info_t constraint() const noexcept;
    std::ostream& display(std::ostream &os) const;
    
    static constexpr std::uint32_t index(bool sign, var_t x) noexcept;
    static constexpr var_t var(info_t l) noexcept;
    static constexpr bool sgn(std::uint32_t l) noexcept;
};

template <typename T>
constexpr Literal<T>::Literal(bool sign, var_t x, T v) noexcept:
        detail::LiteralStorage<T>(index(sign, x), detail::NumericValue<T>(std::move(v))) {}

template <typename T>
constexpr Literal<T>::Literal(bool sign, var_t x, const info_t v) noexcept:
        detail::LiteralStorage<T>(index(sign, x), v) {}

template <typename T>
constexpr bool Literal<T>::sameVariable(const Literal<T> &l) const noexcept{
  return this->isNumeric() == l.isNumeric() and variable() == l.variable();
}

template <typename T>
constexpr bool Literal<T>::operator==(const Literal<T> &l) const noexcept {
    if (this->isNumeric()) {
        return l.isNumeric() and this->id() == l.id() and this->value() == l.value();
    } else {
        return not l.isNumeric() and this->id() == l.id();
    }
}

template <typename T>
constexpr std::uint32_t Literal<T>::index(const bool sign, const var_t x) noexcept {
    return 2 * x + sign;
}

template <typename T>
constexpr var_t Literal<T>::var(info_t l) noexcept { return l / 2; }

template<typename T>
constexpr bool Literal<T>::sgn(std::uint32_t l) noexcept {
    return l & 1;
}

template <typename T>
constexpr Literal<T>::operator std::uint32_t() const noexcept { return this->id(); }

template<typename T>
constexpr var_t Literal<T>::variable() const noexcept{
    return this->id() / 2;
}

template<typename T>
constexpr info_t Literal<T>::constraint() const noexcept {
    return this->semantic() + this->sign();
}

template<typename T>
constexpr bool Literal<T>::isBoolean() const noexcept {
    return not this->isNumeric();
}

template<typename T>
constexpr bool Literal<T>::hasSemantic() const noexcept {
    return this->isNumeric() or this->semantic() != Constant::NoSemantic;
}

template <typename T> Literal<T> operator~(const Literal<T> &l) {
  if (l.isNumeric()) {
    return Literal<T>(not l.sign(), l.variable(), -l.value() - Gap<T>::epsilon());
  } else {
    return Literal<T>(l.id() ^ 1, l.semantic());
  }
}

template <typename T> std::ostream& Literal<T>::display(std::ostream &os) const {
    if(this->id() == Constant::NoVarx) {
      if (this->value() == Constant::Infinity<T>)
        os << "infinity";
      else if (this->value() == -Constant::Infinity<T>)
        os << "-infinity";
      else
        os << "constant: " << this->value();
    } else if(this->isNumeric()) {
        os << (this->sign() ? "x" : "-x") << variable() << " <= " << this->value();
    } else {
      os << (this->sign() ? "b" : "¬b") << variable() << (hasSemantic() ? "*" : "");
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
