/**
* @author Tim Luchterhand
* @date 13.06.24
* @brief
*/

#ifndef TEMPO_LITERALSTORAGE_HPP
#define TEMPO_LITERALSTORAGE_HPP

#include <variant>

#include "Global.hpp"
#include "traits.hpp"

namespace tempo {

    template<typename T>
    struct NumericValue {
        explicit constexpr NumericValue(T val) noexcept: val(val) {}

        T val;
    };

    /**
     * Tagged union class that holds a variable ID and either numeric data or semantic info
     * @tparam T data ype
     */
    template<typename T>
    class LiteralStorage {
        std::uint32_t info;
        union {
            T numericData;
            info_t semanticInfo;

        };

    public:
        /**
         * Ctor
         * constructs a semantic literal
         * @param sign sign of the literal
         * @param literalId identifier of the literal (not of the variable!)
         * @param semanticInfo semantic information
         * @note the id of the corresponding variable should be litId / 2. See also makeSemanticLit
         */
        constexpr LiteralStorage(bool sign, var_t literalId, info_t semanticInfo) noexcept: info(
                (literalId << 2) | sign), semanticInfo(semanticInfo) {}

        /**
         * Ctor
         * constructs a numeric literal
         * @param sign sign of the literal
         * @param literalId identifier of the literal (not of the variable!)
         * @param numericVal numeric value stored in the literal
         * @note see also makeNumericLit
         */
        constexpr LiteralStorage(bool sign, var_t literalId, NumericValue<T> numericVal) noexcept: info(
                (literalId << 2) | sign | 0x0002), numericData(std::move(numericVal.val)) {}

        /**
         * Checks whether the storage holds a numeric value
         * @return true if storage holds a numeric value, false otherwise
         */
        [[nodiscard]] constexpr bool isNumeric() const noexcept {
            return info & 0x0002;
        }

        /**
         * Gets the sign of the literal
         * @return returns the sign of the literal
         */
        [[nodiscard]] constexpr bool sign() const noexcept {
            return info & 0x0001;
        }

        /**
         * Gets the id of the literal
         * @return the literal's id
         */
        [[nodiscard]] constexpr std::uint32_t id() const noexcept {
            return (info & 0xFFFC) >> 2;
        }

        /**
         * Access to the stored value
         * @return the stored value
         * @throws std::bad_variant_access if storage does not contain numeric value
         */
        constexpr T &value() {
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
        constexpr info_t &semantic() {
            if (isNumeric()) {
                throw std::bad_variant_access();
            }

            return semanticInfo;
        }


        /**
         * const value access
         * @copydoc value
         */
        [[nodiscard]] constexpr const T &value() const {
            return const_cast<LiteralStorage *>(this)->value();
        }

        /**
         * const semantic access
         * @copydoc value
         */
        [[nodiscard]] constexpr const info_t &semantic() const {
            return const_cast<LiteralStorage *>(this)->semantic();
        }

    };

    /**
     * Helper function for constructing numeric literal storage
     * @tparam T data type
     * @param sign sign of the literal
     * @param id identifier of the literal
     * @param value numeric value
     * @return Numeric LiteralStorage
     */
    template<typename T>
    constexpr auto makeNumericLit(bool sign, var_t id, T value) noexcept(std::is_nothrow_move_constructible_v<T>) {
        return LiteralStorage<T>(sign, id, NumericValue<T>(std::move(value)));
    }

    /**
     * Helper function for constructing semantic literal storage
     * @tparam T data type
     * @param sign sign of the literal
     * @param id identifier of the literal
     * @param semantic semantic information
     * @return Semantic LiteralStorage
     */
    template<typename T>
    constexpr auto makeSemanticLit(bool sign, var_t id, info_t semantic) noexcept {
        return LiteralStorage<T>(sign, id, semantic);
    }

}

#endif //TEMPO_LITERALSTORAGE_HPP
