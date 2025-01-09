/**
* @author Tim Luchterhand
* @date 09.01.25
* @file Lookup.hpp
* @brief Contiguous lookup table
*/

#ifndef LOOKUP_HPP
#define LOOKUP_HPP

#include <ranges>
#include <functional>
#include <vector>
#include <optional>
#include <algorithm>
#include <Iterators.hpp>

#include "traits.hpp"


namespace tempo {

    /**
     * Requirement for a key-projection
     */
    template<typename P, typename K>
    concept key_projection = std::invocable<P, const K &> and
                             std::integral<std::remove_cvref_t<std::invoke_result_t<P, const K &>>>;

    /**
     * @brief Projection functor that gets the id from a given key
     */
    struct IdProjection {
        template<typename T>
        constexpr auto operator()(const T& t) const {
            return t.id();
        }
    };

    /**
     * @brief Efficient lookup table for non-contiguous key-value pairs
     * @details @copybrief
     * @tparam K key type
     * @tparam V value type
     * @tparam P key projection function
     * @note this class has a memory overhead if keys are non-contiguous.
     */
    template<typename K, typename V, key_projection<K> P = std::identity>
    class Lookup {
        std::vector<V> table{};
        long offset{0};
        P projection{};

        template<concepts::typed_range<K> Keys, std::convertible_to<V> Value = V>
        void initTable(const Keys& keys, const Value &value = {}) {
            if (table.empty()) {
                auto [min, max] = std::ranges::minmax(keys, {}, projection);
                auto prMin = projection(min);
                auto prMax = projection(max);
                assert(prMin <= prMax);
                assert(prMin >= 0 and prMax >= 0);
                table.resize(static_cast<std::size_t>(prMax - prMin) + 1, value);
                offset = static_cast<long>(prMin);
            }
        }

    public:
        Lookup() = default;

        /**
         * Ctor
         * @tparam Keys key range type
         * @tparam Values value range type
         * @param keys range of keys
         * @param defaultValue Value to be inserted into empty places (default: default value of value type)
         * @param values optional range of values (if not given will be default initialized)
         * @param keyBounds optional pair(lower key bound, upper key bound), both inclusive. If not given, will be
         * determined form the key range
         * @param projection optional key projection function (default identity)
         */
        template<concepts::typed_range<K> Keys, concepts::ctyped_range<V> Values = std::vector<V>>
        explicit Lookup(const Keys &keys, const V &defaultValue = {}, const Values &values = {},
                        std::optional<std::pair<long, long>> keyBounds = {}, const P &projection = {})
            : table(keyBounds.has_value() ? static_cast<std::size_t>(keyBounds->second - keyBounds->first) : 0,
                    defaultValue),
              offset(keyBounds.has_value() ? keyBounds->first : 0), projection(projection) {
            initTable(keys, defaultValue);
            for (auto [k, v]: iterators::zip(keys, values)) {
                this->at(k) = v;
            }
        }

        /**
         * Value access without bounds checking
         * @param key Key
         * @return stored value
         */
        decltype(auto) operator[](const K &key) const {
            return table[projection(key) - offset];
        }

        /**
         * @copydoc operator[](const K &key)
         */
        decltype(auto) operator[](const K &key) {
            return table[projection(key) - offset];
        }

        /**
         * Whether a key is contained in the lookup table
         * @param key key to check
         * @return true if key is contained, false otherwise
         */
        bool contains(const K &key) const noexcept {
            auto pr = static_cast<long>(projection(key));
            return pr >= offset and pr <= maxKey();
        }

        /**
         * Value access with bounds checking
         * @param key Key
         * @return stored value
         * @throws std::out_of_range
         */
        decltype(auto) at(const K &key) {
            auto pr = static_cast<long>(projection(key));
            if (pr < offset or pr > maxKey()) {
                throw std::out_of_range(
                    "Lookup::at out of range: key is " + std::to_string(pr) + ", offset is " +
                    std::to_string(offset) + ", max key is " + std::to_string(maxKey()));
            }

            return this->operator[](key);
        }

        /**
         * @copydoc at(const K &key)
         */
        decltype(auto) at(const K &key) const {
            return std::as_const(traits::as_mut(*this).at(key));
        }

        /**
         * Direct access to stored values
         * @return
         */
        auto &data() noexcept {
            return table;
        }

        /**
         * @copydoc data()
         */
        const auto &data() const noexcept {
            return table;
        }

        /**
         * Number of entries
         * @return number of entries
         * @note also includes uninitialized values
         */
        [[nodiscard]] std::size_t size() const noexcept {
            return table.size();
        }

        /**
         * Gets the lowest key stored
         * @return projected value of lowest key stored
         */
        [[nodiscard]] std::size_t keyOffset() const noexcept {
            return offset;
        }

        /**
         * Gets the highest key stored
         * @return projected value of highest key stored
         */
        [[nodiscard]] long maxKey() const noexcept {
            return static_cast<long>(table.size() + offset) - 1;
        }
    };

    template<std::ranges::range Keys, typename T, concepts::ctyped_range<T> Values = std::vector<T>, typename P =
        std::identity>
    Lookup(const Keys &, const T & = {}, const Values & = {}, std::optional<std::pair<long, long>>  = {},
           const P & = {}) -> Lookup<std::ranges::range_value_t<Keys>, T, P>;

}

#endif //LOOKUP_HPP
