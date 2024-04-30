/**
 * @author Tim Luchterhand
 * @date 15.11.22.
 */

#ifndef TEMPO_TRAITS_HPP
#define TEMPO_TRAITS_HPP
#include <type_traits>
#include <functional>
#include "Global.hpp"

/**
 * @brief namespace containing various type traits
 */
namespace tempo::traits {
    template<typename Container, typename = std::void_t<>>
    struct value_type {};

    template<typename Container>
    struct value_type<Container, std::void_t<typename Container::value_type>> {
        using type = typename Container::value_type;
    };

    template<typename T>
    using value_type_t = typename value_type<T>::type;

    template<typename T, typename = std::void_t<>>
    struct is_random_access_iterator {
        static constexpr bool value = false;
    };

    template<typename T>
    struct is_random_access_iterator<T, std::void_t<typename std::iterator_traits<T>::iterator_category>> {
        static constexpr bool value = std::is_base_of_v<std::random_access_iterator_tag,
                typename std::iterator_traits<T>::iterator_category>;
    };

    template<typename T, typename = std::void_t<>>
    struct is_random_accessible : std::false_type {};

    template<typename T>
    struct is_random_accessible<T, std::void_t<typename T::iterator>> {
        static constexpr bool value = is_random_access_iterator<typename T::iterator>::value;
    };

    template<typename T>
    constexpr inline bool is_random_accessible_v = is_random_accessible<T>::value;

    template<typename T, typename = std::void_t<>>
    struct is_iterable : std::false_type {};

    template<typename T>
    struct is_iterable<T, std::void_t<decltype(std::begin(std::declval<T>()), std::end(std::declval<T>()))>>
            : std::true_type {};

    template<typename T>
    constexpr inline bool is_iterable_v = is_iterable<T>::value;

    template<template<typename...> typename Template, typename T>
    struct is_same_template : std::false_type {};

    template<template<typename...> typename Template, typename... Args>
    struct is_same_template<Template, Template<Args...>> : std::true_type {};

    template<template<typename ...> typename Template, typename T>
    constexpr inline bool is_same_template_v = is_same_template<Template, T>::value;
}

namespace tempo::concepts {
    template<typename Functor, typename Return, typename ...Args>
    concept callable_r = std::invocable<Functor, Args...> &&
                         std::convertible_to<std::invoke_result_t<Functor, Args...>, Return>;

    template<typename T>
    concept scalar = std::integral<T> || std::floating_point<T>;

    template<typename E, typename T>
    concept event_dist_fun = callable_r<E, T, event, event>;

    template<typename E>
    concept arbitrary_event_dist_fun = std::invocable<E, event, event> &&
                                       scalar<std::remove_reference_t<std::invoke_result_t<E, event, event>>>;

    template<typename R, typename T>
    concept typed_range = std::ranges::range<R> && std::same_as<std::ranges::range_value_t<R>, T>;

    template<typename R, template<typename> typename T>
    concept ttyped_range = std::ranges::range<R> && traits::is_same_template_v<T, std::ranges::range_value_t<R>>;

    template<typename T>
    concept printable = requires(const T &obj, std::ostream &os) {
        os << obj;
    };
}

#endif //TEMPO_TRAITS_HPP
