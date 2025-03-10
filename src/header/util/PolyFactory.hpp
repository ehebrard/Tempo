/**
* @author Tim Luchterhand
* @date 06.03.25
* @file PolyFactory.hpp
* @brief Implementation of factory pattern for template factories.
* @note This implementation was largely inspired by the answer to this question:
* https://stackoverflow.com/questions/53255776/factory-method-for-template-classes
* The code can be found here
* https://coliru.stacked-crooked.com/a/5102c7a3b83be97f
*/

#ifndef TEMPO_POLYFACTORY_HPP
#define TEMPO_POLYFACTORY_HPP

#include <type_traits>
#include <unordered_map>
#include <concepts>

namespace tempo {
    template<typename...>
    struct types_t {};

    template<typename>
    struct ctor_tag_t {};

    template<typename T>
    inline constexpr ctor_tag_t<T> ctor_tag{};

    namespace detail {

        template<template<typename...> typename PolyType, typename T, typename Tag, typename... CtorArgs>
        struct factory_interface {
            factory_interface() = default;

            factory_interface(const factory_interface &) = delete;

            factory_interface(factory_interface &&) = delete;

            factory_interface &operator=(const factory_interface &) = delete;

            factory_interface &operator=(factory_interface &&) = delete;

            virtual ~factory_interface() = default;

            virtual auto tagged_build(typename Tag::template type<T>, CtorArgs...) const -> PolyType<T> = 0;
        };

        template<template<typename...> typename PolyType, typename SupportedTypes, typename Tag, typename... CtorArgs>
        struct poly_factory_interface_impl;

        template<template<typename...> typename PolyType, typename... SupportedTypes,
            typename Tag, typename... CtorArgs>
        struct poly_factory_interface_impl<PolyType, types_t<SupportedTypes...>, Tag, CtorArgs...>
                : factory_interface<PolyType, SupportedTypes, Tag, CtorArgs...>... {
            using factory_interface<PolyType, SupportedTypes, Tag, CtorArgs...>::tagged_build...;
        };

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
        template<typename CRTP_Base, typename PolyFactoryInterface, template<typename...> typename PolyType, typename T,
            typename Tag, typename... CtorArgs>
        struct factory_impl : PolyFactoryInterface {
            auto tagged_build(typename Tag::template type<T> tag, CtorArgs... args) const -> PolyType<T> final {
                return static_cast<const CRTP_Base *>(this)->template build_impl<T>(
                    std::forward<decltype(tag)>(tag), std::forward<decltype(args)>(args)...);
            }

            using PolyFactoryInterface::build;
        };
#pragma GCC diagnostic pop

        template<template<typename...> typename PolyType, typename SupportedTypes,
            typename Tag, typename... CtorArgs>
        struct poly_factory_interface : poly_factory_interface_impl<PolyType, SupportedTypes, Tag, CtorArgs...> {
            template<typename T>
            auto build(typename Tag::template type<T> tag, CtorArgs... args) const -> PolyType<T> {
                return this->tagged_build(std::forward<decltype(tag)>(tag), std::forward<decltype(args)>(args)...);
            }
        };

        template<typename CRTP_Base, typename PolyFactoryInterface, template<typename...> typename PolyType,
            typename SupportedTypes, typename Tag, typename... CtorArgs>
        struct poly_factory_impl;

        template<typename CRTP_Base, typename PolyFactoryInterface, template<typename...> typename PolyType,
            typename T, typename Tag, typename... CtorArgs>
        struct poly_factory_impl<CRTP_Base, PolyFactoryInterface, PolyType, types_t<T>, Tag, CtorArgs...>
                : factory_impl<CRTP_Base, PolyFactoryInterface, PolyType, T, Tag, CtorArgs...> {
            using factory_impl<CRTP_Base, PolyFactoryInterface, PolyType, T, Tag, CtorArgs...>::tagged_build;
        };

        template<typename CRTP_Base, typename PolyFactoryInterface, template<typename...> typename PolyType,
            typename T1, typename T2, typename... SupportedTypes, typename Tag, typename... CtorArgs>
        struct poly_factory_impl<CRTP_Base, PolyFactoryInterface, PolyType,
                    types_t<T1, T2, SupportedTypes...>, Tag, CtorArgs...>
                : factory_impl<CRTP_Base, poly_factory_impl<CRTP_Base, PolyFactoryInterface, PolyType,
                    types_t<T2, SupportedTypes...>, Tag, CtorArgs...>, PolyType, T1, Tag, CtorArgs...> {
            using factory_impl<CRTP_Base, poly_factory_impl<CRTP_Base, PolyFactoryInterface, PolyType,
                types_t<T2, SupportedTypes...>, Tag, CtorArgs...>, PolyType, T1, Tag, CtorArgs...>::tagged_build;
        };

        template<typename, typename>
        struct contains_type_t : std::false_type {};

        template<typename T, typename... Ts>
        struct contains_type_t<T, types_t<Ts...>> {
            static constexpr bool value = (std::same_as<T, Ts> || ...);
        };

        template<typename T, typename Types>
        concept supported = contains_type_t<T, Types>::value;

        template<typename Tag>
        concept auto_tag = std::is_default_constructible_v<std::remove_cvref_t<Tag>> and
                           (not std::is_lvalue_reference_v<Tag> or std::is_const_v<std::remove_reference_t<Tag>>);

    }


    /**
     * Default value tag to select appropriate template factory
     */
    struct DefaultTag {
        template<typename T>
        using type = ctor_tag_t<T>;
    };

    /**
     * Tag that is passed by value + move
     * @tparam T template to use as tag
     */
    template<template<typename>typename T>
    struct ValueTag {
        template<typename U>
        using type = T<U>;
    };

    /**
     * Tag that is passed by mutable reference
     * @tparam T template to use as tag
     */
    template<template<typename>typename T>
    struct RefTag {
        template<typename U>
        using type = T<U>&;
    };

    /**
     * Tag that is passed as const reference
     * @tparam T template to use as tag
     */
    template<template<typename>typename T>
    struct ConstRefTag {
        template<typename U>
        using type = const T<U>&;
    };

    /**
     * CRTP alias used to implement a factory. Have your implementation derive from this
     * @tparam Impl Concrete implementation
     * @tparam PolyType Comon type returned from the factory
     * @tparam ValueTypes List with all supported types
     * @tparam Tag Type used to identify template factories that is passed as first argument
     * @tparam CtorArgs Argument types passed to the factories
     */
    template<typename Impl, template<typename...> typename PolyType, typename ValueTypes,
        typename Tag, typename... CtorArgs>
    using MakeFactory = detail::poly_factory_impl<Impl, detail::poly_factory_interface<PolyType,
        ValueTypes, Tag, CtorArgs...>, PolyType, ValueTypes, Tag, CtorArgs...>;


    /**
     * Main factory class that dispatches construction to registered factory classes.
     * @tparam PolyType Comon type returned from the factories
     * @tparam SupportedTypes List with all supported types
     * @tparam CtorKey Factory identifier key used for registering factories
     * @tparam Tag Type used to identify template factories that is passed as first argument
     * @tparam CtorArgs Argument types passed to the factories
     */
    template<template<typename...> typename PolyType, typename SupportedTypes, typename CtorKey,
        typename Tag, typename... CtorArgs>
    class PolyFactory {
        using FactoryImpl = detail::poly_factory_interface<PolyType, SupportedTypes, Tag, CtorArgs...>;
        std::unordered_map<CtorKey, const FactoryImpl *> registry;
        template<typename T>
        using TagArg = typename Tag::template type<T>;

        PolyFactory() = default;

    public:
        PolyFactory(const PolyFactory &) = delete;

        PolyFactory(PolyFactory &&) = delete;

        PolyFactory &operator=(const PolyFactory &) = delete;

        PolyFactory &operator=(PolyFactory &&) = delete;

        /**
         * Singleton access
         * @return reference to factory
         */
        static auto get() noexcept -> PolyFactory & {
            static PolyFactory factory;
            return factory;
        }

        /**
         * Register a factory implementation
         * @param key identifier
         * @param factory pointer to factory instance
         * @return true if successful insertion, false if key already exists
         */
        bool registerFactory(const CtorKey &key, const FactoryImpl *factory) {
            return registry.emplace(key, factory).second;
        }

        /**
         * Build method
         * @tparam T value type
         * @param tag tag used to identify the correct template factory
         * @param key factory identifier
         * @param args arguments passed to the factory implementation
         * @return unique_ptr to constructed PolyType
         */
        template<detail::supported<SupportedTypes> T>
        auto build(TagArg<T> tag, const CtorKey &key, CtorArgs... args) const -> PolyType<T> {
            auto *ctor = registry.at(key);
            return ctor->template build<T>(std::forward<decltype(tag)>(tag), std::forward<decltype(args)>(args)...);
        }

        /**
         * Build method that does not take the tag as argument and instead the template type
         * @tparam T value type
         * @param key factory identifier
         * @param args arguments passed to the factory implementation
         * @return unique_ptr to constructed PolyType
         */
        template<detail::supported<SupportedTypes> T> requires(detail::auto_tag<TagArg<T>>)
        auto build(const CtorKey &key, CtorArgs... args) const -> PolyType<T> {
            auto *ctor = registry.at(key);
            return ctor->template build<T>(TagArg<T>{}, std::forward<decltype(args)>(args)...);
        }
    };
}


#endif //TEMPO_POLYFACTORY_HPP
