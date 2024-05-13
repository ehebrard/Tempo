/**
 * @author Tim Luchterhand
 * @date 17.11.23.
 */

#ifndef TEMPO_FACTORY_PATTERN_HPP
#define TEMPO_FACTORY_PATTERN_HPP

#include <string>
#include <variant>

#define GET_MACRO(_1, _2, _3, _4, _5, _6, NAME, ...) NAME
#define TYPE_ENTRY(...) GET_MACRO(__VA_ARGS__, \
TYPE_ENTRY6, TYPE_ENTRY5, TYPE_ENTRY4, TYPE_ENTRY3, TYPE_ENTRY2, TYPE_ENTRY1, TYPE_ENTRY0)(__VA_ARGS__)
#define TYPE_ENTRY1(ARG1) {#ARG1, ARG1##Factory{}}
#define TYPE_ENTRY2(ARG1, ARG2) TYPE_ENTRY1(ARG1), TYPE_ENTRY1(ARG2)
#define TYPE_ENTRY3(ARG1, ARG2, ARG3) TYPE_ENTRY1(ARG1), TYPE_ENTRY2(ARG2, ARG3)
#define TYPE_ENTRY4(ARG1, ARG2, ARG3, ARG4) TYPE_ENTRY1(ARG1), TYPE_ENTRY3(ARG2, ARG3, ARG4)
#define TYPE_ENTRY5(ARG1, ARG2, ARG3, ARG4, ARG5) TYPE_ENTRY1(ARG1), TYPE_ENTRY4(ARG2, ARG3, ARG4, ARG5)
#define TYPE_ENTRY6(FUNC, ARG1, ARG2, ARG3, ARG4, ARG5, ARG6) TYPE_ENTRY1(ARG1), TYPE_ENTRY5(ARG2, ARG3, ARG4, ARG5, ARG6)

#define FACTORY_ENTRY(...) GET_MACRO(__VA_ARGS__, \
FACTORY_TYPE6, FACTORY_TYPE5, FACTORY_TYPE4, FACTORY_TYPE3, FACTORY_TYPE2, FACTORY_TYPE1, FACTORY_TYPE0)(__VA_ARGS__)
#define FACTORY_TYPE1(ARG1) ARG1##Factory
#define FACTORY_TYPE2(ARG1, ARG2) FACTORY_TYPE1(ARG1), FACTORY_TYPE1(ARG2)
#define FACTORY_TYPE3(ARG1, ARG2, ARG3) FACTORY_TYPE1(ARG1), FACTORY_TYPE2(ARG2, ARG3)
#define FACTORY_TYPE4(ARG1, ARG2, ARG3, ARG4) FACTORY_TYPE1(ARG1), FACTORY_TYPE3(ARG2, ARG3, ARG4)
#define FACTORY_TYPE5(ARG1, ARG2, ARG3, ARG4, ARG5) FACTORY_TYPE1(ARG1), FACTORY_TYPE4(ARG2, ARG3, ARG4, ARG5)
#define FACTORY_TYPE6(FUNC, ARG1, ARG2, ARG3, ARG4, ARG5, ARG6) FACTORY_TYPE1(ARG1), FACTORY_TYPE5(ARG2, ARG3, ARG4, ARG5, ARG6)

#define MAKE_POLYMORPHIC_TYPE(TYPE_NAME, ...) using TYPE_NAME = std::variant<__VA_ARGS__>;


#define HOLDS_FOR_ALL(INSTANCE, CONCEPT, ...)                                                                   \
    template<typename>                                                                                          \
    struct __##INSTANCE##_tester__ : std::false_type {};                                                        \
                                                                                                                \
    template<typename ...Args>                                                                                  \
    struct __##INSTANCE##_tester__<std::variant<Args...>> {                                                     \
        static constexpr bool value = (CONCEPT<Args, __VA_ARGS__> && ...);                                      \
    };                                                                                                          \

#define MAKE_T_FACTORY_PATTERN(TYPE_NAME, T_HEADER, CTOR_ARG, ...)                                                     \
    MAKE_POLYMORPHIC_TYPE(TYPE_NAME, __VA_ARGS__)                                                                      \
class TYPE_NAME##Factory final {                                                                                       \
    using TYPE_NAME##FactoryType = std::variant<FACTORY_ENTRY(__VA_ARGS__)>;                                           \
public:                                                                                                                \
    TYPE_NAME##Factory(const TYPE_NAME##Factory&) = delete;                                                            \
    TYPE_NAME##Factory(TYPE_NAME##Factory&&) = delete;                                                                 \
    TYPE_NAME##Factory &operator=(const TYPE_NAME##Factory&) = delete;                                                 \
    TYPE_NAME##Factory &operator=(TYPE_NAME##Factory&&) = delete;                                                      \
    static auto getInstance() noexcept -> const TYPE_NAME##Factory& {                                                  \
        static TYPE_NAME##Factory instance;                                                                            \
        return instance;                                                                                               \
    }                                                                                                                  \
    T_HEADER                                                                                                           \
    auto create(const std::string &typeName, const CTOR_ARG &arguments) const -> TYPE_NAME {                           \
        const auto &constructor = registry.at(typeName);                                                               \
        return std::visit([&arguments](const auto&ctor) -> TYPE_NAME { return ctor.create(arguments); }, constructor); \
    }                                                                                                                  \
private:                                                                                                               \
    TYPE_NAME##Factory() = default;                                                                                    \
    std::unordered_map<std::string, TYPE_NAME##FactoryType> registry{TYPE_ENTRY(__VA_ARGS__)};                         \
};

#define MAKE_DEFAULT_FACTORY(TYPE, ...)                          \
struct TYPE##Factory {                                           \
    static TYPE create(__VA_ARGS__) noexcept { return TYPE{}; }  \
};

#define MAKE_FACTORY_PATTERN(TYPE_NAME, CTOR_ARG, ...) MAKE_T_FACTORY_PATTERN(TYPE_NAME, , CTOR_ARG, __VA_ARGS__)

#define MAKE_FACTORY(TYPE, ...)         \
struct TYPE##Factory {                  \
    static TYPE create(__VA_ARGS__)

#define MAKE_TEMPLATE_FACTORY(TYPE, T_ARG, ARG) \
struct TYPE##Factory {                          \
    template<T_ARG>                             \
    static TYPE create(ARG)                     \

#endif //TEMPO_FACTORY_PATTERN_HPP
