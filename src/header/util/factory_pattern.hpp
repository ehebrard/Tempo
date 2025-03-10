/**
 * @author Tim Luchterhand
 * @date 17.11.23.
 */

#ifndef TEMPO_FACTORY_PATTERN_HPP
#define TEMPO_FACTORY_PATTERN_HPP

#include <string>
#include <variant>

#define ESCAPE(...) __VA_ARGS__

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
#define FACTORY_TYPE7(FUNC, ARG1, ARG2, ARG3, ARG4, ARG5, ARG6, ARG7)          \
  FACTORY_TYPE1(ARG1), FACTORY_TYPE6(ARG2, ARG3, ARG4, ARG5, ARG6, ARG7)

#define MAKE_POLYMORPHIC_TYPE(TYPE_NAME, ...) using TYPE_NAME = std::variant<__VA_ARGS__>;

#define HOLDS_FOR_ALL(INSTANCE, CONCEPT, ...)                                  \
  template <typename> struct __##INSTANCE##_tester__ : std::false_type {};     \
                                                                               \
  template <typename... Args>                                                  \
  struct __##INSTANCE##_tester__<std::variant<Args...>> {                      \
    static constexpr bool value = (CONCEPT<Args, __VA_ARGS__> && ...);         \
  };

#define MAKE_P_FACTORY_PATTERN(NAME, P_TYPE, ...)                         \
  class FACTORY_ENTRY(NAME) final {                                       \
    using NAME##FactoryType = std::variant<FACTORY_ENTRY(__VA_ARGS__)>;   \
                                                                          \
  public:                                                                 \
    FACTORY_ENTRY(NAME)(const FACTORY_ENTRY(NAME) &) = delete;            \
    FACTORY_ENTRY(NAME)(FACTORY_ENTRY(NAME) &&) = delete;                 \
    FACTORY_ENTRY(NAME) &operator=(const FACTORY_ENTRY(NAME) &) = delete; \
    FACTORY_ENTRY(NAME) &operator=(FACTORY_ENTRY(NAME) &&) = delete;      \
    static auto getInstance() noexcept -> const FACTORY_ENTRY(NAME) & {   \
      static FACTORY_ENTRY(NAME) instance;                                \
      return instance;                                                    \
    }                                                                     \
    template<typename ...Args>                                            \
    auto create(const std::string &typeName, Args &&...args ) const       \
        -> P_TYPE {                                                       \
      if (not registry.contains(typeName)) {                              \
        throw std::runtime_error("unknown type " + typeName);             \
      }                                                                   \
      const auto &constructor = registry.at(typeName);                    \
      return std::visit(                                                  \
          [&args...](const auto &ctor) -> P_TYPE {                        \
            return ctor.create(std::forward<Args>(args)...);              \
          },                                                              \
          constructor);                                                   \
    }                                                                     \
                                                                          \
  private:                                                                \
    FACTORY_ENTRY(NAME)() = default;                                      \
    std::unordered_map<std::string, NAME##FactoryType> registry{          \
        TYPE_ENTRY(__VA_ARGS__)};                                         \
  };

#define MAKE_FACTORY_PATTERN(TYPE_NAME, ...) MAKE_P_FACTORY_PATTERN(TYPE_NAME, TYPE_NAME, __VA_ARGS__)

#define MAKE_DEFAULT_FACTORY(TYPE, ...)                          \
struct FACTORY_ENTRY(TYPE) {                                     \
    static TYPE create(__VA_ARGS__) noexcept { return TYPE{}; }  \
};

#define MAKE_DEFAULT_TEMPLATE_FACTORY(TYPE, T_ARGS, ...)         \
struct FACTORY_ENTRY(TYPE) {                                     \
    template<T_ARGS>                                             \
    static TYPE create(__VA_ARGS__) noexcept { return TYPE{}; }  \
};


#define MAKE_P_FACTORY(TYPE, P_TYPE, ...)   \
struct FACTORY_ENTRY(TYPE) {                \
    static P_TYPE create(__VA_ARGS__)

#define MAKE_TEMPLATE_P_FACTORY(TYPE, P_TYPE, T_ARG, ARG)   \
  struct FACTORY_ENTRY(TYPE) {                              \
    template <T_ARG> static P_TYPE create(ARG)

#define MAKE_FACTORY(TYPE, ...) MAKE_P_FACTORY(TYPE, auto, __VA_ARGS__)

#define MAKE_TEMPLATE_FACTORY(TYPE, T_ARG, ...) MAKE_TEMPLATE_P_FACTORY(TYPE, auto, ESCAPE(T_ARG), ESCAPE(__VA_ARGS__))

#define EMPTY

#define DYNAMIC_DISPATCH(FNAME, FARG, LCAPTURE, LARG, CVSPEC)                                               \
decltype(auto) FNAME(FARG) CVSPEC {                                                                         \
    return std::visit([LCAPTURE](CVSPEC auto &impl) -> decltype(auto) { return impl.FNAME(LARG); }, *this); \
}

#define DYNAMIC_DISPATCH_FORWARD(FNAME, CVSPEC)                             \
template<typename ...Args>                                                  \
decltype(auto) FNAME(Args &&...args) CVSPEC {                               \
    return std::visit([&args...](CVSPEC auto &impl) -> decltype(auto) {     \
        return impl.FNAME(std::forward<Args>(args)...);                     \
    }, *this);                                                              \
}

#define DYNAMIC_DISPATCH_VOID(FNAME, CVSPEC) DYNAMIC_DISPATCH(FNAME, EMPTY, EMPTY, EMPTY, CVSPEC)
#endif //TEMPO_FACTORY_PATTERN_HPP
