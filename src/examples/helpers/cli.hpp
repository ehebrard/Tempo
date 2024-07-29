/**
* @author Tim Luchterhand
* @date 24.07.24
* @brief contains command line helpers
*/

#ifndef TEMPO_CLI_HPP
#define TEMPO_CLI_HPP

#include <string>
#include <utility>
#include <concepts>
#include <optional>

#include "util/traits.hpp"
#include "util/Options.hpp"

/**
 * @brief Command line interface helpers
 */
namespace cli {

/**
 * @brief CLI argument specification
 * @tparam T type of value argument
 */
    template<typename T>
    struct ArgSpec {
        /**
         * CTor
         * @param argName long name of the flag
         * @param explanation explanation
         * @param required whether the argument is required
         * @param destination reference to field where the value is stored
         */
        ArgSpec(std::string argName, std::string explanation, bool required, T &destination) noexcept
                : argName(std::move(argName)), explanation(std::move(explanation)), required(required),
                  destination(destination) {}

        /**
         * CTor
         * @tparam U
         * @param argName long name of the flag
         * @param explanation explanation
         * @param required whether the argument is required
         * @param destination reference to field where the value is stored
         * @param defaultVal default value for the argument
         */
        template<typename U>
        ArgSpec(std::string argName, std::string explanation, bool required, T &destination, U &&defaultVal) noexcept
                : argName(std::move(argName)), explanation(std::move(explanation)), required(required),
                  destination(destination), defaultVal(std::forward<U>(defaultVal)) {}

        std::string argName;
        std::string explanation;
        bool required;
        T &destination;
        std::optional<T> defaultVal{};
    };

/**
 * @brief @copydoc ArgSpec
 * @details @copybrief
 * specification for switch args
 */
    template<>
    struct ArgSpec<bool> {
        /**
         * CTor
         * @param argName long name of the flag
         * @param explanation explanation
         * @param destination reference to field where the value is stored
         * @param defaultVal default value for the argument
         */
        ArgSpec(std::string argName, std::string explanation, bool &destination, bool defaultVal) noexcept
                : argName(std::move(argName)), explanation(std::move(explanation)),
                  destination(destination), defaultVal(defaultVal) {}

        std::string argName;
        std::string explanation;
        bool &destination;
        bool defaultVal;
    };

    using SwitchSpec = ArgSpec<bool>;

    namespace detail {
        void configureParser(tempo::Parser &) {}

        template<typename T, typename ...Ts>
        void configureParser(tempo::Parser &parser, const ArgSpec <T> &arg, const Ts &...rest) {
            if constexpr (std::same_as<T, bool>) {
                parser.getCmdLine().add<TCLAP::SwitchArg>(arg.destination, "", arg.argName, arg.explanation,
                                                          arg.defaultVal);
            } else {
                parser.getCmdLine().add<TCLAP::ValueArg<T>>(arg.destination, "", arg.argName, arg.explanation,
                                                            arg.required, arg.defaultVal.value_or(T{}),
                                                            typeid(T).name());
            }

            configureParser(parser, rest...);
        }
    }

/**
 * Parses all default options supported by the solver from from command line. Supports adding additional options
 * @param argc argument counter
 * @param argv cli arguments
 * @param argSpecs additional argument specifications
 * @return parsed options
 */
    template<tempo::concepts::same_template<ArgSpec> ...Ts>
    auto parseOptions(int argc, char **argv, const Ts &...argSpecs) -> tempo::Options {
        tempo::Parser p = tempo::getBaseParser();
        detail::configureParser(p, argSpecs...);
        p.parse(argc, argv);
        return p.getOptions();
    }
}


#endif //TEMPO_CLI_HPP
