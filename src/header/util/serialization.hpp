/**
 * @author Tim Luchterhand
 * @date 21.03.23.
 */

#ifndef TEMPO_SERIALIZATION_HPP
#define TEMPO_SERIALIZATION_HPP

#include <concepts>
#include <nlohmann/json.hpp>
#include <filesystem>
#include <string>
#include <fstream>
#include <Iterators.hpp>
#include <optional>
#include <vector>
#include <sstream>

#include "util/traits.hpp"

#define DELIM ,

#ifndef __JSON_INDENT__
// Should be defined by cmake
#define __JSON_INDENT__
#endif

#define DEFINE_SERIALIZATION(TTYPE, Type, ...)                                              \
template<TTYPE>                                                                             \
void to_json(nlohmann::json& nlohmann_json_j, const Type& nlohmann_json_t) {                \
    NLOHMANN_JSON_EXPAND(NLOHMANN_JSON_PASTE(NLOHMANN_JSON_TO, __VA_ARGS__))                \
}                                                                                           \
template<TTYPE>                                                                             \
void from_json(const nlohmann::json& nlohmann_json_j, Type& nlohmann_json_t) {              \
    Type nlohmann_json_default_obj{};                                                         \
    NLOHMANN_JSON_EXPAND(NLOHMANN_JSON_PASTE(NLOHMANN_JSON_FROM_WITH_DEFAULT, __VA_ARGS__)) \
}

namespace nlohmann {

    template<typename T>
    struct adl_serializer<std::optional<T>> {
        static void to_json(json &j, const std::optional<T> &optional) {
            if (optional.has_value()) {
                j = *optional;
            } else {
                j = nullptr;
            }
        }

        static void from_json(const json &j, std::optional<T> &optional) {
            if (j.is_null()) {
                optional = std::nullopt;
            } else {
                optional = j.get<T>();
            }
        }
    };
}

namespace tempo {
    template<typename T>
    class Scheduler;
}

/**
 * @brief namespace containing serializable datastructures that can be used to save problems and solutions as json
 */
namespace tempo::serialization {

    /**
     * A type that can be serialized by the nlohmann json library
     * @tparam T
     */
    template<typename T>
    concept serializable = requires(const T &instance, nlohmann::json &j) {
        nlohmann::to_json(j, instance);
    };

    using Branch = std::vector<std::pair<var_t, bool>>;

    /**
     * @brief Represents a serializable solution of a scheduling problem
     * @details @copybrief
     * @tparam T timing type
     */
    template<concepts::scalar T>
    struct Solution {
        Solution() = default;

        Solution(unsigned int id, T objective, Branch decisions) : id(id), objective(objective),
                                                                   decisions(std::move(decisions)) {}

        unsigned id;
        T objective;
        Branch decisions;
    };

    /**
     * @brief Represents a serializable intermediate state or partial problem state of a scheduling problem.
     * @details @copybrief
     */
    struct PartialProblem {
        PartialProblem() = default;

        PartialProblem(unsigned int associatedSolution, Branch decisions)
                : associatedSolution(associatedSolution), decisions(std::move(decisions)) {}

        unsigned associatedSolution;
        Branch decisions;
    };

    DEFINE_SERIALIZATION(concepts::scalar T, Solution<T>, id, objective, decisions)
    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(PartialProblem, associatedSolution, decisions)

    /**
     * Converts a serializable object to json and writes it to the specified destination
     * @tparam S Type of object to be serialized
     * @param object Object to be serialized
     * @param destination Destination file
     * @throws std::runtime_error if file could not be opened for writing
     */
    template<serializable S>
    void serializeToFile(const S &object, const std::filesystem::path &destination) {
        std::ofstream file(destination);
        if (not file.is_open()) {
            std::stringstream ss;
            ss << "unable to open file '" << destination << "' for writing";
            throw std::runtime_error(ss.str());
        }

        nlohmann::json j = object;
        file << j.dump(__JSON_INDENT__);
    }


    /**
     * Deserializes an object read from a file
     * @tparam T type of object
     * @param path path to the json file
     * @return deserialized object
     */
    template<typename T>
    T deserializeFromFile(const std::filesystem::path &path) {
        std::ifstream file(path);
        if (not file.is_open()) {
            std::stringstream ss;
            ss << "unable to open file '" << path << "'";
            throw std::runtime_error(ss.str());
        }

        return nlohmann::json::parse(file).get<T>();
    }

} // tempo

#endif //TEMPO_SERIALIZATION_HPP
