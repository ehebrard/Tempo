/**
* @author Tim Luchterhand
* @date 02.08.24
* @brief
*/

#ifndef TEMPO_SHELL_HPP
#define TEMPO_SHELL_HPP

#include <string>
#include <utility>
#include <optional>

namespace shell {
    /**
     * runs a command in the shell
     * @param cmd command to be executed
     * @return pair(exit code, command output)
     * @note command output only supports small outputs
     */
    auto execCommand(const std::string &cmd) -> std::pair<int, std::string>;

    /**
     * Gets the current commit hash
     * @return commit hash or nothing if command failed to run
     */
    auto getCommit() -> std::optional<std::string>;

    /**
     * get the current day
     * @return day in yyyy-mm-dd format
     */
    std::string getTimeStamp();
}



#endif //TEMPO_SHELL_HPP
