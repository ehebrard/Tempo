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
     * get the current day
     * @return day in yyyy-mm-dd format
     */
    std::string getTimeStamp();
}



#endif //TEMPO_SHELL_HPP
