/**
* @author Tim Luchterhand
* @date 02.08.24
* @brief
*/

#include <chrono>
#include <sstream>
#include <stdexcept>
#include <array>
#include <cstdio>

#include "shell.hpp"

namespace shell {

    std::string getTimeStamp() {
        using namespace std::chrono;
        auto now = system_clock::now();
        year_month_day date(floor<days>(now));
        std::stringstream ss;
        // this needs to be this ugly because the std implementation on older g++ versions is broken
        auto d = static_cast<unsigned>(date.day());
        auto m = static_cast<unsigned>(date.month());
        ss << static_cast<int>(date.year()) << "-" << (m < 10 ? "0" : "") << m << "-"
           << (d < 10 ? "0" : "") << d;
        return ss.str();
    }

    // stolen from https://dev.to/aggsol/calling-shell-commands-from-c-8ej
    auto execCommand(const std::string &cmd) -> std::pair<int, std::string> {
        int exitStatus = 0;
        auto pPipe = ::popen(cmd.c_str(), "r");
        if (pPipe == nullptr) {
            throw std::runtime_error("Cannot open pipe");
        }

        std::array<char, 256> buffer{};
        std::string result;
        while (not std::feof(pPipe)) {
            auto bytes = std::fread(buffer.data(), 1, buffer.size(), pPipe);
            result.append(buffer.data(), bytes);
        }

        auto rc = ::pclose(pPipe);
        if (WIFEXITED(rc)) {
            exitStatus = WEXITSTATUS(rc);
        }

        return {exitStatus, std::move(result)};
    }

    auto getCommit() -> std::optional<std::string> {
        auto [status, res] = execCommand("git rev-parse HEAD");
        if (status != 0) {
            return {};
        }

        if (res.ends_with("\n")) {
            res.pop_back();
        }

        return res;
    }
}
