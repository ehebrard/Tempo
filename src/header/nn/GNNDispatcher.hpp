/**
* @author Tim Luchterhand
* @date 14.03.25
* @file GNNDispatcher.hpp
* @brief
*/

#ifndef GNNDISPATCHER_HPP
#define GNNDISPATCHER_HPP

#include <memory>
#include <stdexcept>
#include <iostream>


#include "util/enum.hpp"
#include "util/traits.hpp"
#include "util/SubscribableEvent.hpp"
#include "util/Options.hpp"
#include "Solver.hpp"

namespace tempo::nn {
    PENUM(Dispatch, Full, SingleShot, OnSolution, OnRestart)

    struct BaseDispatcher {
        BaseDispatcher() = default;
        BaseDispatcher(const BaseDispatcher&) = default;
        BaseDispatcher(BaseDispatcher&&) = default;
        BaseDispatcher& operator=(const BaseDispatcher&) = default;
        BaseDispatcher& operator=(BaseDispatcher&&) = default;
        virtual ~BaseDispatcher() = default;

        virtual bool runInference() = 0;
    };

    struct Dispatcher : std::unique_ptr<BaseDispatcher> {
        using std::unique_ptr<BaseDispatcher>::unique_ptr;

        [[nodiscard]] bool runInference() const {
            return this->get()->runInference();
        }
    };

    struct FullGuidance final : BaseDispatcher {
        bool runInference() override;
    };

    struct SingleShotDispatcher : BaseDispatcher {
    protected:
        bool inferenceAllowed = true;
    public:
        bool runInference() override;
    };

    template<bool SolutionOnly>
    class RestartDispatcher final : public SingleShotDispatcher {
        SubscriberHandle restartHandler;
    public:
        RestartDispatcher(const RestartDispatcher&) = delete;
        RestartDispatcher(RestartDispatcher&&) = delete;
        RestartDispatcher& operator=(const RestartDispatcher&) = delete;
        RestartDispatcher& operator=(RestartDispatcher&&) = delete;
        ~RestartDispatcher() override = default;

        template<concepts::scalar T>
        explicit RestartDispatcher(const Solver<T> &solver): restartHandler(
            solver.SearchRestarted.subscribe_handled([this](bool onSolution) {
                this->inferenceAllowed = not SolutionOnly or onSolution;
            })) {}
    };

    template<concepts::scalar T>
    auto make_dispatcher(Dispatch dispatch, const Solver<T> &solver) -> Dispatcher {
        if (solver.getOptions().verbosity >= Options::QUIET) {
            std::cout << "-- Using GNN dispatcher '" << dispatch << "'" << std::endl;
        }

        switch (dispatch) {
            case Dispatch::Full:
                return std::make_unique<FullGuidance>();
            case Dispatch::SingleShot:
                return std::make_unique<SingleShotDispatcher>();
            case Dispatch::OnSolution:
                return std::make_unique<RestartDispatcher<true>>(solver);
            case Dispatch::OnRestart:
                return std::make_unique<RestartDispatcher<false>>(solver);
            default:
                throw std::invalid_argument("unsupported dispatch type " + penum_to_string(dispatch));
        }
    }

}

#endif //GNNDISPATCHER_HPP
