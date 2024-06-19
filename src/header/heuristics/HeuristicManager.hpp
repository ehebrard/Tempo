//
// Created by tim on 16.11.22.
//

#ifndef _TEMPO_HEURISTICMANAGER_HPP
#define _TEMPO_HEURISTICMANAGER_HPP

#include <variant>
#include <optional>
#include <exception>
#include "VSIDS.hpp"
#include "Tightest.hpp"
#include "WeightedDegree.hpp"
//#include "WeightedCriticalPath.hpp"
//#include "EpsilonGreedyBase.hpp"
//#include "util/traits.hpp"
//#include "util/random.hpp"



namespace tempo {
    template<typename T>
    class Scheduler;

template<typename T>
class Solver;
}

/**
 * @brief namespace containing variable (choice point) selection heuristics
 */
namespace tempo::heuristics {
    /**
     * @brief Heuristic factory class that can be used to construct different heuristics and at the same time provides a
     * consistent interface to callers
     */
     template<typename T>
    class HeuristicManager {
       //       using Distance = const RestrictedDistanceMatrix<T> &;
       using Distance = const Scheduler<T> &;
       using Implementations =
           std::variant<Tightest<T>, VSIDS<T>, WeightedDegree<T>
                        //        , EpsilonGreedyVSIDS
                        >;

     public:
       /**
        * Ctor: Internally constructs the heuristic inferred from the given
        * arguments
        * @tparam T type of scheduler
        * @param scheduler scheduler for which to create a heuristic
        * @param options options specifying the type of heuristic and further
        * config values
        * @throws std::runtime_error if an unknown heuristics type was given in
        * options
        */
   
        
        HeuristicManager(Solver<T> &solver, const Options &options) {
          switch (options.choice_point_heuristics) {
              case Options::ChoicePointHeuristics::Tightest: {
                  impl.emplace(std::in_place_type<Tightest<T>>);
                  break;
              }
              case Options::ChoicePointHeuristics::VSIDS: {
                  if (options.learning) {
                      impl.emplace(std::in_place_type<VSIDS<T>>, solver);
                  } else // closest thing if not learning
                      impl.emplace(std::in_place_type<WeightedDegree<T>>, solver);
                  break;
              }
              case Options::ChoicePointHeuristics::WeightedDegree: {
                  impl.emplace(std::in_place_type<WeightedDegree<T>>, solver);
                  break;
              }
//          case Options::ChoicePointHeuristics::WeightedCriticalPath:
//            impl.emplace(std::in_place_type<WeightedDegree<T>>, solver, true);
//            break;
            //                case Options::ChoicePointHeuristics::EG_VSIDS:
            //                  impl.emplace(std::in_place_type<EpsilonGreedyVSIDS>,
            //                               options.vsids_epsilon, scheduler,
            //                               options, scheduler);
            //                  break;
            //
            //#if __TORCH_ENABLED__
            //                case Options::ChoicePointHeuristics::GNN_HeatMap:
            //                    impl.emplace(std::in_place_type<HeatMap>,
            //                    options.gnn_model_location,
            //                    options.feature_extractor_conf,
            //                                 scheduler);
            //                    break;
            //#endif
          default:
            throw std::runtime_error("unknown heuristic type");
          }
        }

        /**
         * Calls the internally stored heuristic with the given arguments
         * @tparam T type of scheduler
         * @param scheduler scheduler for which to select the next choice point
         * @return choice point selected by the internal heuristic
         */
        
        auto nextChoicePoint(Solver<T> &solver) {
            return std::visit([&solver](auto &heuristic) { return heuristic.nextChoicePoint(solver); }, *impl);
        }

    private:
        std::optional<Implementations> impl;
    };
}

#endif //TEMPO_HEURISTICMANAGER_HPP

