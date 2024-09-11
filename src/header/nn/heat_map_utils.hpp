/**
* @author Tim Luchterhand
* @date 11.09.24
* @brief contains functions that calculate different probabilites from edge heat maps
*/

#ifndef TEMPO_HEAT_MAP_UTILS_HPP
#define TEMPO_HEAT_MAP_UTILS_HPP

#include "DistanceConstraint.hpp"
#include "util/Matrix.hpp"
#include "util/SchedulingProblemHelper.hpp"
#include "torch_types.hpp"
#include "util/traits.hpp"
#include "Literal.hpp"
#include "GNN.hpp"

namespace tempo {
    template<typename T>
    class Solver;
}

namespace tempo::nn {
    /**
     * Transforms two evidence masses into a Bayes probability according to https://arxiv.org/pdf/1605.02406.pdf
     * @note If m1 and m2 are Bernoulli probabilities, the result is m1.
     * @tparam F floating point type
     * @param m1 evidence mass for event
     * @param m2 counter evidence mass for event
     * @return Bayes probability of the event
     */
    template<std::floating_point F>
    F dstToBayes(F m1, F m2) {
        assert (m1 >= 0 and m2 >= 0 and m1 <= 1 and m2 <= 1);
        return m1 + F(0.5) * (1 - m1 - m2);
    }

    /**
     * Extracts a probability mass for an edge between two tasks
     * @param taskFrom source node of the edge
     * @param taskTo destination node of the edge
     * @param heatMap heat map containing probability masses
     * @return value in the heat map corresponding to the edge
     * @throws std::runtime_error when heat map does not contain information on the edge
     */
    DataType probabilityMass(unsigned taskFrom, unsigned taskTo, const Matrix<DataType> &heatMap);


    /**
     * Overload of probabilityMass for edges between variables
     * @tparam T timing type
     * @param edge edge between two variables (not tasks)
     * @param heatMap heat map containing probability masses
     * @param mapping mapping from variables to tasks
     * @return value in the heat map corresponding to the edge
     */
    template<concepts::scalar T>
    DataType probabilityMass(const DistanceConstraint<T> &edge, const Matrix<DataType> &heatMap,
                             const VarTaskMapping &mapping) {
        assert(mapping.contains(edge.from));
        assert(mapping.contains(edge.to));
        unsigned from = mapping(edge.from);
        unsigned to = mapping(edge.to);
        return probabilityMass(from, to, heatMap);
    }


    /**
     * Overload of probabilityMass. Uses edge corresponding to given variable
     * @note Calling this function on a numeric variable is undefined behavior!
     * @tparam T timing type
     * @param polarity polarity of the binary variable
     * @param x variable identifier
     * @param heatMap heat map containing probability masses
     * @param solver solver that holds all variables
     * @param mapping mapping from variables to tasks
     * @return value in the heat map corresponding to the edge
     */
    template<concepts::scalar T>
    double probabilityMass(bool polarity, var_t x, const Matrix<DataType> &heatMap, const Solver<T> &solver,
                           const VarTaskMapping &mapping) {
        return probabilityMass(solver.boolean.getEdge(polarity, x), heatMap, mapping);
    }

    /**
     * Overload of probabilityMass. Uses edge corresponding to given literal
     * @tparam T timing type
     * @param lit literal corresponding to an edge
     * @param heatMap heat map containing probability masses
     * @param solver solver that holds all variables
     * @param mapping mapping from variables to tasks
     * @return value in the heat map corresponding to the edge
     */
    template<concepts::scalar T>
    double probabilityMass(Literal<T> lit, const Matrix<DataType> &heatMap, const Solver<T> &solver,
                           const VarTaskMapping &mapping) {
        return probabilityMass(solver.boolean.getEdge(lit), heatMap, mapping);
    }

    /**
     * Calculates a Bernoulli probability that indicates which of two task edges corresponding to the variable x is
     * more likely.
     * @note Calling this function on a numeric variable is undefined behavior!
     * @tparam T timing type
     * @param x variable index
     * @param heatMap heat map containing probability masses
     * @param solver solver that holds all variables
     * @param mapping mapping from variables to tasks
     * @return probability of edge corresponding to positive literal of x being present
     */
    template<concepts::scalar T>
    double pignisticEdgeProbability(var_t x, const Matrix<DataType> &heatMap, const Solver<T> &solver,
                                    const VarTaskMapping &mapping) {
        auto m1 = probabilityMass(true, x, heatMap, solver, mapping);
        auto m2 = probabilityMass(false, x, heatMap, solver, mapping);
        return dstToBayes(m1, m2);
    }

    /**
     * Overload of pignisticEdgeProbability. Extracts edges corresponding to a boolean literal from the solver.
     * @tparam T timing type
     * @param lit literal
     * @param heatMap heat map containing probability masses
     * @param solver solver that holds all variables
     * @param mapping mapping from variables to tasks
     * @return probability of edge corresponding to literal l being present
     */
    template<concepts::scalar T>
    double pignisticEdgeProbability(Literal<T> lit, const Matrix<DataType> &heatMap, const Solver<T> &solver,
                                    const VarTaskMapping &mapping) {
        if (not lit.isBoolean() or not lit.hasSemantic()) {
            throw std::runtime_error("literal needs to be boolean and have a semantic");
        }

        auto m1 = probabilityMass(lit, heatMap, mapping);
        auto m2 = probabilityMass(~lit, heatMap, mapping);
        return dstToBayes(m1, m2);
    }
}

#endif //TEMPO_HEAT_MAP_UTILS_HPP
