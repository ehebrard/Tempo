/**
* @author Tim Luchterhand
* @date 12.11.24
* @brief
*/

#ifndef TEMPO_SOLUTION_HPP
#define TEMPO_SOLUTION_HPP

#include <vector>
#include <ostream>
#include <istream>
#include <Iterators.hpp>

#include "Global.hpp"
#include "Model.hpp"
#include "Solver.hpp"
#include "util/traits.hpp"

namespace tempo {

    template<concepts::scalar T>
    class NumericSolution {
        std::vector<T> lb;
        std::vector<T> ub;
    public:

        NumericSolution() = default;

        NumericSolution(std::vector<T> lower, std::vector<T> upper) noexcept: lb(std::move(lower)),
                                                                              ub(std::move(upper)) {}

        explicit NumericSolution(const Solver<T> &solver) noexcept: NumericSolution(
                solver.numeric.bestSolution(bound::lower),
                solver.numeric.bestSolution(bound::upper)) {}

        [[nodiscard]] std::size_t size() const noexcept {
            return std::min(ub.size(), lb.size());
        }

        T lower(const NumericVar<T> x) const noexcept {
            return -lb[x.id()];
        }

        T lower(const var_t x) const noexcept {
            return -lb[x];
        }

        T upper(const NumericVar<T> x) const noexcept {
            return ub[x.id()];
        }

        T upper(const var_t x) const noexcept {
            return ub[x];
        }

        auto begin() const noexcept {
            return iterators::zip_i(lb.begin(), ub.begin());
        }

        auto end() const noexcept {
            return iterators::zip_i(lb.end(), ub.end());
        }
    };

    class BooleanSolution {
        std::vector<bool> polarities;
    public:
        BooleanSolution() = default;

        explicit BooleanSolution(std::vector<bool> polarities) noexcept;

        [[nodiscard]] std::size_t size() const noexcept;

        [[nodiscard]] auto begin() const noexcept -> std::vector<bool>::const_iterator;

        [[nodiscard]] auto end() const noexcept -> std::vector<bool>::const_iterator;

        template<concepts::scalar T>
        explicit BooleanSolution(const Solver<T> &solver): polarities(solver.boolean.size()) {
            var_t end_bool{static_cast<var_t>(solver.boolean.size())};
            for (var_t i{0}; i < end_bool; ++i) {
                polarities[i] = solver.boolean.isTrue(i);
            }
        }

        template<concepts::scalar T>
        bool value(const BooleanVar<T> x) const noexcept {
            return polarities[x.id()];
        }

        [[nodiscard]] bool value(var_t x) const noexcept;

        template<concepts::scalar T>
        bool consistent(Literal<T> lit) const noexcept {
            return value(lit.variable()) == lit.sign();
        }
    };

    template<concepts::scalar T>
    struct Solution {
        NumericSolution<T> numeric;
        BooleanSolution boolean;

        Solution() = default;

        explicit Solution(const Solver<T> &solver) : numeric(solver), boolean(solver) {}

        bool reachable(const Solver<T> &solver) const {
            return discrepancy(solver) == 0;
        }

        var_t discrepancy(const Solver<T> &solver) const {
            bool canreach{true};
            var_t end_bool{static_cast<var_t>(boolean.size())};
            var_t end_num{static_cast<var_t>(numeric.size())};
            var_t i{1};
            for (; canreach and i < end_bool; ++i) {
                if (boolean.value(i)) {
                    canreach = not solver.boolean.isFalse(i);
                } else {
                    canreach = not solver.boolean.isTrue(i);
                }
            }

            if (i < end_bool) {
                return i;
            }

            i = 1;
            for (; canreach and i < end_num; ++i) {
                canreach = (solver.numeric.lower(i) <= numeric.lower_bound(i) and
                            solver.numeric.upper(i) >= numeric.upper_bound(i));
            }

            if (i < end_num) {
                return -i;
            }

            return 0;
        }

        std::ostream &display(std::ostream &os) const {
            os << boolean.size() << " " << numeric.size();
            for (auto b: boolean) {
                os << " " << b;
            }

            for (auto [lb, ub]: numeric) {
                os << " [" << lb << ", " << ub << "]";
            }

            return os;
        }

        std::istream &load(std::istream &is) {
            std::vector<bool> polarities;
            std::vector<T> ub, lb;
            size_t nb;
            size_t nn;
            bool v;
            is >> nb;
            is >> nn;
            polarities.resize(nb);
            lb.resize(nn);
            ub.resize(nn);
            for (size_t i{0}; i < nb; ++i) {
                is >> v;
                polarities[i] = v;
            }
            for (size_t i{0}; i < nn; ++i) {
                is >> lb[i];
            }
            for (size_t i{0}; i < nn; ++i) {
                is >> ub[i];
            }

            boolean = BooleanSolution(std::move(polarities));
            numeric = NumericSolution(std::move(lb), std::move(ub));
            return is;
        }
    };

    template <typename T>
    std::ostream &operator<<(std::ostream &os, const Solution<T> &x) {
        return x.display(os);
    }

    template <typename T>
    std::istream &operator>>(std::istream &is, Solution<T> &x) {
        return x.load(is);
    }
}

#endif //TEMPO_SOLUTION_HPP
