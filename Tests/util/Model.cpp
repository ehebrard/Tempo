/**
* @author Tim Luchterhand
* @date 22.10.24
* @brief
*/

#include <gtest/gtest.h>

#include "Solver.hpp"
#include "Literal.hpp"
#include "Model.hpp"
#include "util/Matrix.hpp"

struct TestResource: std::vector<tempo::Interval<int>> {
    tempo::Matrix<tempo::Literal<int>> literals;

    explicit TestResource(decltype(literals) literals) noexcept : literals(std::move(literals)) {}

    [[nodiscard]] const auto &getDisjunctiveLiterals() const noexcept {
        return literals;
    }
};

TEST(model, booleanVarsFromResources) {
    using namespace tempo;
    Matrix<Literal<int>> literals(2, 2, Contradiction<int>);
    literals.at(0, 0) = makeBooleanLiteral<int>(true, 15);
    literals.at(0, 1) = makeBooleanLiteral<int>(false, 1, 17);
    literals.at(1, 0) = makeBooleanLiteral<int>(true, 1, 17);
    TestResource resource(std::move(literals));
    auto vars = booleanVarsFromResources(resource);
    ASSERT_EQ(vars.size(), 2);
    std::ranges::sort(vars, {}, [](const auto &var) { return var.id(); });
    EXPECT_EQ(vars.front(), BooleanVar(1, 17));
    EXPECT_EQ(vars.back(), BooleanVar(15, 0));
    literals.fill(3, 3, Contradiction<int>);
    literals.at(1, 2) = makeBooleanLiteral<int>(false, 15, 0);
    literals.at(2, 2) = makeBooleanLiteral<int>(false, 3, 14);
    std::vector resources{std::move(resource), TestResource(std::move(literals))};
    vars = booleanVarsFromResources(resources);
    ASSERT_EQ(vars.size(), 3);
    std::ranges::sort(vars, {}, [](const auto &var) { return var.id(); });
    EXPECT_EQ(vars.front(), BooleanVar(1, 17));
    EXPECT_EQ(vars.at(1), BooleanVar(3, 14));
    EXPECT_EQ(vars.back(), BooleanVar(15, 0));
}