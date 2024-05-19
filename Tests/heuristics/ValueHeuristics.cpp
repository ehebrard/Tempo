/**
 * @author Tim Luchterhand
 * @date 13.05.24
 * @brief
 */

#include <gtest/gtest.h>
#include <vector>

#include "Scheduler.hpp"
#include "heuristics/SolutionGuided.hpp"
#include "heuristics/TightestValue.hpp"
#include "util/traits.hpp"

struct TestValueHeuristic
    : public tempo::heuristics::BaseValueHeuristic<TestValueHeuristic> {
  bool called = false;

  explicit TestValueHeuristic(double epsilon)
      : tempo::heuristics::BaseValueHeuristic<TestValueHeuristic>(epsilon) {}

  template <typename Sched> auto choose(tempo::var x, const Sched &) {
    called = true;
    return tempo::POS(x);
  }
};

TEST(value_heuristics, base_value_heuristic) {
  using namespace tempo;
  TestValueHeuristic h(0);
  EXPECT_EQ(h.choosePolarity(0, 0), POS(0));
  EXPECT_TRUE(h.called);
  h = TestValueHeuristic(1);
  h.choosePolarity(0, 0);
  EXPECT_FALSE(h.called);
  EXPECT_THROW(TestValueHeuristic(2), std::runtime_error);
}

struct DummyScheduler {
  bool invert = false;
  std::vector<bool> solution{};
  bool hasSol = false;

  bool hasSolution() const { return hasSol; }

  const auto &getSolution() const { return solution; }

  auto getEdge(tempo::lit lit) const -> tempo::DistanceConstraint<int> {
    return lit % 2 == 0 and not invert ? tempo::DistanceConstraint{10, 20, 6}
                                       : tempo::DistanceConstraint{11, 21, 2};
  }

  int upper(tempo::event e) const { return e % 2 == 0 ? 20 : 10; }

  int lower(tempo::event e) const { return e % 2 == 0 ? 0 : 8; }
};

TEST(value_heuristics, TightestValue) {
  using namespace tempo;
  using namespace tempo::heuristics;
  EXPECT_TRUE((value_heuristic<TightestValue, Scheduler<int>>));
  DummyScheduler sched{.invert = false};
  EXPECT_EQ(TightestValue::choose(1, sched), NEG(1));
  sched.invert = true;
  EXPECT_EQ(TightestValue::choose(1, sched), POS(1));
}

TEST(value_heuristics, SolutionGuided) {
  using namespace tempo;
  using namespace tempo::heuristics;
  EXPECT_TRUE((value_heuristic<SolutionGuided, Scheduler<int>>));
  DummyScheduler sched;
  SolutionGuided h(0);
  EXPECT_EQ(h.choose(1, sched), TightestValue::choose(1, sched));
  sched.invert = true;
  EXPECT_EQ(h.choose(1, sched), TightestValue::choose(1, sched));
  sched.solution = {true, false, false};
  sched.hasSol = true;
  EXPECT_EQ(h.choose(0, sched), POS(0));
  EXPECT_EQ(h.choose(1, sched), NEG(1));
  EXPECT_EQ(h.choose(2, sched), NEG(2));
}