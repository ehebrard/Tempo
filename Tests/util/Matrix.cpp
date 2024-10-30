/**
 * @author Tim Luchterhand
 * @date 21.03.23.
 */

#include <gtest/gtest.h>
#include <ranges>


#include "util/Matrix.hpp"

TEST(util, Matrix_default_ctor) {
    using namespace tempo;
    Matrix<int> m;
    EXPECT_EQ(m.numRows(), 0);
    EXPECT_EQ(m.numColumns(), 0);
    EXPECT_EQ(m.storageLayout(), Layout::RowMajor);
}

TEST(util, Matrix_ctor) {
    using namespace tempo;
    Matrix<int> m(2, 3);
    EXPECT_EQ(m.numRows(), 2);
    EXPECT_EQ(m.numColumns(), 3);
    m.for_each([](const auto &val) { EXPECT_EQ(val, 0); });
}

TEST(util, Matrix_initial_value) {
    using namespace tempo;
    Matrix<int> m(2, 3, 4);
    m.for_each([](const auto &val) { EXPECT_EQ(val, 4); });
}

TEST(util, Matrix_Ctor_range) {
    using namespace tempo;
    Matrix<int> mFromInitList(2, 3, {1, 2, 3, 4, 5, 6});
    const std::array vals{1, 2, 3, 4, 5, 6};
    Matrix<int> mFromRange(2, 3, vals);
    EXPECT_EQ(mFromInitList.rawData(), mFromRange.rawData());
    EXPECT_TRUE(std::ranges::equal(mFromRange.rawData(), vals));
}

TEST(util, Matrix_Ctor_range_exception) {
    using namespace tempo;
    EXPECT_THROW((Matrix<int>(2, 3, {1, 2, 3, 4, 5})), std::runtime_error);
    const std::array vals{1, 2, 3, 4, 5};
    EXPECT_THROW((Matrix<int>(2, 3, vals)), std::runtime_error);
}

TEST(util, Matrix_Ctor_functor) {
    using namespace tempo;
    Matrix<int> source(2, 3, {1, 2, 3, 4, 5, 6});
    Matrix<int> destination(2, 2, source);
    EXPECT_EQ(destination, Matrix<int>(2, 2, {1, 2, 4, 5}));
    EXPECT_EQ(destination.storageLayout(), Layout::RowMajor);
    destination = {2, 2, source, Layout::ColMajor};
    EXPECT_EQ(destination, Matrix<int>(2, 2, {1, 2, 4, 5}));
    EXPECT_EQ(destination.storageLayout(), Layout::ColMajor);
}

TEST(util, Matrix_move) {
    using namespace tempo;
    Matrix<int> matrix(2, 3, {1, 2, 3, 4, 5, 6});
    auto copy = matrix;
    copy.at(1, 2) = 17;
    EXPECT_NE(matrix, copy);
    const auto * const dataPtr = matrix.rawData().data();
    auto moved = std::move(matrix);
    EXPECT_EQ(dataPtr, moved.rawData().data());
    EXPECT_EQ(matrix.numColumns(), 0);
    EXPECT_EQ(matrix.numRows(), 0);
    matrix.resize(3, 3);
    matrix = std::move(moved);
    EXPECT_EQ(dataPtr, matrix.rawData().data());
    EXPECT_EQ(moved.numColumns(), 0);
    EXPECT_EQ(moved.numRows(), 0);
    EXPECT_EQ(matrix.at(1, 2), 6);
    matrix = copy;
    EXPECT_EQ(matrix.at(1, 2), 17);
    matrix.at(1, 2) = 0;
    EXPECT_EQ(copy.at(1, 2), 17);
}

TEST(util, Matrix_layout) {
    using namespace tempo;
    Matrix<int> m(2, 3, 4, Layout::ColMajor);
    EXPECT_EQ(m.storageLayout(), Layout::ColMajor);
}

TEST(util, Matrix_value_access) {
    using namespace tempo;
    Matrix<int> m(2, 3, 4);
    EXPECT_EQ(m(1, 2), 4);
    m(0, 1) = -2;
    EXPECT_EQ(m(0, 1), -2);
    EXPECT_EQ(m(0, 1), m.at(0, 1));
    EXPECT_THROW(m.at(2, 1), std::out_of_range);
}

TEST(util, Matrix_for_each) {
    using namespace tempo;
    Matrix<int> m(2, 2, 4, tempo::Layout::RowMajor);
    m.for_each([i = 0](auto &val) mutable { val = i++; });
    EXPECT_EQ(m(0, 0), 0);
    EXPECT_EQ(m(0, 1), 1);
    EXPECT_EQ(m(1, 0), 2);
    EXPECT_EQ(m(1, 1), 3);
}

TEST(util, Matrix_resize) {
    using namespace tempo;
    Matrix<int> m(2, 2, 4);
    m.resize(3, 4);
    EXPECT_EQ(m.numRows(), 3);
    EXPECT_EQ(m.numColumns(), 4);
    m.at(2, 3) = 17;
    EXPECT_EQ(m.at(2, 3), 17);
    m.resize(7, 10, 198);
    EXPECT_EQ(m.numRows(), 7);
    EXPECT_EQ(m.numColumns(), 10);
    EXPECT_EQ(m.at(6, 8), 198);
}

TEST(util, Matrix_fill) {
    using namespace tempo;
    Matrix<int> m(2, 2, 4);
    m.fill(19);
    EXPECT_EQ(m.numRows(), 2);
    EXPECT_EQ(m.numColumns(), 2);
    m.for_each([](auto val) { EXPECT_EQ(val, 19); });
    m.fill(1, 2, -3);
    EXPECT_EQ(m.numRows(), 1);
    EXPECT_EQ(m.numColumns(), 2);
    m.for_each([](auto val) { EXPECT_EQ(val, -3); });
    m.fill(3, 3, 1337);
    EXPECT_EQ(m.numRows(), 3);
    EXPECT_EQ(m.numColumns(), 3);
    m.for_each([](auto val) { EXPECT_EQ(val, 1337); });
}

TEST(util, Matrix_change_layout) {
    using namespace tempo;
    Matrix<int> m(2, 3, 4, tempo::Layout::RowMajor);
    m.for_each([i = 0](auto &val) mutable { val = i++; });
    EXPECT_EQ(m(0, 0), 0);
    EXPECT_EQ(m(0, 1), 1);
    EXPECT_EQ(m(0, 2), 2);
    EXPECT_EQ(m(1, 0), 3);
    EXPECT_EQ(m(1, 1), 4);
    EXPECT_EQ(m(1, 2), 5);
    m.changeLayout(tempo::Layout::ColMajor);
    EXPECT_EQ(m(0, 0), 0);
    EXPECT_EQ(m(0, 1), 2);
    EXPECT_EQ(m(0, 2), 4);
    EXPECT_EQ(m(1, 0), 1);
    EXPECT_EQ(m(1, 1), 3);
    EXPECT_EQ(m(1, 2), 5);
}

TEST(util, Matrix_equality) {
    using namespace tempo;
    Matrix<int> m1(2, 3, {1, 2, 3, 4, 5, 6});
    Matrix<int> m2(2, 3, {1, 2, 3, 4, 5, 6});
    EXPECT_EQ(m1, m2);
    m2.at(1, 1) = 17;
    EXPECT_NE(m1, m2);
    m2.at(1, 1) = 5;
    m2.resize(2, 2);
    EXPECT_NE(m1, m2);
    m2 = {2, 3, {1, 4, 2, 5, 3, 6}, Layout::ColMajor};
    EXPECT_EQ(m1, m2);
}

TEST(util, Matrix_row) {
    using namespace tempo;
    namespace r = std::ranges;
    Matrix<int> m(3, 3, {1, 2, 3, 4, 5, 6, 7, 8, 9});
    const auto &cm = m;
    auto first = m.row(0);
    auto last = cm.row(2);
    EXPECT_TRUE(r::equal(first, std::vector{1, 2, 3}));
    EXPECT_TRUE(r::equal(last, std::vector{7, 8, 9}));
    EXPECT_TRUE(std::is_const_v<std::remove_reference_t<r::range_reference_t<decltype(last)>>>);
    first[1] = 17;
    EXPECT_EQ(m.at(0, 1), 17);
}

TEST(util, Matrix_col) {
    using namespace tempo;
    namespace r = std::ranges;
    Matrix<int> m(3, 3, {1, 2, 3, 4, 5, 6, 7, 8, 9});
    const auto &cm = m;
    auto first = m.column(0);
    auto last = cm.column(2);
    EXPECT_TRUE(r::equal(first, std::vector{1, 4, 7}));
    EXPECT_TRUE(r::equal(last, std::vector{3, 6, 9}));
    EXPECT_TRUE(std::is_const_v<std::remove_reference_t<r::range_reference_t<decltype(last)>>>);
    first[1] = 17;
    EXPECT_EQ(m.at(1, 0), 17);
}
