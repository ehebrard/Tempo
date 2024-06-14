/**
* @author Tim Luchterhand
* @date 13.06.24
* @brief
*/

#include <gtest/gtest.h>

#include "util/LiteralStorage.hpp"


TEST(util, LiteralStorage_Ctor) {
    using namespace tempo;
    LiteralStorage ls(true, 17, NumericValue(0.25f));
    EXPECT_TRUE(ls.isNumeric());
    EXPECT_TRUE(ls.sign());
    EXPECT_EQ(ls.id(), 17);
    EXPECT_EQ(ls.value(), 0.25f);
    ls = {false, 169, 14};
    EXPECT_FALSE(ls.isNumeric());
    EXPECT_FALSE(ls.sign());
    EXPECT_EQ(ls.id(), 169);
    EXPECT_EQ(ls.semantic(), 14);
    LiteralStorage<std::uint16_t> smallLit(true, 29, NumericValue<std::uint16_t>(14));
    EXPECT_TRUE(smallLit.sign());
    EXPECT_EQ(smallLit.id(), 29);
    EXPECT_EQ(smallLit.value(), 14);

}

TEST(util, LiteralStorage_memory) {
    using namespace tempo;
    EXPECT_EQ(sizeof(LiteralStorage<float>), 8);
    EXPECT_EQ(alignof(LiteralStorage<float>), 4);
}

TEST(util, LiteralStorage_access) {
    using namespace tempo;
    LiteralStorage ls(true, 13, NumericValue(1.4f));
    EXPECT_TRUE(ls.isNumeric());
    EXPECT_EQ(ls.value(), 1.4f);
    ls.value() = 0.75;
    EXPECT_EQ(ls.id(), 13);
    EXPECT_EQ(ls.value(), 0.75);
    ls = {false, 18, 14};
    EXPECT_FALSE(ls.isNumeric());
    EXPECT_EQ(ls.semantic(), 14);
    ls.semantic() = 17;
    EXPECT_EQ(ls.semantic(), 17);
    EXPECT_EQ(ls.id(), 18);
}

TEST(util, LiteralStorage_invalid_access) {
    using namespace tempo;
    LiteralStorage ls(true, 13, NumericValue(0.17f));
    EXPECT_THROW(ls.semantic(), std::bad_variant_access);
    ls = {true, 13, 14};
    EXPECT_THROW(ls.value(), std::bad_variant_access);
}

TEST(util, LiteralStorage_creator_functions) {
    using namespace tempo;
    auto ls = makeSemanticLit<float>(true, 17, 4);
    EXPECT_TRUE(ls.sign());
    EXPECT_FALSE(ls.isNumeric());
    EXPECT_EQ(ls.semantic(), 4);
    EXPECT_EQ(ls.id(), 17);
    ls = makeNumericLit(false, 13, 0.3f);
    EXPECT_FALSE(ls.sign());
    EXPECT_TRUE(ls.isNumeric());
    EXPECT_EQ(ls.value(), 0.3f);
    EXPECT_EQ(ls.id(), 13);
}