/**
* @author Tim Luchterhand
* @date 13.06.24
* @brief
*/

#include <gtest/gtest.h>

#include "Literal.hpp"


TEST(util, LiteralStorage_Ctor) {
    using namespace tempo::detail;
    LiteralStorage ls(17, 0.25, Numeric{});
    EXPECT_TRUE(ls.isNumeric());
    EXPECT_TRUE(ls.sign());
    EXPECT_EQ(ls.id(), 17);
    EXPECT_EQ(ls.value(), 0.25f);
    ls = {1168623114, 14, Boolean{}};
    EXPECT_FALSE(ls.isNumeric());
    EXPECT_FALSE(ls.sign());
    EXPECT_EQ(ls.id(), 1168623114);
    EXPECT_EQ(ls.semantic(), 14);
    LiteralStorage<std::uint16_t> smallLit(29, 14, Numeric{});
    EXPECT_TRUE(smallLit.sign());
    EXPECT_EQ(smallLit.id(), 29);
    EXPECT_EQ(smallLit.value(), 14);

}

TEST(util, LiteralStorage_memory) {
    using namespace tempo::detail;
    EXPECT_EQ(sizeof(LiteralStorage<float>), 8);
    EXPECT_EQ(alignof(LiteralStorage<float>), 4);
}

TEST(util, LiteralStorage_access) {
    using namespace tempo::detail;
    LiteralStorage ls(13, 1.4f, Numeric{});
    EXPECT_TRUE(ls.isNumeric());
    EXPECT_EQ(ls.value(), 1.4f);
    ls.setValue(0.75f);
    EXPECT_EQ(ls.id(), 13);
    EXPECT_EQ(ls.value(), 0.75);
    ls = {18, 14, Boolean{}};
    EXPECT_FALSE(ls.isNumeric());
    EXPECT_EQ(ls.semantic(), 14);
    ls.setValue(0.34f);
    EXPECT_TRUE(ls.isNumeric());
    EXPECT_EQ(ls.value(), 0.34f);
}

TEST(util, LiteralStorage_invalid_access) {
    using namespace tempo::detail;
    LiteralStorage ls(13, 0.17f, Numeric{});
    EXPECT_THROW((void)ls.semantic(), std::bad_variant_access);
    ls = {13, 14, Boolean{}};
    EXPECT_THROW((void)ls.value(), std::bad_variant_access);
}

TEST(util, LiteralStorage_creator_functions) {
    using namespace tempo::detail;
    auto ls = makeSemanticLit<float>(17, 4);
    EXPECT_TRUE(ls.sign());
    EXPECT_FALSE(ls.isNumeric());
    EXPECT_EQ(ls.semantic(), 4);
    EXPECT_EQ(ls.id(), 17);
    ls = makeNumericLit(12, 0.3f);
    EXPECT_FALSE(ls.sign());
    EXPECT_TRUE(ls.isNumeric());
    EXPECT_EQ(ls.value(), 0.3f);
    EXPECT_EQ(ls.id(), 12);
}

template<tempo::concepts::scalar T>
void testSerializationRt(tempo::Literal<T> lit) {
    nlohmann::json j;
    j = lit;
    auto rt = j.get<tempo::Literal<T>>();
    EXPECT_EQ(rt, lit);
}

TEST(util, LiteralStorage_serialization) {
    using namespace tempo;
    auto lss = makeBooleanLiteral<float>(true, 4, 17);
    auto lsn = makeNumericLiteral(false, 7, 31);
    auto lsn1 = makeNumericLiteral(true, 19, 0.4f);
    testSerializationRt(lss);
    testSerializationRt(lsn);
    testSerializationRt(lsn1);
}