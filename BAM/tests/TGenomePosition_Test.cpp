//
// Created by madleina on 02.08.20.
//

#include "../TGenomePosition.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

TEST(TGenomePositionTest, constructor_noParams){
    BAM::TGenomePosition pos;
    EXPECT_EQ(pos.refID(), 0);
    EXPECT_EQ(pos.position(), 0);
}

TEST(TGenomePositionTest, constructor_refID_pos){
    BAM::TGenomePosition pos(10, 100);
    EXPECT_EQ(pos.refID(), 10);
    EXPECT_EQ(pos.position(), 100);
}

TEST(TGenomePositionTest, constructor_otherObject){
    BAM::TGenomePosition pos1(10, 100);
    BAM::TGenomePosition pos2(pos1);
    EXPECT_EQ(pos2.refID(), 10);
    EXPECT_EQ(pos2.position(), 100);
}

TEST(TGenomePositionTest, move_refID_pos){
    BAM::TGenomePosition pos(10, 100);
    pos.move(1, 50);
    EXPECT_EQ(pos.refID(), 1);
    EXPECT_EQ(pos.position(), 50);
}

TEST(TGenomePositionTest, move_otherObject){
    BAM::TGenomePosition pos1(10, 100);
    BAM::TGenomePosition pos2(1, 50);
    pos1.move(pos2);
    EXPECT_EQ(pos1.refID(), 1);
    EXPECT_EQ(pos1.position(), 50);
    EXPECT_EQ(pos2.refID(), 1);
    EXPECT_EQ(pos2.position(), 50);
}

TEST(TGenomePositionTest, operator_Eq){
    BAM::TGenomePosition pos2(1, 50);
    BAM::TGenomePosition pos1(2, 100);
    pos1 = pos2;
    EXPECT_EQ(pos1.refID(), 1);
    EXPECT_EQ(pos1.position(), 50);
    EXPECT_EQ(pos2.refID(), 1);
    EXPECT_EQ(pos2.position(), 50);
}

TEST(TGenomePositionTest, operator_Plus){
    BAM::TGenomePosition pos1(10, 100);
    BAM::TGenomePosition pos2 = pos1 + 10;
    EXPECT_EQ(pos2.refID(), 10);
    EXPECT_EQ(pos2.position(), 110);
}

TEST(TGenomePositionTest, operator_Minus_length){
    BAM::TGenomePosition pos1(10, 100);
    // this is ok
    BAM::TGenomePosition pos2 = pos1 - 10;
    EXPECT_EQ(pos2.refID(), 10);
    EXPECT_EQ(pos2.position(), 90);
    // this is negative, so set it to zero
    pos2 = pos1 - 200;
    EXPECT_EQ(pos2.refID(), 10);
    EXPECT_EQ(pos2.position(), 0);
}

TEST(TGenomePositionTest, operator_Minus_otherObject){
    BAM::TGenomePosition pos1(10, 100);

    // same chromosome -> ok
    BAM::TGenomePosition pos2(10, 50);
    uint32_t position = pos1 - pos2;
    EXPECT_EQ(position, 50);

    // same chromosome, but pos2 is larger than pos1 -> set to zero
    pos2.move(10, 200);
    position = pos1 - pos2;
    EXPECT_EQ(position, 0);

    // other chromosomes -> throw
    pos2.move(1, 100);
    EXPECT_THROW(pos1 - pos2, std::runtime_error);
}

TEST(TGenomePositionTest, operator_PlusEq){
    BAM::TGenomePosition pos(10, 100);
    pos += 10;
    EXPECT_EQ(pos.refID(), 10);
    EXPECT_EQ(pos.position(), 110);
}

TEST(TGenomePositionTest, operator_MinusEq){
    BAM::TGenomePosition pos(10, 100);
    // this is ok
    pos -= 10;
    EXPECT_EQ(pos.refID(), 10);
    EXPECT_EQ(pos.position(), 90);

    // this is negative, so set it to zero
    pos -= 100;
    EXPECT_EQ(pos.refID(), 10);
    EXPECT_EQ(pos.position(), 0);
}

TEST(TGenomePositionTest, operator_PlusPlus){
    BAM::TGenomePosition pos(10, 100);
    ++pos;
    EXPECT_EQ(pos.refID(), 10);
    EXPECT_EQ(pos.position(), 101);
}

TEST(TGenomePositionTest, operator_MinusMinus){
    BAM::TGenomePosition pos(10, 100);
    // this is ok
    --pos;
    EXPECT_EQ(pos.refID(), 10);
    EXPECT_EQ(pos.position(), 99);

    // this is negative, so set it to zero
    pos.move(10, 0);
    --pos;
    EXPECT_EQ(pos.refID(), 10);
    EXPECT_EQ(pos.position(), 0);
}

TEST(TGenomePositionTest, sameChr){
    BAM::TGenomePosition pos1(10, 100);
    BAM::TGenomePosition pos2(10, 99);
    // same chr
    EXPECT_TRUE(pos1.sameChr(pos2));

    // other chr
    pos2.move(11, 99);
    EXPECT_FALSE(pos1.sameChr(pos2));
}

TEST(TGenomePositionTest, operator_equalEqual){
    BAM::TGenomePosition pos1(10, 100);
    BAM::TGenomePosition pos2(10, 100);
    // same pos and chr
    EXPECT_TRUE(pos1 == pos2);

    // same pos, other chr
    pos2.move(11, 100);
    EXPECT_FALSE(pos1 == pos2);

    // other pos, same chr
    pos2.move(10, 99);
    EXPECT_FALSE(pos1 == pos2);

    // other pos, other chr
    pos2.move(11, 99);
    EXPECT_FALSE(pos1 == pos2);
}

TEST(TGenomePositionTest, operator_smaller){
    BAM::TGenomePosition pos1(10, 100);
    BAM::TGenomePosition pos2(10, 200);
    // same chr, smaller pos
    EXPECT_TRUE(pos1 < pos2);

    // same chr, larger pos
    pos2.move(10, 1);
    EXPECT_FALSE(pos1 < pos2);

    // same chr, equal pos
    pos2.move(10, 100);
    EXPECT_FALSE(pos1 < pos2);

    // smaller chr
    pos2.move(11, 100);
    EXPECT_TRUE(pos1 < pos2);

    // larger chr
    pos2.move(9, 100);
    EXPECT_FALSE(pos1 < pos2);
}


TEST(TGenomePositionTest, operator_larger){
    BAM::TGenomePosition pos1(10, 200);
    BAM::TGenomePosition pos2(10, 100);
    // same chr, larger pos
    EXPECT_TRUE(pos1 > pos2);

    // same chr, smaller pos
    pos2.move(300, 1);
    EXPECT_FALSE(pos1 > pos2);

    // same chr, equal pos
    pos2.move(10, 200);
    EXPECT_FALSE(pos1 > pos2);

    // larger chr
    pos2.move(9, 100);
    EXPECT_TRUE(pos1 > pos2);

    // smaller chr
    pos2.move(11, 100);
    EXPECT_FALSE(pos1 > pos2);
}

TEST(TGenomePositionTest, operator_smallerEqual){
    BAM::TGenomePosition pos1(10, 100);
    BAM::TGenomePosition pos2(10, 200);
    // same chr, smaller pos
    EXPECT_TRUE(pos1 <= pos2);

    // same chr, larger pos
    pos2.move(10, 1);
    EXPECT_FALSE(pos1 <= pos2);

    // same chr, equal pos
    pos2.move(10, 100);
    EXPECT_TRUE(pos1 <= pos2);

    // smaller chr
    pos2.move(11, 100);
    EXPECT_TRUE(pos1 <= pos2);

    // larger chr
    pos2.move(9, 100);
    EXPECT_FALSE(pos1 <= pos2);
}

TEST(TGenomePositionTest, operator_largerEqual){
    BAM::TGenomePosition pos1(10, 200);
    BAM::TGenomePosition pos2(10, 100);
    // same chr, larger pos
    EXPECT_TRUE(pos1 >= pos2);

    // same chr, smaller pos
    pos2.move(300, 1);
    EXPECT_FALSE(pos1 >= pos2);

    // same chr, equal pos
    pos2.move(10, 200);
    EXPECT_TRUE(pos1 >= pos2);

    // larger chr
    pos2.move(9, 100);
    EXPECT_TRUE(pos1 >= pos2);

    // smaller chr
    pos2.move(11, 100);
    EXPECT_FALSE(pos1 >= pos2);
}

TEST(TGenomePositionTest, operator_stream){
    std::ostringstream os;
    BAM::TGenomePosition pos(10, 100);

    os << pos;
    EXPECT_EQ(os.str(), "10:100");
}

// TODO: test > and < against TGenomeWindow