//
// Created by madleina on 02.08.20.
//

#include "TGenomePosition.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

//-----------------------------------------------------
// TGenomePosition
//-----------------------------------------------------

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

TEST(TGenomePositionTest, operator_smallerWindow){
    BAM::TGenomeWindow window(1, 10, 20);

    // pos is larger than _to
    BAM::TGenomePosition pos(1, 25);
    EXPECT_TRUE(pos > window);

    // pos is equal to _to
    pos.move(1, 20);
    EXPECT_TRUE(pos > window);

    // pos is smaller than _to
    pos.move(1, 15);
    EXPECT_FALSE(pos > window);

    // pos is smaller than _from
    pos.move(1, 5);
    EXPECT_FALSE(pos > window);

    // pos is a smaller chromosome
    pos.move(0, 25);
    EXPECT_FALSE(pos > window);

    // pos is a larger chromosome
    pos.move(2, 25);
    EXPECT_TRUE(pos > window);
}

TEST(TGenomePositionTest, operator_largerWindow){
    BAM::TGenomeWindow window(1, 10, 20);

    // pos is smaller than _from
    BAM::TGenomePosition pos(1, 5);
    EXPECT_TRUE(pos < window);

    // pos is equal to _from
    pos.move(1, 10);
    EXPECT_FALSE(pos < window);

    // pos is larger than _from
    pos.move(1, 15);
    EXPECT_FALSE(pos < window);

    // pos is larger than _to
    pos.move(1, 25);
    EXPECT_FALSE(pos < window);

    // pos is a smaller chromosome
    pos.move(0, 25);
    EXPECT_TRUE(pos < window);

    // pos is a larger chromosome
    pos.move(2, 25);
    EXPECT_FALSE(pos < window);
}

TEST(TGenomePositionTest, operator_stream){
    std::ostringstream os;
    BAM::TGenomePosition pos(10, 100);

    os << pos;
    EXPECT_EQ(os.str(), "10:100");
}

//-----------------------------------------------------
// TGenomeWindow
//-----------------------------------------------------

TEST(TGenomeWindowTest, constructor_noParams){
    BAM::TGenomeWindow window;
    EXPECT_EQ(window.refID(), 0);
    EXPECT_EQ(window.fromOnChr(), 0);
    EXPECT_EQ(window.toOnChr(), 0);
}

TEST(TGenomeWindowTest, constructor_refID_start_end){
    BAM::TGenomeWindow window(10, 0, 100);
    EXPECT_EQ(window.refID(), 10);
    EXPECT_EQ(window.fromOnChr(), 0);
    EXPECT_EQ(window.toOnChr(), 100);
}

TEST(TGenomeWindowTest, constructor_refID_start){
    BAM::TGenomeWindow window(10, 100);
    EXPECT_EQ(window.refID(), 10);
    EXPECT_EQ(window.fromOnChr(), 100);
    EXPECT_EQ(window.toOnChr(), 101);
}


TEST(TGenomeWindowTest, constructor_positionFrom){
    BAM::TGenomePosition pos(10, 100);
    BAM::TGenomeWindow window(pos);
    EXPECT_EQ(window.refID(), 10);
    EXPECT_EQ(window.fromOnChr(), 100);
    EXPECT_EQ(window.toOnChr(), 101);
}

TEST(TGenomeWindowTest, constructor_positionFrom_positionTo){
    BAM::TGenomePosition from(10, 100);
    BAM::TGenomePosition to(10, 200);
    BAM::TGenomeWindow window(from, to);
    EXPECT_EQ(window.refID(), 10);
    EXPECT_EQ(window.fromOnChr(), 100);
    EXPECT_EQ(window.toOnChr(), 200);
}

TEST(TGenomeWindowTest, constructor_otherWindow){
    BAM::TGenomeWindow window(10, 100, 200);
    BAM::TGenomeWindow window2(window);
    EXPECT_EQ(window.refID(), 10);
    EXPECT_EQ(window.fromOnChr(), 100);
    EXPECT_EQ(window.toOnChr(), 200);
    EXPECT_EQ(window2.refID(), 10);
    EXPECT_EQ(window2.fromOnChr(), 100);
    EXPECT_EQ(window2.toOnChr(), 200);
}

TEST(TGenomeWindowTest, clear){
    BAM::TGenomeWindow window(10, 10, 100);
    window.clear();
    EXPECT_EQ(window.refID(), 0);
    EXPECT_EQ(window.fromOnChr(), 0);
    EXPECT_EQ(window.toOnChr(), 0);
}

TEST(TGenomeWindowTest, getters){
    BAM::TGenomeWindow window(10, 10, 100);

    BAM::TGenomePosition from = window.from();
    EXPECT_EQ(from.refID(), 10);
    EXPECT_EQ(from.position(), 10);
    EXPECT_EQ(window.fromOnChr(), 10);

    BAM::TGenomePosition to = window.to();
    EXPECT_EQ(to.refID(), 10);
    EXPECT_EQ(to.position(), 100);
    EXPECT_EQ(window.toOnChr(), 100);

    EXPECT_EQ(window.size(), 90);
}

TEST(TGenomeWindowTest, move_refId_from_to){
    BAM::TGenomeWindow window(10, 0, 100);

    // move ok
    window.move(1, 10, 11);
    EXPECT_EQ(window.refID(), 1);
    EXPECT_EQ(window.fromOnChr(), 10);
    EXPECT_EQ(window.toOnChr(), 11);

    // to == from -> throw
    EXPECT_THROW(window.move(1, 10, 10), std::runtime_error);
    // to > from -> throw
    EXPECT_THROW(window.move(1, 20, 10), std::runtime_error);
}

TEST(TGenomeWindowTest, move_from_length){
    BAM::TGenomeWindow window(10, 0, 100);

    // move ok
    BAM::TGenomePosition pos(1, 10);
    window.move(pos, 20);
    EXPECT_EQ(window.refID(), 1);
    EXPECT_EQ(window.fromOnChr(), 10);
    EXPECT_EQ(window.toOnChr(), 30);
}

TEST(TGenomeWindowTest, move_from_to){
    BAM::TGenomeWindow window(10, 0, 100);

    // move ok
    BAM::TGenomePosition from(1, 10);
    BAM::TGenomePosition to(1, 11);

    window.move(from, to);
    EXPECT_EQ(window.refID(), 1);
    EXPECT_EQ(window.fromOnChr(), 10);
    EXPECT_EQ(window.toOnChr(), 11);

    // not same chr -> throw
    to.move(2, 11);
    EXPECT_THROW(window.move(from, to), std::runtime_error);

    // to == from -> throw
    to.move(1, 10);
    EXPECT_THROW(window.move(from, to), std::runtime_error);

    // to > from -> throw
    to.move(1, 9);
    EXPECT_THROW(window.move(from, to), std::runtime_error);
}

TEST(TGenomeWindowTest, move_otherWindow){
    BAM::TGenomeWindow window1(10, 0, 100);
    BAM::TGenomeWindow window2(1, 10, 20);

    window1.move(window2);
    EXPECT_EQ(window1.refID(), 1);
    EXPECT_EQ(window1.fromOnChr(), 10);
    EXPECT_EQ(window1.toOnChr(), 20);
    EXPECT_EQ(window2.refID(), 1);
    EXPECT_EQ(window2.fromOnChr(), 10);
    EXPECT_EQ(window2.toOnChr(), 20);
}

TEST(TGenomeWindowTest, operator_equal){
    BAM::TGenomeWindow window1;
    BAM::TGenomeWindow window2(1, 10, 20);

    window1 = window2;
    EXPECT_EQ(window1.refID(), 1);
    EXPECT_EQ(window1.fromOnChr(), 10);
    EXPECT_EQ(window1.toOnChr(), 20);
    EXPECT_EQ(window2.refID(), 1);
    EXPECT_EQ(window2.fromOnChr(), 10);
    EXPECT_EQ(window2.toOnChr(), 20);
}

TEST(TGenomeWindowTest, operator_plus){
    BAM::TGenomeWindow window1(1, 10, 20);
    BAM::TGenomeWindow window2 = window1 + 10;

    EXPECT_EQ(window1.refID(), 1);
    EXPECT_EQ(window1.fromOnChr(), 10);
    EXPECT_EQ(window1.toOnChr(), 20);
    EXPECT_EQ(window2.refID(), 1);
    EXPECT_EQ(window2.fromOnChr(), 20);
    EXPECT_EQ(window2.toOnChr(), 30);
}

TEST(TGenomeWindowTest, operator_minus){
    BAM::TGenomeWindow window1(1, 10, 20);

    // this is ok
    BAM::TGenomeWindow window2 = window1 - 5;
    EXPECT_EQ(window1.refID(), 1); // window1 doesn't change
    EXPECT_EQ(window1.fromOnChr(), 10);
    EXPECT_EQ(window1.toOnChr(), 20);
    EXPECT_EQ(window2.refID(), 1);
    EXPECT_EQ(window2.fromOnChr(), 5);
    EXPECT_EQ(window2.toOnChr(), 15);

    // overshoot _from -> set _from to zero
    window2 = window1 - 11;
    EXPECT_EQ(window2.refID(), 1);
    EXPECT_EQ(window2.fromOnChr(), 0);
    EXPECT_EQ(window2.toOnChr(), 9);

    // overshoot _from and _to -> set _from to 0 and _to to 1
    window2 = window1 - 21;
    EXPECT_EQ(window2.refID(), 1);
    EXPECT_EQ(window2.fromOnChr(), 0);
    EXPECT_EQ(window2.toOnChr(), 1);
}

TEST(TGenomeWindowTest, sameChr){
    BAM::TGenomeWindow window1(1, 10, 20);
    BAM::TGenomeWindow window2(1, 100, 200);
    // same chr
    EXPECT_TRUE(window1.sameChr(window2));

    // other chr
    window2.move(11, 100, 200);
    EXPECT_FALSE(window1.sameChr(window2));
}

TEST(TGenomeWindowTest, contains_position){
    BAM::TGenomeWindow window(1, 10, 20);

    // these are ok
    BAM::TGenomePosition pos(1, 15);
    EXPECT_TRUE(window.contains(pos));
    pos.move(1, 10);
    EXPECT_TRUE(window.contains(pos));

    // position is outside
    pos.move(1, 20);
    EXPECT_FALSE(window.contains(pos));
    pos.move(1, 9);
    EXPECT_FALSE(window.contains(pos));

    // chromosome is different
    pos.move(2, 15);
    EXPECT_FALSE(window.contains(pos));
}

TEST(TGenomeWindowTest, contains_window){
    BAM::TGenomeWindow window1(1, 10, 20);

    // these are ok
    BAM::TGenomeWindow window2(1, 11, 19);
    EXPECT_TRUE(window1.contains(window2));
    window2.move(1, 10, 19);
    EXPECT_TRUE(window1.contains(window2));
    window2.move(1, 11, 20);
    EXPECT_TRUE(window1.contains(window2));
    window2.move(1, 10, 20);
    EXPECT_TRUE(window1.contains(window2));

    // positions are outside
    window2.move(1, 9, 19);
    EXPECT_FALSE(window1.contains(window2));
    window2.move(1, 10, 21);
    EXPECT_FALSE(window1.contains(window2));

    // chromosome is different
    window2.move(2, 11, 19);
    EXPECT_FALSE(window1.contains(window2));
}

TEST(TGenomeWindowTest, overlaps){
    BAM::TGenomeWindow window1(1, 10, 20);

    // fully within
    BAM::TGenomeWindow window2(1, 11, 19);
    EXPECT_TRUE(window1.overlaps(window2));
    window2.move(1, 10, 19);
    EXPECT_TRUE(window1.overlaps(window2));
    window2.move(1, 11, 20);
    EXPECT_TRUE(window1.overlaps(window2));
    window2.move(1, 10, 20);
    EXPECT_TRUE(window1.overlaps(window2));

    // overlaps on left side
    window2.move(1, 5, 11);
    EXPECT_TRUE(window1.overlaps(window2));

    // overlaps on right side
    window2.move(1, 19, 25);
    EXPECT_TRUE(window1.overlaps(window2));

    // overlaps on both sides
    window2.move(1, 9, 21);
    EXPECT_TRUE(window1.overlaps(window2));

    // positions don't overlap
    window2.move(1, 20, 25); // extends, but doesnt overlap
    EXPECT_FALSE(window1.overlaps(window2));
    window2.move(1, 5, 10);
    EXPECT_FALSE(window1.overlaps(window2)); // extends, but doesnt overlap
    window2.move(1, 5, 7);
    EXPECT_FALSE(window1.overlaps(window2));
    window2.move(1, 21, 25);
    EXPECT_FALSE(window1.overlaps(window2));

    // chromosome is different
    window2.move(2, 11, 19);
    EXPECT_FALSE(window1.overlaps(window2));
}

TEST(TGenomeWindowTest, overlapsOrExtends){
    BAM::TGenomeWindow window1(1, 10, 20);

    // fully within
    BAM::TGenomeWindow window2(1, 11, 19);
    EXPECT_TRUE(window1.overlapsOrExtends(window2));
    window2.move(1, 10, 19);
    EXPECT_TRUE(window1.overlapsOrExtends(window2));
    window2.move(1, 11, 20);
    EXPECT_TRUE(window1.overlapsOrExtends(window2));
    window2.move(1, 10, 20);
    EXPECT_TRUE(window1.overlapsOrExtends(window2));

    // overlaps on left side
    window2.move(1, 5, 11);
    EXPECT_TRUE(window1.overlapsOrExtends(window2));
    window2.move(1, 5, 10);
    EXPECT_TRUE(window1.overlapsOrExtends(window2));

    // overlaps on right side
    window2.move(1, 19, 25);
    EXPECT_TRUE(window1.overlapsOrExtends(window2));
    window2.move(1, 20, 25); // extends
    EXPECT_TRUE(window1.overlapsOrExtends(window2));

    // overlaps on both sides
    window2.move(1, 9, 21);
    EXPECT_TRUE(window1.overlapsOrExtends(window2));

    // positions don't overlap
    window2.move(1, 5, 7);
    EXPECT_FALSE(window1.overlapsOrExtends(window2));
    window2.move(1, 21, 25);
    EXPECT_FALSE(window1.overlapsOrExtends(window2));

    // chromosome is different
    window2.move(2, 11, 19);
    EXPECT_FALSE(window1.overlapsOrExtends(window2));
}

TEST(TGenomeWindowTest, mergeWith){
    BAM::TGenomeWindow window1(1, 10, 20);

    // fully within
    BAM::TGenomeWindow window2(1, 11, 19);
    EXPECT_TRUE(window1.mergeWith(window2));
    EXPECT_EQ(window1.refID(), 1);
    EXPECT_EQ(window1.fromOnChr(), 10);
    EXPECT_EQ(window1.toOnChr(), 20);
    EXPECT_EQ(window2.refID(), 1); // window2 shouldn't change
    EXPECT_EQ(window2.fromOnChr(), 11);
    EXPECT_EQ(window2.toOnChr(), 19);

    // overlaps on left side
    window1.move(1, 10, 20);
    window2.move(1, 5, 15);
    EXPECT_TRUE(window1.mergeWith(window2));
    EXPECT_EQ(window1.refID(), 1);
    EXPECT_EQ(window1.fromOnChr(), 5);
    EXPECT_EQ(window1.toOnChr(), 20);

    // extends on left side
    window1.move(1, 10, 20);
    window2.move(1, 5, 10);
    EXPECT_TRUE(window1.mergeWith(window2));
    EXPECT_EQ(window1.refID(), 1);
    EXPECT_EQ(window1.fromOnChr(), 5);
    EXPECT_EQ(window1.toOnChr(), 20);

    // overlaps on right side
    window1.move(1, 10, 20);
    window2.move(1, 15, 25);
    EXPECT_TRUE(window1.mergeWith(window2));
    EXPECT_EQ(window1.refID(), 1);
    EXPECT_EQ(window1.fromOnChr(), 10);
    EXPECT_EQ(window1.toOnChr(), 25);

    // extends on right side
    window1.move(1, 10, 20);
    window2.move(1, 20, 25);
    EXPECT_TRUE(window1.mergeWith(window2));
    EXPECT_EQ(window1.refID(), 1);
    EXPECT_EQ(window1.fromOnChr(), 10);
    EXPECT_EQ(window1.toOnChr(), 25);

    // overlaps on both sides
    window1.move(1, 10, 20);
    window2.move(1, 9, 21);
    EXPECT_TRUE(window1.mergeWith(window2));
    EXPECT_EQ(window1.refID(), 1);
    EXPECT_EQ(window1.fromOnChr(), 9);
    EXPECT_EQ(window1.toOnChr(), 21);

    // positions don't overlap -> nothing happens
    window1.move(1, 10, 20);
    window2.move(1, 5, 9);
    EXPECT_FALSE(window1.mergeWith(window2));
    EXPECT_EQ(window1.refID(), 1);
    EXPECT_EQ(window1.fromOnChr(), 10);
    EXPECT_EQ(window1.toOnChr(), 20);

    // chromosome is different -> nothing happens
    window1.move(1, 10, 20);
    window2.move(2, 19, 25);
    EXPECT_FALSE(window1.mergeWith(window2));
    EXPECT_EQ(window1.refID(), 1);
    EXPECT_EQ(window1.fromOnChr(), 10);
    EXPECT_EQ(window1.toOnChr(), 20);
}

TEST(TGenomeWindowTest, operator_plusEqual){
    BAM::TGenomeWindow window(1, 10, 20);
    window += 10;

    EXPECT_EQ(window.refID(), 1);
    EXPECT_EQ(window.fromOnChr(), 20);
    EXPECT_EQ(window.toOnChr(), 30);
}

TEST(TGenomeWindowTest, operator_minusEqual){
    BAM::TGenomeWindow window(1, 10, 20);

    // this is ok
    window -= 5;
    EXPECT_EQ(window.refID(), 1);
    EXPECT_EQ(window.fromOnChr(), 5);
    EXPECT_EQ(window.toOnChr(), 15);

    // overshoot _from -> set _from to zero
    window.move(1, 10, 20); // reset
    window -= 11;
    EXPECT_EQ(window.refID(), 1);
    EXPECT_EQ(window.fromOnChr(), 0);
    EXPECT_EQ(window.toOnChr(), 9);

    // overshoot _from and _to -> set _from to 0 and _to to 1
    window.move(1, 10, 20); // reset
    window -= 21;
    EXPECT_EQ(window.refID(), 1);
    EXPECT_EQ(window.fromOnChr(), 0);
    EXPECT_EQ(window.toOnChr(), 1);
}

TEST(TGenomeWindowTest, resize){
    BAM::TGenomeWindow window(1, 10, 20);
    window.resize(20);

    EXPECT_EQ(window.refID(), 1);
    EXPECT_EQ(window.fromOnChr(), 10);
    EXPECT_EQ(window.toOnChr(), 30);
}

TEST(TGenomeWindowTest, operator_smallerPosition){
    BAM::TGenomeWindow window(1, 10, 20);

    // pos is larger than _to
    BAM::TGenomePosition pos(1, 25);
    EXPECT_TRUE(window < pos);

    // pos is equal to _to
    pos.move(1, 20);
    EXPECT_TRUE(window < pos);

    // pos is smaller than _to
    pos.move(1, 15);
    EXPECT_FALSE(window < pos);

    // pos is smaller than _from
    pos.move(1, 5);
    EXPECT_FALSE(window < pos);

    // pos is a smaller chromosome
    pos.move(0, 25);
    EXPECT_FALSE(window < pos);

    // pos is a larger chromosome
    pos.move(2, 25);
    EXPECT_TRUE(window < pos);
}

TEST(TGenomeWindowTest, operator_largerPosition){
    BAM::TGenomeWindow window(1, 10, 20);

    // pos is smaller than _from
    BAM::TGenomePosition pos(1, 5);
    EXPECT_TRUE(window > pos);

    // pos is equal to _from
    pos.move(1, 10);
    EXPECT_FALSE(window > pos);

    // pos is larger than _from
    pos.move(1, 15);
    EXPECT_FALSE(window > pos);

    // pos is larger than _to
    pos.move(1, 25);
    EXPECT_FALSE(window > pos);

    // pos is a smaller chromosome
    pos.move(0, 25);
    EXPECT_TRUE(window > pos);

    // pos is a larger chromosome
    pos.move(2, 25);
    EXPECT_FALSE(window > pos);
}

TEST(TGenomeWindowTest, operator_smallerWindow){
    BAM::TGenomeWindow window1(1, 10, 20);

    // window1.from is smaller than window2.from
    BAM::TGenomeWindow window2(1, 15, 18);
    EXPECT_TRUE(window1 < window2);

    // window1.from is equal to window2.from
    window2.move(1, 10, 15);
    EXPECT_FALSE(window1 < window2);

    // window1.from is larger than window2.from
    window2.move(1, 5, 15);
    EXPECT_FALSE(window1 < window2);

    // window1 is on a larger chromosome than window2
    window2.move(0, 25, 30);
    EXPECT_FALSE(window1 < window2);

    // window1 is on a smaller chromosome than window2
    window2.move(2, 25, 30);
    EXPECT_TRUE(window1 < window2);
}

TEST(TGenomeWindowTest, operator_largerWindow){
    BAM::TGenomeWindow window1(1, 10, 20);

    // window1.from is larger than window2.from
    BAM::TGenomeWindow window2(1, 5, 15);
    EXPECT_TRUE(window1 > window2);

    // window1.from is equal to window2.from
    window2.move(1, 10, 15);
    EXPECT_FALSE(window1 > window2);

    // window1.from is smaller than window2.from
    window2.move(1, 12, 15);
    EXPECT_FALSE(window1 > window2);

    // window1 is on a smaller chromosome than window2
    window2.move(2, 25, 30);
    EXPECT_FALSE(window1 > window2);

    // window1 is on a larger chromosome than window2
    window2.move(0, 25, 30);
    EXPECT_TRUE(window1 > window2);
}

TEST(TGenomeWindowTest, merge){
    BAM::TGenomeWindow window1(1, 10, 20);

    // fully within
    BAM::TGenomeWindow window2(1, 11, 19);
    BAM::TGenomeWindow res = merge(window1, window2);
    EXPECT_EQ(res.refID(), 1);
    EXPECT_EQ(res.fromOnChr(), 10);
    EXPECT_EQ(res.toOnChr(), 20);

    // overlaps on left side
    window2.move(1, 5, 15);
    res = merge(window1, window2);
    EXPECT_EQ(res.refID(), 1);
    EXPECT_EQ(res.fromOnChr(), 5);
    EXPECT_EQ(res.toOnChr(), 20);

    // extends on left side
    window2.move(1, 5, 10);
    res = merge(window1, window2);
    EXPECT_EQ(res.refID(), 1);
    EXPECT_EQ(res.fromOnChr(), 5);
    EXPECT_EQ(res.toOnChr(), 20);

    // overlaps on right side
    window2.move(1, 15, 25);
    res = merge(window1, window2);
    EXPECT_EQ(res.refID(), 1);
    EXPECT_EQ(res.fromOnChr(), 10);
    EXPECT_EQ(res.toOnChr(), 25);

    // extends on right side
    window2.move(1, 20, 25);
    res = merge(window1, window2);
    EXPECT_EQ(res.refID(), 1);
    EXPECT_EQ(res.fromOnChr(), 10);
    EXPECT_EQ(res.toOnChr(), 25);

    // overlaps on both sides
    window2.move(1, 9, 21);
    res = merge(window1, window2);
    EXPECT_EQ(res.refID(), 1);
    EXPECT_EQ(res.fromOnChr(), 9);
    EXPECT_EQ(res.toOnChr(), 21);

    // positions don't overlap -> nothing happens
    window2.move(1, 5, 9);
    EXPECT_THROW(merge(window1, window2), std::runtime_error);

    // chromosome is different -> nothing happens
    window2.move(2, 19, 25);
    EXPECT_THROW(merge(window1, window2), std::runtime_error);
}

TEST(TGenomeWindowTest, operator_stream){
    std::ostringstream os;
    BAM::TGenomeWindow window(10, 100, 200);

    os << window;
    EXPECT_EQ(os.str(), "10:100-10:200");
}
