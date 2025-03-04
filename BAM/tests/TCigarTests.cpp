#include "gtest/gtest.h"
#include "TCigar.h"

using BAM::TCigar;

TEST(TCigarTests, properties) {
	TCigar cigar;
	cigar.add('S', 1);
	cigar.add('M', 2);
	cigar.add('I', 3);
	cigar.add('=', 4);
	cigar.add('D', 5);
	cigar.add('X', 6);
	cigar.add('N', 7);
	cigar.add('S', 8);

	EXPECT_EQ(cigar.compileString(), "1S2M3I4=5D6X7N8S");

	// one type
	EXPECT_EQ(cigar.lengthSoftClippedLeft(), 1);
	EXPECT_EQ(cigar.lengthSoftClippedRight(), 8);
	EXPECT_EQ(cigar.lengthInserted(), 3);
	EXPECT_EQ(cigar.lengthDeleted(), 5);
	EXPECT_EQ(cigar.lengthSkipped(), 7);

	// combined
	EXPECT_EQ(cigar.lengthSoftClipped(), 1 + 8);
	EXPECT_EQ(cigar.lengthAligned(), 2 + 4 + 6);
	EXPECT_EQ(cigar.lengthMapped(), 2 + 4 + 6 + 5 + 7);
	EXPECT_EQ(cigar.lengthRead(), 1 + 2 + 3 + 4 + 6 + 8);
	EXPECT_EQ(cigar.lengthSequenced(), 2 + 3 + 4 + 6);

}

TEST(TCigarTests, remove) {
	TCigar cigar;
	cigar.add('S', 5);
	cigar.add('M', 2);
	cigar.add('I', 3);
	cigar.add('=', 4);
	cigar.add('D', 5);
	cigar.add('X', 6);
	cigar.add('N', 7);
	cigar.add('S', 8);

	cigar.trimSoftClips(5);
	// one type
	EXPECT_EQ(cigar.lengthSoftClippedLeft(), 5);
	EXPECT_EQ(cigar.lengthSoftClippedRight(), 5);
	EXPECT_EQ(cigar.lengthInserted(), 3);
	EXPECT_EQ(cigar.lengthDeleted(), 5);
	EXPECT_EQ(cigar.lengthSkipped(), 7);

	// combined
	EXPECT_EQ(cigar.lengthSoftClipped(), 5 + 5);
	EXPECT_EQ(cigar.lengthAligned(), 2 + 4 + 6);
	EXPECT_EQ(cigar.lengthMapped(), 2 + 4 + 6 + 5 + 7);
	EXPECT_EQ(cigar.lengthRead(), 5 + 2 + 3 + 4 + 6 + 5);
	EXPECT_EQ(cigar.lengthSequenced(), 2 + 3 + 4 + 6);

	cigar.trimSoftClips(2);
	// one type
	EXPECT_EQ(cigar.lengthSoftClippedLeft(), 2);
	EXPECT_EQ(cigar.lengthSoftClippedRight(), 2);
	EXPECT_EQ(cigar.lengthInserted(), 3);
	EXPECT_EQ(cigar.lengthDeleted(), 5);
	EXPECT_EQ(cigar.lengthSkipped(), 7);

	// combined
	EXPECT_EQ(cigar.lengthSoftClipped(), 2 + 2);
	EXPECT_EQ(cigar.lengthAligned(), 2 + 4 + 6);
	EXPECT_EQ(cigar.lengthMapped(), 2 + 4 + 6 + 5 + 7);
	EXPECT_EQ(cigar.lengthRead(), 2 + 2 + 3 + 4 + 6 + 2);
	EXPECT_EQ(cigar.lengthSequenced(), 2 + 3 + 4 + 6);

	cigar.removeSoftClips();

	// one type
	EXPECT_EQ(cigar.lengthSoftClippedLeft(), 0);
	EXPECT_EQ(cigar.lengthSoftClippedRight(), 0);
	EXPECT_EQ(cigar.lengthInserted(), 3);
	EXPECT_EQ(cigar.lengthDeleted(), 5);
	EXPECT_EQ(cigar.lengthSkipped(), 7);

	// combined
	EXPECT_EQ(cigar.lengthSoftClipped(), 0);
	EXPECT_EQ(cigar.lengthAligned(), 2 + 4 + 6);
	EXPECT_EQ(cigar.lengthMapped(), 2 + 4 + 6 + 5 + 7);
	EXPECT_EQ(cigar.lengthRead(), 2 + 3 + 4 + 6);
	EXPECT_EQ(cigar.lengthSequenced(), 2 + 3 + 4 + 6);
}

TEST(TCigarTests, add) {
	TCigar cigar;
	cigar.add('S', 1);
	cigar.add('M', 2);
	cigar.add('I', 3);
	cigar.add('=', 4);
	cigar.add('D', 5);
	cigar.add('X', 6);
	cigar.add('N', 7);
	cigar.add('S', 8);

	EXPECT_EQ(cigar.compileString(), "1S2M3I4=5D6X7N8S");

	cigar.addSoftClipsLeft(1);
	EXPECT_EQ(cigar.compileString(), "2S1M3I4=5D6X7N8S");

	cigar.addSoftClipsRight(2);
	EXPECT_EQ(cigar.compileString(), "2S1M3I4=5D6X5N8S");

	cigar.addSoftClipsLeft(3);
	EXPECT_EQ(cigar.compileString(), "8S2=5D6X5N8S");

	cigar.addSoftClipsRight(4);
	EXPECT_EQ(cigar.compileString(), "8S2=5D6X1N8S");

	cigar.addSoftClipsLeft(5);
	EXPECT_EQ(cigar.compileString(), "10S2D6X1N8S");

	cigar.addSoftClipsRight(6);
	EXPECT_EQ(cigar.compileString(), "10S2D1X13S");

	EXPECT_ANY_THROW(cigar.addSoftClipsLeft(7));
	EXPECT_ANY_THROW(cigar.addSoftClipsRight(7));
}
