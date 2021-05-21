/*
 * TTypesTests.cpp
 *
 *  Created on: May 10, 2021
 *      Author: phaentu
 */


#include "../../commonutilities/tests/TestCase.h"
#include "Types.h"

using namespace testing;
using namespace BAM;

//------------------------------------------------
// Base
//------------------------------------------------

TEST(TypesTests, Base){
    //construct from char
	EXPECT_EQ(Base('A'), BaseEnum::A);
	EXPECT_EQ(Base('C'), BaseEnum::C);
	EXPECT_EQ(Base('G'), BaseEnum::G);
	EXPECT_EQ(Base('T'), BaseEnum::T);
	EXPECT_EQ(Base('N'), BaseEnum::N);

	EXPECT_EQ(Base('a'), BaseEnum::A);
	EXPECT_EQ(Base('c'), BaseEnum::C);
	EXPECT_EQ(Base('g'), BaseEnum::G);
	EXPECT_EQ(Base('t'), BaseEnum::T);
	EXPECT_EQ(Base('n'), BaseEnum::N);

	EXPECT_EQ(Base('F'), BaseEnum::N);
	EXPECT_EQ(Base('q'), BaseEnum::N);

	//construct from enum
	EXPECT_EQ(Base(BaseEnum::A), BaseEnum::A);
	EXPECT_EQ(Base(BaseEnum::C), BaseEnum::C);
	EXPECT_EQ(Base(BaseEnum::G), BaseEnum::G);
	EXPECT_EQ(Base(BaseEnum::T), BaseEnum::T);
	EXPECT_EQ(Base(BaseEnum::N), BaseEnum::N);

	//cast to char
	EXPECT_EQ((char) Base(BaseEnum::A), 'A');
	EXPECT_EQ((char) Base(BaseEnum::C), 'C');
	EXPECT_EQ((char) Base(BaseEnum::G), 'G');
	EXPECT_EQ((char) Base(BaseEnum::T), 'T');
	EXPECT_EQ((char) Base(BaseEnum::N), 'N');
};

//------------------------------------------------
// Genotype
//------------------------------------------------

TEST(TypesTests, Genotype){
	//construct from Base
	Base A_(BaseEnum::A);
	Base C_(BaseEnum::C);
	Base G_(BaseEnum::G);
	Base T_(BaseEnum::T);
	Base N_(BaseEnum::N);

	EXPECT_EQ(Genotype(A_, A_), GenotypeEnum::AA);
	EXPECT_EQ(Genotype(A_, C_), GenotypeEnum::AC);
	EXPECT_EQ(Genotype(A_, G_), GenotypeEnum::AG);
	EXPECT_EQ(Genotype(A_, T_), GenotypeEnum::AT);

	EXPECT_EQ(Genotype(C_, A_), GenotypeEnum::AC);
	EXPECT_EQ(Genotype(C_, C_), GenotypeEnum::CC);
	EXPECT_EQ(Genotype(C_, G_), GenotypeEnum::CG);
	EXPECT_EQ(Genotype(C_, T_), GenotypeEnum::CT);

	EXPECT_EQ(Genotype(G_, A_), GenotypeEnum::AG);
	EXPECT_EQ(Genotype(G_, C_), GenotypeEnum::CG);
	EXPECT_EQ(Genotype(G_, G_), GenotypeEnum::GG);
	EXPECT_EQ(Genotype(G_, T_), GenotypeEnum::GT);

	EXPECT_EQ(Genotype(T_, A_), GenotypeEnum::AT);
	EXPECT_EQ(Genotype(T_, C_), GenotypeEnum::CT);
	EXPECT_EQ(Genotype(T_, G_), GenotypeEnum::GT);
	EXPECT_EQ(Genotype(T_, T_), GenotypeEnum::TT);

	EXPECT_EQ(Genotype(N_, N_), GenotypeEnum::NN);
	EXPECT_EQ(Genotype(A_, N_), GenotypeEnum::NN);
	EXPECT_EQ(Genotype(N_, C_), GenotypeEnum::NN);

	//construct from BaseEnum
	EXPECT_EQ(Genotype(BaseEnum::A, BaseEnum::A), GenotypeEnum::AA);
	EXPECT_EQ(Genotype(BaseEnum::A, BaseEnum::C), GenotypeEnum::AC);
	EXPECT_EQ(Genotype(BaseEnum::A, BaseEnum::G), GenotypeEnum::AG);
	EXPECT_EQ(Genotype(BaseEnum::A, BaseEnum::T), GenotypeEnum::AT);

	EXPECT_EQ(Genotype(BaseEnum::C, BaseEnum::A), GenotypeEnum::AC);
	EXPECT_EQ(Genotype(BaseEnum::C, BaseEnum::C), GenotypeEnum::CC);
	EXPECT_EQ(Genotype(BaseEnum::C, BaseEnum::G), GenotypeEnum::CG);
	EXPECT_EQ(Genotype(BaseEnum::C, BaseEnum::T), GenotypeEnum::CT);

	EXPECT_EQ(Genotype(BaseEnum::G, BaseEnum::A), GenotypeEnum::AG);
	EXPECT_EQ(Genotype(BaseEnum::G, BaseEnum::C), GenotypeEnum::CG);
	EXPECT_EQ(Genotype(BaseEnum::G, BaseEnum::G), GenotypeEnum::GG);
	EXPECT_EQ(Genotype(BaseEnum::G, BaseEnum::T), GenotypeEnum::GT);

	EXPECT_EQ(Genotype(BaseEnum::T, BaseEnum::A), GenotypeEnum::AT);
	EXPECT_EQ(Genotype(BaseEnum::T, BaseEnum::C), GenotypeEnum::CT);
	EXPECT_EQ(Genotype(BaseEnum::T, BaseEnum::G), GenotypeEnum::GT);
	EXPECT_EQ(Genotype(BaseEnum::T, BaseEnum::T), GenotypeEnum::TT);

	EXPECT_EQ(Genotype(BaseEnum::N, BaseEnum::N), GenotypeEnum::NN);
	EXPECT_EQ(Genotype(BaseEnum::A, BaseEnum::N), GenotypeEnum::NN);
	EXPECT_EQ(Genotype(BaseEnum::N, BaseEnum::C), GenotypeEnum::NN);

    //cast to string
	EXPECT_EQ( (std::string) Genotype(BaseEnum::A, BaseEnum::A), std::string("AA"));
	EXPECT_EQ( (std::string) Genotype(BaseEnum::A, BaseEnum::C), std::string("AC"));
	EXPECT_EQ( (std::string) Genotype(BaseEnum::A, BaseEnum::G), std::string("AG"));
	EXPECT_EQ( (std::string) Genotype(BaseEnum::A, BaseEnum::T), std::string("AT"));
	EXPECT_EQ( (std::string) Genotype(BaseEnum::C, BaseEnum::A), std::string("AC"));

	EXPECT_EQ( (std::string) Genotype(BaseEnum::C, BaseEnum::A), std::string("AC"));
	EXPECT_EQ( (std::string) Genotype(BaseEnum::C, BaseEnum::C), std::string("CC"));
	EXPECT_EQ( (std::string) Genotype(BaseEnum::C, BaseEnum::G), std::string("CG"));
	EXPECT_EQ( (std::string) Genotype(BaseEnum::C, BaseEnum::T), std::string("CT"));

	EXPECT_EQ( (std::string) Genotype(BaseEnum::G, BaseEnum::A), std::string("AG"));
	EXPECT_EQ( (std::string) Genotype(BaseEnum::G, BaseEnum::C), std::string("CG"));
	EXPECT_EQ( (std::string) Genotype(BaseEnum::G, BaseEnum::G), std::string("GG"));
	EXPECT_EQ( (std::string) Genotype(BaseEnum::G, BaseEnum::T), std::string("GT"));

	EXPECT_EQ( (std::string) Genotype(BaseEnum::T, BaseEnum::A), std::string("AT"));
	EXPECT_EQ( (std::string) Genotype(BaseEnum::T, BaseEnum::C), std::string("CT"));
	EXPECT_EQ( (std::string) Genotype(BaseEnum::T, BaseEnum::G), std::string("GT"));
	EXPECT_EQ( (std::string) Genotype(BaseEnum::T, BaseEnum::T), std::string("TT"));

	EXPECT_EQ( (std::string) Genotype(BaseEnum::N, BaseEnum::N), std::string("NN"));
	EXPECT_EQ( (std::string) Genotype(BaseEnum::A, BaseEnum::N), std::string("NN"));
	EXPECT_EQ( (std::string) Genotype(BaseEnum::N, BaseEnum::C), std::string("NN"));
};

//------------------------------------------------
// BaseContext
//------------------------------------------------

TEST(TypesTests, BaseContext){
	//construct from Base
	Base A_(BaseEnum::A);
	Base C_(BaseEnum::C);
	Base G_(BaseEnum::G);
	Base T_(BaseEnum::T);
	Base N_(BaseEnum::N);

	EXPECT_EQ(BaseContext(A_, A_), BaseContextEnum::cAA);
	EXPECT_EQ(BaseContext(A_, C_), BaseContextEnum::cAC);
	EXPECT_EQ(BaseContext(A_, G_), BaseContextEnum::cAG);
	EXPECT_EQ(BaseContext(A_, T_), BaseContextEnum::cAT);

	EXPECT_EQ(BaseContext(C_, A_), BaseContextEnum::cCA);
	EXPECT_EQ(BaseContext(C_, C_), BaseContextEnum::cCC);
	EXPECT_EQ(BaseContext(C_, G_), BaseContextEnum::cCG);
	EXPECT_EQ(BaseContext(C_, T_), BaseContextEnum::cCT);

	EXPECT_EQ(BaseContext(G_, A_), BaseContextEnum::cGA);
	EXPECT_EQ(BaseContext(G_, C_), BaseContextEnum::cGC);
	EXPECT_EQ(BaseContext(G_, G_), BaseContextEnum::cGG);
	EXPECT_EQ(BaseContext(G_, T_), BaseContextEnum::cGT);

	EXPECT_EQ(BaseContext(T_, A_), BaseContextEnum::cTA);
	EXPECT_EQ(BaseContext(T_, C_), BaseContextEnum::cTC);
	EXPECT_EQ(BaseContext(T_, G_), BaseContextEnum::cTG);
	EXPECT_EQ(BaseContext(T_, T_), BaseContextEnum::cTT);

	EXPECT_EQ(BaseContext(N_, A_), BaseContextEnum::cNA);
	EXPECT_EQ(BaseContext(N_, C_), BaseContextEnum::cNC);
	EXPECT_EQ(BaseContext(N_, G_), BaseContextEnum::cNG);
	EXPECT_EQ(BaseContext(N_, T_), BaseContextEnum::cNT);

	EXPECT_EQ(BaseContext(A_, N_), BaseContextEnum::cNN);
	EXPECT_EQ(BaseContext(C_, N_), BaseContextEnum::cNN);
	EXPECT_EQ(BaseContext(G_, N_), BaseContextEnum::cNN);
	EXPECT_EQ(BaseContext(T_, N_), BaseContextEnum::cNN);
	EXPECT_EQ(BaseContext(N_, N_), BaseContextEnum::cNN);

	//cast to string
	EXPECT_EQ( (std::string) BaseContext(BaseEnum::A, BaseEnum::A), std::string("AA"));
	EXPECT_EQ( (std::string) BaseContext(BaseEnum::A, BaseEnum::C), std::string("AC"));
	EXPECT_EQ( (std::string) BaseContext(BaseEnum::A, BaseEnum::G), std::string("AG"));
	EXPECT_EQ( (std::string) BaseContext(BaseEnum::A, BaseEnum::T), std::string("AT"));

	EXPECT_EQ( (std::string) BaseContext(BaseEnum::C, BaseEnum::A), std::string("CA"));
	EXPECT_EQ( (std::string) BaseContext(BaseEnum::C, BaseEnum::C), std::string("CC"));
	EXPECT_EQ( (std::string) BaseContext(BaseEnum::C, BaseEnum::G), std::string("CG"));
	EXPECT_EQ( (std::string) BaseContext(BaseEnum::C, BaseEnum::T), std::string("CT"));

	EXPECT_EQ( (std::string) BaseContext(BaseEnum::G, BaseEnum::A), std::string("GA"));
	EXPECT_EQ( (std::string) BaseContext(BaseEnum::G, BaseEnum::C), std::string("GC"));
	EXPECT_EQ( (std::string) BaseContext(BaseEnum::G, BaseEnum::G), std::string("GG"));
	EXPECT_EQ( (std::string) BaseContext(BaseEnum::G, BaseEnum::T), std::string("GT"));

	EXPECT_EQ( (std::string) BaseContext(BaseEnum::T, BaseEnum::A), std::string("TA"));
	EXPECT_EQ( (std::string) BaseContext(BaseEnum::T, BaseEnum::C), std::string("TC"));
	EXPECT_EQ( (std::string) BaseContext(BaseEnum::T, BaseEnum::G), std::string("TG"));
	EXPECT_EQ( (std::string) BaseContext(BaseEnum::T, BaseEnum::T), std::string("TT"));

	EXPECT_EQ( (std::string) BaseContext(BaseEnum::A, BaseEnum::N), std::string("NN"));
	EXPECT_EQ( (std::string) BaseContext(BaseEnum::C, BaseEnum::N), std::string("NN"));
	EXPECT_EQ( (std::string) BaseContext(BaseEnum::G, BaseEnum::N), std::string("NN"));
	EXPECT_EQ( (std::string) BaseContext(BaseEnum::T, BaseEnum::N), std::string("NN"));
	EXPECT_EQ( (std::string) BaseContext(BaseEnum::N, BaseEnum::N), std::string("NN"));
};

//------------------------------------------------
// BaseContext
//------------------------------------------------

TEST(TypesTests, AllelicCombination){
	//construct from Base
	Base A_(BaseEnum::A);
	Base C_(BaseEnum::C);
	Base G_(BaseEnum::G);
	Base T_(BaseEnum::T);
	Base N_(BaseEnum::N);

	EXPECT_EQ(AllelicCombination(A_, C_), AllelicCombinationEnum::alleleicCombinationAC);
	EXPECT_EQ(AllelicCombination(A_, G_), AllelicCombinationEnum::alleleicCombinationAG);
	EXPECT_EQ(AllelicCombination(A_, T_), AllelicCombinationEnum::alleleicCombinationAT);
	EXPECT_EQ(AllelicCombination(C_, G_), AllelicCombinationEnum::alleleicCombinationCG);
	EXPECT_EQ(AllelicCombination(C_, T_), AllelicCombinationEnum::alleleicCombinationCT);
	EXPECT_EQ(AllelicCombination(G_, T_), AllelicCombinationEnum::alleleicCombinationGT);

	EXPECT_EQ(AllelicCombination(C_, A_), AllelicCombinationEnum::alleleicCombinationAC);
	EXPECT_EQ(AllelicCombination(G_, A_), AllelicCombinationEnum::alleleicCombinationAG);
	EXPECT_EQ(AllelicCombination(T_, A_), AllelicCombinationEnum::alleleicCombinationAT);
	EXPECT_EQ(AllelicCombination(G_, C_), AllelicCombinationEnum::alleleicCombinationCG);
	EXPECT_EQ(AllelicCombination(T_, C_), AllelicCombinationEnum::alleleicCombinationCT);
	EXPECT_EQ(AllelicCombination(T_, G_), AllelicCombinationEnum::alleleicCombinationGT);

	EXPECT_EQ(AllelicCombination(A_, N_), AllelicCombinationEnum::alleleicCombinationNN);
	EXPECT_EQ(AllelicCombination(C_, N_), AllelicCombinationEnum::alleleicCombinationNN);
	EXPECT_EQ(AllelicCombination(G_, N_), AllelicCombinationEnum::alleleicCombinationNN);
	EXPECT_EQ(AllelicCombination(T_, N_), AllelicCombinationEnum::alleleicCombinationNN);
	EXPECT_EQ(AllelicCombination(N_, A_), AllelicCombinationEnum::alleleicCombinationNN);
	EXPECT_EQ(AllelicCombination(N_, C_), AllelicCombinationEnum::alleleicCombinationNN);
	EXPECT_EQ(AllelicCombination(N_, G_), AllelicCombinationEnum::alleleicCombinationNN);
	EXPECT_EQ(AllelicCombination(N_, T_), AllelicCombinationEnum::alleleicCombinationNN);
	EXPECT_EQ(AllelicCombination(N_, N_), AllelicCombinationEnum::alleleicCombinationNN);

	EXPECT_ANY_THROW(AllelicCombination(A_, A_));
	EXPECT_ANY_THROW(AllelicCombination(C_, C_));
	EXPECT_ANY_THROW(AllelicCombination(G_, G_));
	EXPECT_ANY_THROW(AllelicCombination(T_, T_));

	EXPECT_EQ(AllelicCombination(A_, C_).firstAllele(), BaseEnum::A);
	EXPECT_EQ(AllelicCombination(A_, C_).secondAllele(), BaseEnum::C);
	EXPECT_EQ(AllelicCombination(A_, G_).firstAllele(), BaseEnum::A);
	EXPECT_EQ(AllelicCombination(A_, G_).secondAllele(), BaseEnum::G);
	EXPECT_EQ(AllelicCombination(A_, T_).firstAllele(), BaseEnum::A);
	EXPECT_EQ(AllelicCombination(A_, T_).secondAllele(), BaseEnum::T);
	EXPECT_EQ(AllelicCombination(C_, G_).firstAllele(), BaseEnum::C);
	EXPECT_EQ(AllelicCombination(C_, G_).secondAllele(), BaseEnum::G);
	EXPECT_EQ(AllelicCombination(C_, T_).firstAllele(), BaseEnum::C);
	EXPECT_EQ(AllelicCombination(C_, T_).secondAllele(), BaseEnum::T);
	EXPECT_EQ(AllelicCombination(G_, T_).firstAllele(), BaseEnum::G);
	EXPECT_EQ(AllelicCombination(G_, T_).secondAllele(), BaseEnum::T);

	EXPECT_EQ(AllelicCombination(A_, C_).homoFirst(), GenotypeEnum::AA);
	EXPECT_EQ(AllelicCombination(A_, C_).het(), GenotypeEnum::AC);
	EXPECT_EQ(AllelicCombination(A_, C_).homoSecond(), GenotypeEnum::CC);

	EXPECT_EQ(AllelicCombination(G_, C_).homoFirst(), GenotypeEnum::CC);
	EXPECT_EQ(AllelicCombination(G_, C_).het(), GenotypeEnum::CG);
	EXPECT_EQ(AllelicCombination(G_, C_).homoSecond(), GenotypeEnum::GG);

	EXPECT_EQ(AllelicCombination(T_, A_).homoFirst(), GenotypeEnum::AA);
	EXPECT_EQ(AllelicCombination(T_, A_).het(), GenotypeEnum::AT);
	EXPECT_EQ(AllelicCombination(T_, A_).homoSecond(), GenotypeEnum::TT);
};

//------------------------------------------------
// ErroRate, PhredErrorRate and Quality
//------------------------------------------------

TEST(TypesTests, ErrorRate){
	//from / to quality
	for(uint8_t q = 33; q < 127 ; ++q){
		BaseQuality qual(q);
		PhredIntErrorRate phredInt(qual);
		EXPECT_EQ(qual, BaseQuality(phredInt));

		HighPrecisionPhredIntErrorRate highPhredInt(qual);
		EXPECT_EQ(qual, BaseQuality(highPhredInt));

		PhredErrorRate phred(qual);
		EXPECT_EQ(qual, BaseQuality(phred));

		ErrorRate err(qual);
		EXPECT_EQ(qual, BaseQuality(err));

		LogErrorRate lerr(qual);
		EXPECT_EQ(qual, BaseQuality(lerr));
	}

	//from / to phredInt
	for(uint8_t p = PhredIntErrorRate::min(); p < PhredIntErrorRate::max(); ++p){
		PhredIntErrorRate phredInt(p);

		HighPrecisionPhredIntErrorRate highPhredInt(phredInt);
		EXPECT_EQ(phredInt, PhredIntErrorRate(highPhredInt));

		PhredErrorRate phred(phredInt);
		EXPECT_EQ(phredInt, PhredIntErrorRate(phred));

		ErrorRate err(phredInt);
		EXPECT_EQ(phredInt, PhredIntErrorRate(err));

		LogErrorRate lerr(phredInt);
		EXPECT_EQ(phredInt, PhredIntErrorRate(lerr));
	}

	//from / to highPhredInt
	//use phredInt to populate
	for(uint8_t p = PhredIntErrorRate::min(); p < PhredIntErrorRate::max(); ++p){
		PhredIntErrorRate phredInt(p);
		HighPrecisionPhredIntErrorRate highPhredInt(phredInt);

		PhredErrorRate phred(highPhredInt);
		EXPECT_EQ(highPhredInt, HighPrecisionPhredIntErrorRate(phred));

		ErrorRate err(highPhredInt);
		EXPECT_EQ(highPhredInt, HighPrecisionPhredIntErrorRate(err));

		LogErrorRate lerr(highPhredInt);
		EXPECT_EQ(highPhredInt, HighPrecisionPhredIntErrorRate(lerr));
	}

	//from / to phred
	//use phredInt to populate
	for(uint8_t p = 0; p < 255; ++p){
		PhredIntErrorRate phredInt(p);
		PhredErrorRate phred(phredInt);

		ErrorRate err(phred);
		EXPECT_DOUBLE_EQ((double) phred, (double) PhredErrorRate(err));

		LogErrorRate lerr(phred);
		EXPECT_DOUBLE_EQ((double) phred, (double) PhredErrorRate(lerr));
	}

	//from / to error
	//use phredInt to populate
	for(uint8_t p = 0; p < 255; ++p){
		PhredIntErrorRate phredInt(p);
		ErrorRate error(phredInt);

		LogErrorRate lerr(error);
		EXPECT_NEAR( (double) error, (double) ErrorRate(lerr), 10E-15);
	}
};





