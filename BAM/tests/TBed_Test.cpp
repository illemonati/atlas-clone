//
// Created by caduffm on 8/4/20.
//
#include "genometools/BED/TBed.h"
#include "gtest/gtest.h"

#include <stdio.h>
#include <map>
#include <ostream>
#include <stdexcept>
#include <string>

#include "genometools/GenomePositions/TChromosomes.h"
#include "coretools/Files/TOutputFile.h"
#include "genometools/GenomePositions/TGenomePosition.h"

//-------------------------------------------------------------
// TBedChromosome
//-------------------------------------------------------------
//
TEST(TBedChromosomeTest, constructor){
    genometools::TBedChromosome chr(1, "1");

    EXPECT_EQ(chr.refId, 1);
    EXPECT_EQ(chr.name, "1");
}

TEST(TBedChromosomeTest, operator_smaller){
    genometools::TBedChromosome chr1(1, "chr1");
    genometools::TBedChromosome chr2(2, "chr2");

    // bool operator<(const TBedChromosome & Other) const
    EXPECT_TRUE(chr1 < chr2);

    // bool operator<(const std::string & Name) const
    EXPECT_TRUE(chr1 < "chr2");

    // bool operator<(const std::string & Name, const TBedChromosome & Chr)
    EXPECT_TRUE("chr1" < chr2);
}

/*
//-------------------------------------------------------------
// TBed_base
//-------------------------------------------------------------

TEST(TBed_baseTest, addChromosome){
    genometools::TBed_base bedBase;

    // add chromosome and check if name is correct
    auto it = bedBase.addChromosome("chr1");
    EXPECT_EQ(it->name, "chr1");
    EXPECT_EQ(it->refId, 0);

    // add another chromosome
    it = bedBase.addChromosome("chr5");
    EXPECT_EQ(it->name, "chr5");
    EXPECT_EQ(it->refId, 1);

    // add another chromosome, not ordered
    it = bedBase.addChromosome("chr3");
    EXPECT_EQ(it->name, "chr3");
    EXPECT_EQ(it->refId, 2);

    // add chromosome that already exists
    it = bedBase.addChromosome("chr5");
    EXPECT_EQ(it->name, "chr5");
    EXPECT_EQ(it->refId, 1);
}

TEST(TBed_baseTest, addChromosomes){
    genometools::TBed_base bedBase;

    // create TChromosomes
    genometools::TChromosomes chrs;
    chrs.appendChromosome("chr1", 10);
    chrs.appendChromosome("chr5", 10);
    chrs.appendChromosome("chr3", 10);

    // add it to bedBase
    bedBase.addChromosomes(chrs);

    // get RefID and check if it is correct
    EXPECT_EQ(bedBase.getRefID("chr1"), 0);
    EXPECT_EQ(bedBase.getRefID("chr5"), 1);
    EXPECT_EQ(bedBase.getRefID("chr3"), 2);

    // get name and check if it is correct
    EXPECT_EQ(bedBase.getChromosomeName(0), "chr1");
    EXPECT_EQ(bedBase.getChromosomeName(1), "chr5");
    EXPECT_EQ(bedBase.getChromosomeName(2), "chr3");
}

TEST(TBed_baseTest, addChromosomes_throw){
    genometools::TBed_base bedBase;

    // create TChromosomes
    genometools::TChromosomes chrs;
    chrs.appendChromosome("chr1", 10);
    chrs.appendChromosome("chr5", 10);
    chrs.appendChromosome("chr3", 10);
    chrs.appendChromosome("chr3", 10); // add chr3 twice -> different refID's

    // add it to bedBase
    EXPECT_THROW(bedBase.addChromosomes(chrs), coretools::err::TDevError);
}

TEST(TBed_baseTest, addChromosomes_sameRefID_dontthrow){
    genometools::TBed_base bedBase;

    // first add single chromosome
    auto it = bedBase.addChromosome("chr1");
    EXPECT_EQ(it->name, "chr1");
    EXPECT_EQ(it->refId, 0);

    // create TChromosomes
    genometools::TChromosomes chrs;
    chrs.appendChromosome("chr1", 10);
    chrs.appendChromosome("chr5", 10);
    chrs.appendChromosome("chr3", 10);

    // add it to bedBase -> chr1 has same refID -> shouldn't throw
    bedBase.addChromosomes(chrs);

    // get RefID and check if it is correct
    EXPECT_EQ(bedBase.getRefID("chr1"), 0);
    EXPECT_EQ(bedBase.getRefID("chr5"), 1);
    EXPECT_EQ(bedBase.getRefID("chr3"), 2);

    // get name and check if it is correct
    EXPECT_EQ(bedBase.getChromosomeName(0), "chr1");
    EXPECT_EQ(bedBase.getChromosomeName(1), "chr5");
    EXPECT_EQ(bedBase.getChromosomeName(2), "chr3");
}

TEST(TBed_baseTest, getters){
    genometools::TBed_base bedBase;

    // add chromosomes
    auto it = bedBase.addChromosome("chr1");
    it = bedBase.addChromosome("chr5");
    it = bedBase.addChromosome("chr3");

    // get RefID
    EXPECT_EQ(bedBase.getRefID("chr1"), 0);
    EXPECT_EQ(bedBase.getRefID("chr5"), 1);
    EXPECT_EQ(bedBase.getRefID("chr3"), 2);

    // get name
    EXPECT_EQ(bedBase.getChromosomeName(0), "chr1");
    EXPECT_EQ(bedBase.getChromosomeName(1), "chr5");
    EXPECT_EQ(bedBase.getChromosomeName(2), "chr3");

    // has windows on chromosomes
    EXPECT_TRUE(bedBase.hasWindowsOnChr("chr1"));
    EXPECT_FALSE(bedBase.hasWindowsOnChr("chr10"));
    EXPECT_TRUE(bedBase.hasWindowsOnChr(0));
    EXPECT_FALSE(bedBase.hasWindowsOnChr(10));

    // invalid requests
    EXPECT_THROW(bedBase.getRefID("1"), coretools::err::TDevError);
    EXPECT_THROW(bedBase.getChromosomeName(5), coretools::err::TDevError);
}
*/

//-------------------------------------------------------------
// TBed
//-------------------------------------------------------------

TEST(TBedTest, add_TGenomeWindow_otherChr){
    // create bed and add chromosome
    genometools::TBed bed;
    bed.addChromosome("chr1");

    // create window which has a different chromosome and add -> throws
    EXPECT_THROW(bed.add(genometools::TGenomeWindow(1, 10, 20)), coretools::err::TDevError);
}

TEST(TBedTest, add_TGenomeWindow_overlapsOneDownstream){
    // create bed and add chromosome
    genometools::TBed bed;
    bed.addChromosome("chr1");

    //        |----------------------|
    //        10                     20
    //                     |---------------------------------|
    //                     15                                30

    bed.add(genometools::TGenomeWindow(0, 10, 20));
    bed.add(genometools::TGenomeWindow(0, 15, 30));

    // check if windows were correctly merged
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 10, 30)));
}

TEST(TBedTest, add_TGenomeWindow_overlapsOneUpstream){
    // create bed and add chromosome
    genometools::TBed bed;
    bed.addChromosome("chr1");

    //                    |----------------------|
    //                    10                     20
    //         |----------------------|
    //         5                     15

    bed.add(genometools::TGenomeWindow(0, 10, 20));
    bed.add(genometools::TGenomeWindow(0, 5, 15));

    // check if windows were correctly merged
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 5, 20)));
}

TEST(TBedTest, add_TGenomeWindow_overlapsOneExactly){
    // create bed and add chromosome
    genometools::TBed bed;
    bed.addChromosome("chr1");

    //                    |----------------------|
    //                    10                     20
    //                    |----------------------|
    //                    10                     20
    bed.add(genometools::TGenomeWindow(0, 10, 20));
    bed.add(genometools::TGenomeWindow(0, 10, 20));

    // check if windows were correctly merged
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 10, 20)));
}

TEST(TBedTest, add_TGenomeWindow_overlapsOneOverhangs){
    // create bed and add chromosome
    genometools::TBed bed;
    bed.addChromosome("chr1");

    //                    |----------------------|
    //                    10                     20
    //          |--------------------------------------------|
    //          5                                            25

    bed.add(genometools::TGenomeWindow(0, 10, 20));
    bed.add(genometools::TGenomeWindow(0, 5, 25));

    // check if windows were correctly merged
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 5, 25)));
}

TEST(TBedTest, add_TGenomeWindow_noOverlap){
    // create bed and add chromosome
    genometools::TBed bed;
    bed.addChromosome("chr1");

    //                    |----------------------|
    //                    10                     20
    //                                                                 |----------------------|
    //                                                                 30                     40

    bed.add(genometools::TGenomeWindow(0, 10, 20));
    bed.add(genometools::TGenomeWindow(0, 30, 40));

    // check if windows were not merged
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 10, 20)));
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 30, 40)));
}

TEST(TBedTest, add_TGenomeWindow_overlapsTwoDownstream){
    // create bed and add chromosome
    genometools::TBed bed;
    bed.addChromosome("chr1");

    //                    |----------------------|                           |----------------------|
    //                    10                     20                          30                     40
    //                                |------------------------------------------------------------------------------|
    //                                15                                                                             50

    bed.add(genometools::TGenomeWindow(0, 10, 20));
    bed.add(genometools::TGenomeWindow(0, 30, 40));
    bed.add(genometools::TGenomeWindow(0, 15, 50));

    // check if windows were merged
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 10, 50)));
}

TEST(TBedTest, add_TGenomeWindow_overlapsTwoUpstream){
    // create bed and add chromosome
    genometools::TBed bed;
    bed.addChromosome("chr1");

    //                    |----------------------|                           |----------------------|
    //                    10                     20                          30                     40
    //      |------------------------------------------------------------------------------|
    //      5                                                                              35

    bed.add(genometools::TGenomeWindow(0, 10, 20));
    bed.add(genometools::TGenomeWindow(0, 30, 40));
    bed.add(genometools::TGenomeWindow(0, 5, 35));

    // check if windows were merged
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 5, 40)));
}

TEST(TBedTest, add_TGenomeWindow_overlapsTwoExactly){
    // create bed and add chromosome
    genometools::TBed bed;
    bed.addChromosome("chr1");

    //                    |----------------------|                           |----------------------|
    //                    10                     20                          30                     40
    //                    |-------------------------------------------------------------------------|
    //                    10                                                                        40

    bed.add(genometools::TGenomeWindow(0, 10, 20));
    bed.add(genometools::TGenomeWindow(0, 30, 40));
    bed.add(genometools::TGenomeWindow(0, 10, 40));

    // check if windows were merged
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 10, 40)));
}

TEST(TBedTest, add_TGenomeWindow_overlapsTwoOverhang){
    // create bed and add chromosome
    genometools::TBed bed;
    bed.addChromosome("chr1");

    //                    |----------------------|                           |----------------------|
    //                    10                     20                          30                     40
    //         |-----------------------------------------------------------------------------------------------------|
    //         5                                                                                                     50

    bed.add(genometools::TGenomeWindow(0, 10, 20));
    bed.add(genometools::TGenomeWindow(0, 30, 40));
    bed.add(genometools::TGenomeWindow(0, 5, 50));

    // check if windows were merged
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 5, 50)));
}

TEST(TBedTest, add_TGenomeWindow_overlapsThreeDownstream){
    // create bed and add chromosome
    genometools::TBed bed;
    bed.addChromosome("chr1");

    //                    |----------------------|                           |----------------------|                           |----------------------|
    //                    10                     20                          30                     40                          50                     60
    //                                |--------------------------------------------------------------------------------------------------------------------------|
    //                                15                                                                                                                         65

    bed.add(genometools::TGenomeWindow(0, 10, 20));
    bed.add(genometools::TGenomeWindow(0, 30, 40));
    bed.add(genometools::TGenomeWindow(0, 50, 60));
    bed.add(genometools::TGenomeWindow(0, 15, 65));

    // check if windows were merged
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 10, 65)));
}

TEST(TBedTest, add_TGenomeWindow_overlapsThreeUpstream){
    // create bed and add chromosome
    genometools::TBed bed;
    bed.addChromosome("chr1");

    //                    |----------------------|                           |----------------------|                           |----------------------|
    //                    10                     20                          30                     40                          50                     60
    //          |---------------------------------------------------------------------------------------------------------------------------|
    //          5                                                                                                                           55

    bed.add(genometools::TGenomeWindow(0, 10, 20));
    bed.add(genometools::TGenomeWindow(0, 30, 40));
    bed.add(genometools::TGenomeWindow(0, 50, 60));
    bed.add(genometools::TGenomeWindow(0, 5, 55));

    // check if windows were merged
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 5, 60)));
}

TEST(TBedTest, add_TGenomeWindow_twoChromosomes){
    // create bed and add chromosomes
    genometools::TBed bed;
    bed.addChromosome("chr1");
    bed.addChromosome("chr2");

    //                    |----------------------|
    //                    chr1:10                chr1:20
    //         |----------------------|
    //         chr2:5                chr2:15

    bed.add(genometools::TGenomeWindow(0, 10, 20));
    bed.add(genometools::TGenomeWindow(1, 5, 15));

    // check if windows were NOT merged
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 10, 20)));
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(1, 5, 15)));
}

TEST(TBedTest, add_TGenomePosition){
    // create bed and add chromosome
    genometools::TBed bed;
    bed.addChromosome("chr1");

    bed.add(genometools::TGenomeWindow(0, 10, 20));
    bed.add(genometools::TGenomeWindow(0, 30, 40));
    bed.add(genometools::TGenomeWindow(0, 50, 60));
    bed.add(genometools::TGenomePosition(0, 5));
    bed.add(genometools::TGenomePosition(0, 20));
    bed.add(genometools::TGenomePosition(0, 41));

    // check if windows were merged
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 5, 6)));
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 10, 21)));
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 30, 40)));
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 41, 42)));
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 50, 60)));
}

TEST(TBedTest, add_RefID_pos){
    // create bed and add chromosome
    genometools::TBed bed;
    bed.addChromosome("chr1");

    bed.add(genometools::TGenomeWindow(0, 10, 20));
    bed.add(genometools::TGenomeWindow(0, 30, 40));
    bed.add(genometools::TGenomeWindow(0, 50, 60));
    bed.add(0, 5);
    bed.add(0, 20);
    bed.add(0, 41);
    EXPECT_THROW(bed.add(1, 10), coretools::err::TDevError);

    // check if windows were merged
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 5, 6)));
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 10, 21)));
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 30, 40)));
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 41, 42)));
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 50, 60)));
}

TEST(TBedTest, add_Chr_pos){
    // create bed and add chromosome
    genometools::TBed bed;
    bed.addChromosome("chr1");

    bed.add(genometools::TGenomeWindow(0, 10, 20));
    bed.add(genometools::TGenomeWindow(0, 30, 40));
    bed.add(genometools::TGenomeWindow(0, 50, 60));
    bed.add("chr1", 5);
    bed.add("chr1", 20);
    bed.add("chr1", 41);
    bed.add("chr2", 1);

    // check if windows were merged
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 5, 6)));
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 10, 21)));
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 30, 40)));
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 41, 42)));
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 50, 60)));
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(1, 1, 2)));
}

void writeInputBed(){
    coretools::TOutputFile out("bed.bed", 3);
    out.writeln("chr1", 10, 20);
    out.writeln("chr1", 30, 40);
    out.writeln("chr1", 50, 60);
    out.writeln("chr1", 15, 35);
    out.writeln("chr1", 1, 2);
    out.writeln("chr2", 55, 100);
    out.close();
}

TEST(TBedTest, add_Filename){
    // create bed and add chromosome
    genometools::TBed bed;

    // write input bed file
    writeInputBed();

    bed.add("bed.bed");

    // check if windows were merged
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 1, 2)));
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 10, 40)));
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 50, 60)));
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(1, 55, 100)));

    // remove file
    remove("bed.bed");
}

TEST(TBedTest, add_Filename_Chromosomes){
    // create bed and add chromosome
    genometools::TBed bed;

    // create TChromosomes
    genometools::TChromosomes chrs;
    chrs.appendChromosome("chr1", 10);
    chrs.appendChromosome("chr5", 10);
    chrs.appendChromosome("chr2", 10);

    // write input bed file
    writeInputBed();

    bed.add("bed.bed", chrs);

    // check if windows were merged
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 1, 2)));
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 10, 40)));
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(0, 50, 60)));
    EXPECT_TRUE(bed.exists(genometools::TGenomeWindow(2, 55, 100))); // refID is 2, because chr2 has refID=2 in TChromosomes chrs

    // remove file
    remove("bed.bed");
}

TEST(TBedTest, write){
    // create bed and add chromosome
    genometools::TBed bed;
    bed.addChromosome("chr1");
    bed.addChromosome("chr2");

    bed.add(genometools::TGenomeWindow(0, 10, 20));
    bed.add(genometools::TGenomeWindow(0, 30, 40));
    bed.add(genometools::TGenomeWindow(0, 50, 60));
    bed.add(genometools::TGenomeWindow(0, 15, 35));
    bed.add(genometools::TGenomeWindow(0, 1, 2));
    bed.add(genometools::TGenomeWindow(1, 55, 100));

    // write file
    bed.write("bed.bed");

    // create new bed object that reads this file as input
    genometools::TBed bed2;
    bed2.add("bed.bed");

    // check if windows were merged
    EXPECT_TRUE(bed2.exists(genometools::TGenomeWindow(0, 1, 2)));
    EXPECT_TRUE(bed2.exists(genometools::TGenomeWindow(0, 10, 40)));
    EXPECT_TRUE(bed2.exists(genometools::TGenomeWindow(0, 50, 60)));
    EXPECT_TRUE(bed2.exists(genometools::TGenomeWindow(1, 55, 100)));

    // remove file
    remove("bed.bed");
}

TEST(TBedTest, size_length){
    // create bed and add chromosome
    genometools::TBed bed;
    bed.addChromosome("chr1");
    bed.addChromosome("chr2");

    bed.add(genometools::TGenomeWindow(0, 10, 20));
    bed.add(genometools::TGenomeWindow(0, 30, 40));
    bed.add(genometools::TGenomeWindow(0, 50, 60));
    bed.add(genometools::TGenomeWindow(0, 15, 35));
    bed.add(genometools::TGenomeWindow(0, 1, 2));
    bed.add(genometools::TGenomeWindow(1, 55, 100));

    EXPECT_EQ(bed.size(), 4); // 4 windows
    EXPECT_EQ(bed.length(), 86);
}

TEST(TBedTest, numChromosomesWithWindows){
    // create bed and add chromosome
    genometools::TBed bed;
    bed.addChromosome("chr1");
    bed.addChromosome("chr2");
    bed.addChromosome("chr4");
    bed.addChromosome("chr3");

    bed.add(genometools::TGenomeWindow(0, 10, 20));
    bed.add(genometools::TGenomeWindow(0, 30, 40));
    bed.add(genometools::TGenomeWindow(0, 50, 60));
    bed.add(genometools::TGenomeWindow(0, 15, 35));
    bed.add(genometools::TGenomeWindow(0, 1, 2));
    bed.add(genometools::TGenomeWindow(1, 55, 100));
    bed.add(genometools::TGenomeWindow(3, 10, 100));

    EXPECT_EQ(bed.numChromosomesWithWindows(), 3); // chr3 does not have windows
}

TEST(TBedTest, lower_bound){
    // create bed and add chromosome
    genometools::TBed bed;
    bed.addChromosome("chr1");
    bed.addChromosome("chr2");

    bed.add(genometools::TGenomeWindow(0, 10, 20));
    bed.add(genometools::TGenomeWindow(0, 30, 40));
    bed.add(genometools::TGenomeWindow(0, 50, 60));
    bed.add(genometools::TGenomeWindow(0, 1, 2));
    bed.add(genometools::TGenomeWindow(1, 55, 100));

    auto it = bed.lower_bound(genometools::TGenomeWindow(0, 45, 48));
    EXPECT_EQ(it->refID(), 0);
    EXPECT_EQ(it->fromOnChr(), 50);
    EXPECT_EQ(it->toOnChr(), 60);

    it = bed.lower_bound(genometools::TGenomeWindow(0, 2, 15));
    EXPECT_EQ(it->refID(), 0);
    EXPECT_EQ(it->fromOnChr(), 10);
    EXPECT_EQ(it->toOnChr(), 20);

    it = bed.lower_bound(genometools::TGenomeWindow(1, 10, 70));
    EXPECT_EQ(it->refID(), 1);
    EXPECT_EQ(it->fromOnChr(), 55);
    EXPECT_EQ(it->toOnChr(), 100);
}

TEST(TBedTest, begin){
    // create bed and add chromosome
    genometools::TBed bed;
    bed.addChromosome("chr1");
    bed.addChromosome("chr2");

    bed.add(genometools::TGenomeWindow(0, 10, 20));
    bed.add(genometools::TGenomeWindow(0, 30, 40));
    bed.add(genometools::TGenomeWindow(0, 50, 60));
    bed.add(genometools::TGenomeWindow(0, 1, 2));
    bed.add(genometools::TGenomeWindow(1, 55, 100));

    auto it = bed.begin();
    EXPECT_EQ(it->refID(), 0);
    EXPECT_EQ(it->fromOnChr(), 1);
    EXPECT_EQ(it->toOnChr(), 2);
}

TEST(TBedTest, end){
    // create bed and add chromosome
    genometools::TBed bed;
    bed.addChromosome("chr1");
    bed.addChromosome("chr2");

    bed.add(genometools::TGenomeWindow(0, 10, 20));
    bed.add(genometools::TGenomeWindow(0, 30, 40));
    bed.add(genometools::TGenomeWindow(0, 50, 60));
    bed.add(genometools::TGenomeWindow(0, 1, 2));
    bed.add(genometools::TGenomeWindow(1, 55, 100));

    auto it = bed.end();
    it--;
    EXPECT_EQ(it->refID(), 1);
    EXPECT_EQ(it->fromOnChr(), 55);
    EXPECT_EQ(it->toOnChr(), 100);
}

TEST(TBedTest, begin_refID){
    // create bed and add chromosome
    genometools::TBed bed;
    bed.addChromosome("chr1");
    bed.addChromosome("chr2");

    bed.add(genometools::TGenomeWindow(0, 10, 20));
    bed.add(genometools::TGenomeWindow(0, 30, 40));
    bed.add(genometools::TGenomeWindow(0, 50, 60));
    bed.add(genometools::TGenomeWindow(0, 1, 2));
    bed.add(genometools::TGenomeWindow(1, 55, 100));

    auto it = bed.begin(0);
    EXPECT_EQ(it->refID(), 0);
    EXPECT_EQ(it->fromOnChr(), 1);
    EXPECT_EQ(it->toOnChr(), 2);

    it = bed.begin(1);
    EXPECT_EQ(it->refID(), 1);
    EXPECT_EQ(it->fromOnChr(), 55);
    EXPECT_EQ(it->toOnChr(), 100);
}

TEST(TBedTest, end_refID){
    // create bed and add chromosome
    genometools::TBed bed;
    bed.addChromosome("chr1");
    bed.addChromosome("chr2");

    bed.add(genometools::TGenomeWindow(0, 10, 20));
    bed.add(genometools::TGenomeWindow(0, 30, 40));
    bed.add(genometools::TGenomeWindow(0, 50, 60));
    bed.add(genometools::TGenomeWindow(0, 1, 2));
    bed.add(genometools::TGenomeWindow(1, 55, 100));

    auto it = bed.end(0);
    it--;
    EXPECT_EQ(it->refID(), 0);
    EXPECT_EQ(it->fromOnChr(), 50);
    EXPECT_EQ(it->toOnChr(), 60);

    it = bed.end(1);
    it--;
    EXPECT_EQ(it->refID(), 1);
    EXPECT_EQ(it->fromOnChr(), 55);
    EXPECT_EQ(it->toOnChr(), 100);
}

//-----------------------------------------------------
// TGenomeWindowList
//-----------------------------------------------------

TEST(TGenomeWindowListTest, add_TGenomeWindow_nonOverlapping){
    // create genomeWindowList
    genometools::TGenomeWindowList genomeWindowList;

    //        |----------------------|
    //        10                     20
    //                                                 |---------------------------------|
    //                                                 25                                50

    genomeWindowList.add(genometools::TGenomeWindow(0, 10, 20));
    genomeWindowList.add(genometools::TGenomeWindow(0, 25, 50));

    // check if windows were NOT merged
    EXPECT_TRUE(genomeWindowList.exists(genometools::TGenomeWindow(0, 10, 20)));
    EXPECT_TRUE(genomeWindowList.exists(genometools::TGenomeWindow(0, 25, 50)));
}

TEST(TGenomeWindowListTest, add_TGenomeWindow_overlapsOneDownstream){
    // create genomeWindowList
    genometools::TGenomeWindowList genomeWindowList;

    //        |----------------------|
    //        10                     20
    //                     |---------------------------------|
    //                     15                                30

    genomeWindowList.add(genometools::TGenomeWindow(0, 10, 20));
    genomeWindowList.add(genometools::TGenomeWindow(0, 15, 30));

    // check if windows were NOT merged
    EXPECT_TRUE(genomeWindowList.exists(genometools::TGenomeWindow(0, 10, 20)));
    EXPECT_TRUE(genomeWindowList.exists(genometools::TGenomeWindow(0, 15, 30)));
}

TEST(TGenomeWindowListTest, add_TGenomeWindow_overlapsOneUpstream){
    // create genomeWindowList
    genometools::TGenomeWindowList genomeWindowList;

    //                     |----------------------|
    //                      10                     20
    //        |---------------------------|
    //        5                          15

    genomeWindowList.add(genometools::TGenomeWindow(0, 10, 20));
    genomeWindowList.add(genometools::TGenomeWindow(0, 5, 15));

    // check if windows were NOT merged
    EXPECT_TRUE(genomeWindowList.exists(genometools::TGenomeWindow(0, 10, 20)));
    EXPECT_TRUE(genomeWindowList.exists(genometools::TGenomeWindow(0, 5, 15)));
}

TEST(TGenomeWindowListTest, add_TGenomeWindow_differentChr){
    // create genomeWindowList
    genometools::TGenomeWindowList genomeWindowList;

    //                     |----------------------|
    //                  chr1:10                chr1:20
    //        |---------------------------|
    //     chr2:5                        chr2:15

    genomeWindowList.add(genometools::TGenomeWindow(0, 10, 20));
    genomeWindowList.add(genometools::TGenomeWindow(1, 5, 15));

    // check if windows were NOT merged
    EXPECT_TRUE(genomeWindowList.exists(genometools::TGenomeWindow(0, 10, 20)));
    EXPECT_TRUE(genomeWindowList.exists(genometools::TGenomeWindow(1, 5, 15)));
}

TEST(TGenomeWindowListTest, add_Filename_Chromosomes){
    // create genomeWindowList
    genometools::TGenomeWindowList genomeWindowList;

    // create TChromosomes
    genometools::TChromosomes chrs;
    chrs.appendChromosome("chr1", 10);
    chrs.appendChromosome("chr5", 10);
    chrs.appendChromosome("chr2", 10);

    // write input bed file
    writeInputBed();

    genomeWindowList.add("bed.bed", chrs);

    // check if windows were merged
    EXPECT_TRUE(genomeWindowList.exists(genometools::TGenomeWindow(0, 1, 2)));
    EXPECT_TRUE(genomeWindowList.exists(genometools::TGenomeWindow(0, 10, 20)));
    EXPECT_TRUE(genomeWindowList.exists(genometools::TGenomeWindow(0, 30, 40)));
    EXPECT_TRUE(genomeWindowList.exists(genometools::TGenomeWindow(0, 50, 60)));
    EXPECT_TRUE(genomeWindowList.exists(genometools::TGenomeWindow(0, 15, 35)));
    EXPECT_TRUE(genomeWindowList.exists(genometools::TGenomeWindow(2, 55, 100))); // refID is 2, because chr2 has refID=2 in TChromosomes chrs

    // remove file
    remove("bed.bed");
}


TEST(TGenomeWindowListTest, write){
    // create genomeWindowList
    genometools::TGenomeWindowList genomeWindowList;

    // create TChromosomes
    genometools::TChromosomes chrs;
    chrs.appendChromosome("chr1", 10);
    chrs.appendChromosome("chr5", 10);
    chrs.appendChromosome("chr2", 10);

    genomeWindowList.add(genometools::TGenomeWindow(0, 10, 20));
    genomeWindowList.add(genometools::TGenomeWindow(0, 30, 40));
    genomeWindowList.add(genometools::TGenomeWindow(0, 50, 60));
    genomeWindowList.add(genometools::TGenomeWindow(0, 15, 35));
    genomeWindowList.add(genometools::TGenomeWindow(0, 1, 2));
    genomeWindowList.add(genometools::TGenomeWindow(1, 55, 100));

    genomeWindowList.write("bed.bed", chrs);

    // open a second genomeWindowList and read bed as input
    genometools::TGenomeWindowList genomeWindowList2;

    genomeWindowList2.add("bed.bed", chrs);

    // check if windows were merged
    EXPECT_TRUE(genomeWindowList2.exists(genometools::TGenomeWindow(0, 1, 2)));
    EXPECT_TRUE(genomeWindowList2.exists(genometools::TGenomeWindow(0, 10, 20)));
    EXPECT_TRUE(genomeWindowList2.exists(genometools::TGenomeWindow(0, 30, 40)));
    EXPECT_TRUE(genomeWindowList2.exists(genometools::TGenomeWindow(0, 50, 60)));
    EXPECT_TRUE(genomeWindowList2.exists(genometools::TGenomeWindow(0, 15, 35)));
    EXPECT_TRUE(genomeWindowList2.exists(genometools::TGenomeWindow(1, 55, 100)));

    // remove file
    remove("bed.bed");
}

TEST(TGenomeWindowListTest, size_length){
    // create genomeWindowList
    genometools::TGenomeWindowList genomeWindowList;

    genomeWindowList.add(genometools::TGenomeWindow(0, 10, 20));
    genomeWindowList.add(genometools::TGenomeWindow(0, 30, 40));
    genomeWindowList.add(genometools::TGenomeWindow(0, 50, 60));
    genomeWindowList.add(genometools::TGenomeWindow(0, 15, 35));
    genomeWindowList.add(genometools::TGenomeWindow(0, 1, 2));
    genomeWindowList.add(genometools::TGenomeWindow(1, 55, 100));

    EXPECT_EQ(genomeWindowList.size(), 6); // 6 windows
    EXPECT_EQ(genomeWindowList.length(), 96);
}

TEST(TGenomeWindowListTest, getters){
    // create genomeWindowList
    genometools::TGenomeWindowList genomeWindowList;

    genomeWindowList.add(genometools::TGenomeWindow(0, 10, 20));
    genomeWindowList.add(genometools::TGenomeWindow(0, 30, 40));
    genomeWindowList.add(genometools::TGenomeWindow(0, 50, 60));
    genomeWindowList.add(genometools::TGenomeWindow(0, 15, 35));
    genomeWindowList.add(genometools::TGenomeWindow(0, 1, 2));
    genomeWindowList.add(genometools::TGenomeWindow(1, 55, 100));
    genomeWindowList.add(genometools::TGenomeWindow(3, 10, 100));

    // numChromosomesWithWindows
    EXPECT_EQ(genomeWindowList.numChromosomesWithWindows(), 3);

    // hasWindowsOnChr
    EXPECT_TRUE(genomeWindowList.hasWindowsOnChr(0));
    EXPECT_TRUE(genomeWindowList.hasWindowsOnChr(1));
    EXPECT_TRUE(genomeWindowList.hasWindowsOnChr(3));
    EXPECT_FALSE(genomeWindowList.hasWindowsOnChr(2));

    // numWindowsOnChr
    EXPECT_EQ(genomeWindowList.numWindowsOnChr(0), 5);
    EXPECT_EQ(genomeWindowList.numWindowsOnChr(1), 1);
    EXPECT_EQ(genomeWindowList.numWindowsOnChr(3), 1);
    EXPECT_EQ(genomeWindowList.numWindowsOnChr(2), 0);
}

TEST(TGenomeWindowListTest, empty){
    // create genomeWindowList
    genometools::TGenomeWindowList genomeWindowList;
    EXPECT_TRUE(genomeWindowList.empty());

    genomeWindowList.add(genometools::TGenomeWindow(0, 10, 20));
    EXPECT_FALSE(genomeWindowList.empty());
}

TEST(TGenomeWindowListTest, begin){
    // create bed and add chromosome
    genometools::TGenomeWindowList genomeWindowList;

    genomeWindowList.add(genometools::TGenomeWindow(3, 10, 100));
    genomeWindowList.add(genometools::TGenomeWindow(0, 10, 20));
    genomeWindowList.add(genometools::TGenomeWindow(0, 30, 40));
    genomeWindowList.add(genometools::TGenomeWindow(0, 50, 60));
    genomeWindowList.add(genometools::TGenomeWindow(1, 55, 100));
    genomeWindowList.add(genometools::TGenomeWindow(0, 15, 35));
    genomeWindowList.add(genometools::TGenomeWindow(0, 1, 2));

    auto it = genomeWindowList.begin();
    EXPECT_EQ(it->refID(), 0);
    EXPECT_EQ(it->fromOnChr(), 1);
    EXPECT_EQ(it->toOnChr(), 2);
}

TEST(TGenomeWindowListTest, end){
    // create bed and add chromosome
    genometools::TGenomeWindowList genomeWindowList;

    genomeWindowList.add(genometools::TGenomeWindow(3, 10, 100));
    genomeWindowList.add(genometools::TGenomeWindow(0, 10, 20));
    genomeWindowList.add(genometools::TGenomeWindow(0, 30, 40));
    genomeWindowList.add(genometools::TGenomeWindow(0, 50, 60));
    genomeWindowList.add(genometools::TGenomeWindow(1, 55, 100));
    genomeWindowList.add(genometools::TGenomeWindow(0, 15, 35));
    genomeWindowList.add(genometools::TGenomeWindow(0, 1, 2));

    auto it = genomeWindowList.end();
    it--;
    EXPECT_EQ(it->refID(), 3);
    EXPECT_EQ(it->fromOnChr(), 10);
    EXPECT_EQ(it->toOnChr(), 100);
}


TEST(TGenomeWindowListTest, loop){
    // create bed and add chromosome
    genometools::TGenomeWindowList genomeWindowList;

    genomeWindowList.add(genometools::TGenomeWindow(3, 10, 100));
    genomeWindowList.add(genometools::TGenomeWindow(0, 10, 20));
    genomeWindowList.add(genometools::TGenomeWindow(1, 55, 100));
    genomeWindowList.add(genometools::TGenomeWindow(0, 1, 2));

    int i = 0;
    for (auto it = genomeWindowList.begin(); it != genomeWindowList.end(); it++, i++){
        // just check if entries are correct (different for each loop)
        if (i == 0){
            EXPECT_EQ(it->refID(), 0);
            EXPECT_EQ(it->fromOnChr(), 1);
            EXPECT_EQ(it->toOnChr(), 2);
        } else if (i == 1){
            EXPECT_EQ(it->refID(), 0);
            EXPECT_EQ(it->fromOnChr(), 10);
            EXPECT_EQ(it->toOnChr(), 20);
        } else if (i == 2){
            EXPECT_EQ(it->refID(), 1);
            EXPECT_EQ(it->fromOnChr(), 55);
            EXPECT_EQ(it->toOnChr(), 100);
        } else {
            EXPECT_EQ(it->refID(), 3);
            EXPECT_EQ(it->fromOnChr(), 10);
            EXPECT_EQ(it->toOnChr(), 100);
        }
    }
}


