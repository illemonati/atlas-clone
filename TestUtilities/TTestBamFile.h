/*
 * TTestBamFile.h
 *
 *  Created on: Aug 10, 2020
 *      Author: wegmannd
 */

#ifndef TESTUTILITIES_TTESTBAMFILE_H_
#define TESTUTILITIES_TTESTBAMFILE_H_

#include "TBamFile.h"
#include "globalConstants.h"
#include "stringFunctions.h"

namespace TestUtilities{

//--------------------------------------
// TTestBamFile
//--------------------------------------
class TTestBamFile{
protected:
	//header info
	BAM::TSamHeader _header;
	BAM::TChromosomes _chromosomes;
	BAM::TReadGroups _readGroups;

	//BAM file for writing
	std::string _filename;
	BAM::TOutputBamFile _bamFile;

	//tmp vars for dummy alignments
	uint32_t _counter;
	std::string _dummySequence;
	std::string _dummyQualities;
	uint16_t _dummyMapQual;
	uint32_t _dummyBasePos, _dummyQualPos;
	uint32_t _dummyLength, _dummyMinLength, _dummyMaxLength;
	bool _dummyIsReverseStrand;
	uint32_t _dummyReadGroup;
	std::string _dummyCigarChars; uint32_t _dummyCigarPos;
	BAM::TSamFlags _dummyFlag;

    std::string _constructSequence(uint32_t length);
    std::string _constructSequenceQualities(uint32_t length);

    void _iterateReadGroupAndReverseStrand();
    void _iterateCigar(BAM::TCigar & cigar, uint32_t length);
    void _iterateLength();
    virtual void _iterateFlags();
    void _iterateMappingQuality();

    BAM::TAlignment _constructAlignment(const std::string & sequence, const std::string & qualities, const BAM::TGenomePosition & position, const BAM::TCigar & cigar, const uint32_t & readGroup, const bool & isReverseStrand, const bool & complicatedSamFlag);

        //other
	BAM::TQualityMap _qualMap;
	GenotypeLikelihoods::TGenotypeMap _genoMap;

	//storage of written alignments
	std::vector<BAM::TAlignment> _writtenAlignments;
    void _storeAlignment(const BAM::TAlignment & alignment);

    void _initialize(const std::vector<uint32_t> ChrLength, const uint32_t & NumReadGroups);
	void _initializeChromosomes(const std::vector<uint32_t> ChrLength);
	void _initializeReadGroups(const uint32_t & NumReadGroups);
    uint32_t _computeDistanceBetweenAlignments(const uint32_t & numAlignments);

public:
	TTestBamFile(const std::vector<uint32_t> ChrLength, const uint32_t & NumReadGroups);
	TTestBamFile(const std::string & Filename, const std::vector<uint32_t> ChrLength, const uint32_t & NumReadGroups);

	void openOutput(const std::string & Filename);
	void closeOutput();
	void writeAlignment(const BAM::TAlignment & alignment);
	// write dummy alignments where sequence and qualities are shuffled all the time
	virtual void writeDummyAlignment(const BAM::TGenomePosition & position, const BAM::TCigar & cigar, const uint32_t & readGroup, const bool & isReverseStrand, const bool & complicatedSamFlag = false);
	void writeDummyAlignment(const BAM::TGenomePosition & position, const BAM::TCigar & cigar, const bool & complicatedSamFlag = false);
	void writeDummyAlignment(const BAM::TGenomePosition & position, const uint32_t & length, const bool & complicatedSamFlag = false);
	void writeDummyAlignment(const BAM::TGenomePosition & position, const bool & complicatedSamFlag = false);
	virtual void writeDummyAlignments(const uint32_t & numAlignments, const bool & complicatedSamFlag = false);

    // write dummy alignments where sequence and qualities are same within one alignment
    void writeDummyAlignment(const char& oneBase, const char& oneQual, const BAM::TGenomePosition & position, const BAM::TCigar & cigar, const uint32_t & readGroup, const bool & isReverseStrand);
    void writeDummyAlignment(const char& oneBase, const char& oneQual, const BAM::TGenomePosition & position, const BAM::TCigar & cigar);
    void writeDummyAlignment(const char& oneBase, const char& oneQual, const BAM::TGenomePosition & position, const uint32_t & length);
    void writeDummyAlignment(const char& oneBase, const char& oneQual, const BAM::TGenomePosition & position);

    //getters
	std::string filename()const { return _filename; };
	BAM::TSamHeader& header(){ return _header; };
	BAM::TChromosomes& chromosomes(){ return _chromosomes; };
	BAM::TReadGroups& readGroups(){ return _readGroups; };
	std::vector<BAM::TAlignment>::iterator beginWrittenAlignments(){ return _writtenAlignments.begin(); };
	std::vector<BAM::TAlignment>::iterator endWrittenAlignments(){ return _writtenAlignments.end(); };
};

//--------------------------------------
// TTestBamFilePairedEnd
//--------------------------------------

class TTestBamFilePairedEnd : public TTestBamFile {
protected:
    void _iterateFlags() override;

    BAM::TAlignment & _pickFirstMate(std::vector<bool> & used);
    BAM::TAlignment & _pickSecondMate(uint32_t refIDMate1, std::vector<bool> & used);

public:
    TTestBamFilePairedEnd(const std::vector<uint32_t> ChrLength, const uint32_t & NumReadGroups);
    TTestBamFilePairedEnd(const std::string & Filename, const std::vector<uint32_t> ChrLength, const uint32_t & NumReadGroups);

    void writeDummyAlignment(const BAM::TGenomePosition & position, const BAM::TCigar & cigar, const uint32_t & readGroup, const bool & isReverseStrand, const bool & complicatedSamFlag = false) override;
    void writeDummyAlignments(const uint32_t & numAlignments, const bool & complicatedSamFlag) override;
};

}; //end namespace

#endif /* TESTUTILITIES_TTESTBAMFILE_H_ */
