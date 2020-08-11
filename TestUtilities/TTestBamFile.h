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
private:
	//header info
	BAM::TSamHeader _header;
	BAM::TChromosomes _chromosomes;
	BAM::TReadGroups _readGroups;

	//BAM file for writing
	std::string _filename;
	BAM::TOutputBamFile _bamFile;

	//tmp vars for dummy alignments
	std::string _dummySequence;
	std::string _dummyQualities;
	uint32_t _dummyBasePos, _dummyQualPos;
	uint32_t _dummyLength, _dummyMinLength, _dummyMaxLength;
	bool _dummyIsReverseStrand;
	uint32_t _dummyReadGroup;
	BAM::TGenomePosition _dummyPosition;
	std::string _dummyCigarChars; uint32_t _dummyCigarPos;

	//other
	BAM::TQualityMap _qualMap;
	GenotypeLikelihoods::TGenotypeMap _genoMap;

	//storage of written alignments
	std::vector<BAM::TAlignment> _writtenAlignments;

	void _initialize(const std::vector<uint32_t> ChrLength, const uint32_t & NumReadGroups);
	void _initializeChromosomes(const std::vector<uint32_t> ChrLength);
	void _initializeReadGroups(const uint32_t & NumReadGroups);

public:
	TTestBamFile(const std::vector<uint32_t> ChrLength, const uint32_t & NumReadGroups);
	TTestBamFile(const std::string & Filename, const std::vector<uint32_t> ChrLength, const uint32_t & NumReadGroups);

	void openOutput(const std::string & Filename);
	void closeOutput();
	void writeAlignment(const BAM::TAlignment & alignment);
	void writeDummyAlignment(const BAM::TGenomePosition & position, const uint32_t & length, const BAM::TCigar & cigar, const uint32_t & readGroup, const bool & isReverseStrand);
	void writeDummyAlignment(const BAM::TGenomePosition & position, const uint32_t & length, const BAM::TCigar & cigar);
    void writeDummyAlignment(const BAM::TGenomePosition & position, const uint32_t & length);
    void writeDummyAlignment(const BAM::TGenomePosition & position);
	void writeDummyAlignments(const uint32_t & numAlignments);

	//getters
	std::string filename()const { return _filename; };
	BAM::TSamHeader& header(){ return _header; };
	BAM::TChromosomes& chromosomes(){ return _chromosomes; };
	BAM::TReadGroups& readGroups(){ return _readGroups; };
	std::vector<BAM::TAlignment>::iterator beginWrittenAlignments(){ return _writtenAlignments.begin(); };
	std::vector<BAM::TAlignment>::iterator endWrittenAlignments(){ return _writtenAlignments.end(); };
};

}; //end namespace

#endif /* TESTUTILITIES_TTESTBAMFILE_H_ */
