/*
 * TBamFile.h
 *
 *  Created on: May 18, 2020
 *      Author: wegmannd
 */

#ifndef BAM_TBAMFILE_H_
#define BAM_TBAMFILE_H_

#include "../bamtools/api/BamReader.h"
#include "../bamtools/api/BamWriter.h"
#include "TBamFilter.h"
#include "TSamHeader.h"
#include "TAlignment.h"
#include "TNumericRange.h"
#include "globalConstants.h"

namespace BAM{

//-----------------------------------------------------
//TBamFile
//-----------------------------------------------------
class TOutputBamFile; //forward declaration

class TBamFile{
	friend class TOutputBamFile;
private:
	//BAM file
	std::string _filename;
	BamTools::BamReader _bamReader;
	BamTools::BamRegion _bamRegion;
 	BamTools::SamHeader _bamHeader;
 	bool _open;
 	int64_t _fileSize;

 	//header
 	TChromosomes _chromosomes;
 	std::vector<TChromosome>::iterator _curChromosome;
 	TReadGroups _readGroups;
 	TSamHeader _samHeader;

 	//counters
 	uint64_t _numAlignmentRead;
 	uint64_t _numAlignmentsPassedQC;
 	bool _limitNumReads;
 	uint64_t _maxNumReadsToRead;

 	//current alignment
 	BamTools::BamAlignment _curBamAlignment;
 	uint16_t _curReadGroupID;
 	TCigar _curCigar;
 	TGenomePosition _curAlignmentPosition, _previousAlignmentPosition;
 	bool _chrChanged;

	//alignment filters
 	bool _QCFiltersPassed;
 	uint16_t _maxReadLength;
 	TAlignmentList _blacklist;
 	TBamFileFilterRange<uint32_t> _readLengthFilter;
 	bool _allowTooLongReads;
 	bool _keepAll;
 	TBamFileFilterBool _duplicateFilter;
 	TBamFileFilterBool _softClippedFilter;
 	TBamFileFilterBool _improperPairsFilter;
 	TBamFileFilterBool _unmappedFilter;
 	TBamFileFilterBool _failedQCFilter;
 	TBamFileFilterBool _secondaryFilter;
 	TBamFileFilterBool _supplementaryFilter;
 	TBamFileFilterBool _longerThanFragmentFilter;
 	TBamFileFilterBool _readGroupFilter;
 	TBamFileFilterBool _fwdStrandFilter;
 	TBamFileFilterBool _revStrandFilter;
 	TBamFileFilterBool _firstMateFilter;
 	TBamFileFilterBool _secondMateFilter;
 	TBamFileFilterBool _blacklistFilter;
 	TBamFileFilterRange<uint8_t> _mappingQualityFilter;
 	TBamFileFilterRange<uint32_t> _fragmentLengthFilter;
 	TBamFileFilter _externalFilter;

	void _fillSamHeader(TSamHeader & SamHeader);
	void _fillChromosomes(TChromosomes & chromosomes);
	void _fillReadGroups(TReadGroups & readGroups);
 	void _applyFilters();

	//output filtered reads
	bool _updateLog;
	TBamFileLog* _bamLog;

	//report progress
	TLog* _logfile;
	coretools::TTimer _timer;
	uint32_t _progressFrequency;
	uint64_t _lastProgressPrinted;

	std::string _millionReadsRead(){ return coretools::str::to_string_with_precision((double) _numAlignmentRead / 1000000.0, 1); };
	void _openForWriting(BamTools::BamWriter & bamWriter, const std::string filename);

public:
	TBamFile();

	//access header info READ ONLY
	const TChromosomes& chromosomes() const{ return _chromosomes; };
	const TReadGroups& readGroups() const { return _readGroups; };
	const TSamHeader samHeader() const{ return _samHeader; };

	//modify header info: know what you do!
	TReadGroups& readGroupsMutable(){ return _readGroups; };

	//filters
	void setFilters(TParameters & params, TLog* logfile);
	void setLimits(TParameters & params, TLog* logfile);
	void setKeepAll();
	void curFilterOut();
	void filterOut(const std::string & alignmentName, const bool & isReverseStrand);
	void setExternalFilterReason(const std::string reason);
	void openBamLog(TParameters & params, TLog* logfile);
	void writeToBamLog(const std::string & alignmentName, const bool & isReverseStrand, const std::string & reason);

	//get filter status
 	bool duplicateFilterEnabled() const{ return _duplicateFilter.filters(); };
 	bool softClippedFilterEnabled() const{ return _softClippedFilter.filters(); };
 	bool improperPairsFilterEnabled() const{ return _improperPairsFilter.filters(); };
 	bool unmappedFilterEnabled() const{ return _unmappedFilter.filters(); };
 	bool failedQCFilterEnabled() const{ return _failedQCFilter.filters(); };
 	bool secondaryFilterEnabled() const{ return _secondaryFilter.filters(); };
 	bool supplementaryFilterEnabled() const{ return _supplementaryFilter.filters(); };
 	bool longerThanFragmentFilterEnabled() const{ return _longerThanFragmentFilter.filters(); };
 	bool readGroupFilterEnabled() const{ return _readGroupFilter.filters(); };
 	bool fwdStrandFilterEnabled() const{ return _fwdStrandFilter.filters(); };
 	bool revStrandFilterEnabled() const{ return _revStrandFilter.filters(); };
 	bool firstMateFilterEnabled() const{ return _firstMateFilter.filters(); };
 	bool secondMateFilterEnabled() const{ return _secondMateFilter.filters(); };
 	bool blacklistFilterEnabled() const{ return _blacklistFilter.filters(); };
 	bool mappingQualityFilterEnabled() const{ return _mappingQualityFilter.filters(); };
 	bool fragmentLengthfilterEnabled() const{ return _fragmentLengthFilter.filters(); };
 	bool externalFilterEnabled() const{ return _externalFilter.filters(); };

	//reading
	void open(const std::string Filename, const bool IndexNotRequired, TLog* Logfile);
	bool isOpen() const{ return _open; };
	void close();
	bool readNextAlignment();
	bool readNextAlignmentThatPassesFilters();
	bool readNextAlignment(TAlignment & alignment);
	bool readNextAlignmentThatPassesFilters(TAlignment & alignment);
	void fill(TAlignment & alignment) const;

	bool jump(const TGenomePosition Position);
	void rewind();

	//writing
	void writeCurAlignment(TOutputBamFile & out);

	//getters for cur alignment
	const std::string curName() const{ return _curBamAlignment.Name; };
	TGenomePosition curPosition() const { return _curAlignmentPosition; };
	uint32_t refID() const{ return _curChromosome->refID(); };
	const TChromosome& curChromosome() const{ return *_curChromosome; };
	const TCigar& curCIGAR() const{ return _curCigar; };
	uint16_t curReadGroupID() const{ return _curReadGroupID; };
	bool chrChanged() const{ return _chrChanged; };
	bool curPassedQC() const{ return _QCFiltersPassed; };
	uint16_t curReadLength() const{ return _curCigar.lengthRead(); };
	uint16_t curUsableSequence() const{ return _curCigar.lengthSequenced(); };
	uint16_t curFragmentLength() const;
	uint16_t curUsableAlignedLength(TQualityFilter & qualFilter) const;
	uint16_t curMappingQuality() const{ return _curBamAlignment.MapQuality; };
	bool curIsPaired() const{ return _curBamAlignment.IsPaired(); };
	bool curIsProperPair() const{ return _curBamAlignment.IsProperPair(); };
	bool curIsReverseStrand() const{ return _curBamAlignment.IsReverseStrand(); };
	bool curIsDuplicate() const{ return _curBamAlignment.IsDuplicate(); };
	bool curIsMapped() const{ return _curBamAlignment.IsMapped(); };
	bool curIsFailedQC() const { return _curBamAlignment.IsFailedQC(); };
	bool curIsSecondary() const{ return !_curBamAlignment.IsPrimaryAlignment(); };
	bool curIsSupplementary() const{ return _curBamAlignment.IsSupplementary(); };
	bool curIsLongerThanFragment() const {return _curBamAlignment.IsProperPair() && _curBamAlignment.InsertSize < _curCigar.lengthAligned(); };
	bool curIsFirstMate() const{ return _curBamAlignment.IsFirstMate(); };
	bool curIsSecondMate() const{ return _curBamAlignment.IsSecondMate(); };
	std::string curQuerySequence(const uint16_t start, const uint16_t length) const;

	//modify cur alignment
	void curSetNewReadGroup(const uint16_t id);
	void curAddSamField();
	void curAddSamField(const std::string tag, const std::string value);
	void curAddSamField(const std::string tag, const float value);

	//other getters
	bool isOpen(){ return _open; };
	std::string filename() const{ return _filename; };
	uint16_t maxReadLength(){ return _maxReadLength; };
	uint64_t numAlignmentsRead(){ return _numAlignmentRead; };
	double positionInFile(){ return (double) _bamReader.tell() / (double) _fileSize; };
	uint16_t numReadGroups() const{ return _readGroups.size(); };

	//progress reporting
	//TODO: try to make general and include in common utilities
	void printSummaryNoEndIndent();
	void printSummary();
	void startProgressReporting(uint32_t Frequency=1000000);
	void printProgress();
	void printEndWithSummary();
	void printEndNoEndIndent();

	//comparisons
	bool operator==(const TGenomePosition & Position) const{ return _curAlignmentPosition == Position; };
	bool operator<(const TGenomePosition & Position) const{ return _curAlignmentPosition < Position; };
	bool operator>(const TGenomePosition & Position) const{ return _curAlignmentPosition > Position; };
	bool operator<(const TGenomeWindow & Window) const{ return _curAlignmentPosition < Window; };
	bool operator>(const TGenomeWindow & Window) const{ return _curAlignmentPosition > Window; };
	bool operator<(const TChromosome & Chromosome) const{ return _curAlignmentPosition < Chromosome.chrStart; };
	bool operator>(const TChromosome & Chromosome) const{ return _curAlignmentPosition > Chromosome.chrEnd; };
};

//------------------------------------------------
// TQualityAdjusterForWriting
// Manages the printing of quality scores when writing BAM files
//------------------------------------------------
class TQualityAdjusterForWriting{
private:
	bool _adjust;
	bool _binIllumina;
	bool _limitRange;
	BaseQuality _minQual {BaseQuality::min()};
	BaseQuality _maxQual {BaseQuality::max()};

	char _adjustOneQuality(BaseQuality qual) const;

public:
	TQualityAdjusterForWriting();

	bool adjusts() const { return _adjust; };
	void binQualitiesIllumina();
	void limitRange(const BaseQuality & min, const BaseQuality & max);
	void limitRange(const TNumericRange<uint8_t> & Range);
	std::string rangeString();
	void adjustQualities(std::string & qualities) const;
};

//----------------------------------------------------
//TOutputBamFile
//----------------------------------------------------
class TOutputBamFile{
	friend TBamFile;

private:
 	std::string _outputFilename;
 	BamTools::BamWriter _bamWriter;
 	bool _openForWriting;
 	const TReadGroups* _readGroups;

 	std::multiset<TAlignment, std::less<>> _futureAlignments;

 	//quality output transformations
 	TQualityAdjusterForWriting _qualityAdjuster;

 	void _writeAlignment(const TAlignment & alignment);
 	void _writeAlignment(BamTools::BamAlignment & alignment);

public:
 	TOutputBamFile();
 	TOutputBamFile(const std::string filename, const TBamFile & original);
 	~TOutputBamFile();

 	void open(const std::string Filename, const TSamHeader & Header, const TChromosomes & Chromosomes, const TReadGroups & ReadGroups);
	void open(const std::string Filename, const TBamFile & Original);
	void open(TParameters & params, TLog* logfile, const std::string Filename, const TSamHeader & Header, const TChromosomes & Chromosomes, const TReadGroups & ReadGroups);
	void setQualityAdjusterForWriting(TParameters & params, TLog* logfile);
	bool isOpen() const{ return _openForWriting; };
	void close(TLog* logfile);
	void close();
	void closeNoIndex();
	void writeAlignment(const TAlignment & alignment);
};

}; //end namespace


#endif /* BAM_TBAMFILE_H_ */
