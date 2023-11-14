/*
 * TBamFile.h
 *
 *  Created on: May 18, 2020
 *      Author: wegmannd
 */

#ifndef BAM_TBAMFILE_H_
#define BAM_TBAMFILE_H_

#include <cstdint>
#include <functional>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "api/BamAlignment.h"
#include "api/BamAux.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/SamHeader.h"

#include "coretools/Math/TNumericRange.h"
#include "coretools/Math/counters.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/TTimer.h"
#include "genometools/GenomePositions/TChromosomes.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "genometools/PhredProbabilityTypes.h"

#include "TAlignment.h"
#include "TAlignmentList.h"
#include "TBamFilter.h"
#include "TCigar.h"
#include "TReadGroups.h"
#include "TSamHeader.h"

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
	size_t _stepSizeFindLastAlignment = 100;
 	size_t _fileSize;

 	//header
 	genometools::TChromosomes _chromosomes;
 	std::vector<genometools::TChromosome>::iterator _curChromosome;
 	TReadGroups _readGroups;
 	TSamHeader _samHeader;

 	//counters
 	size_t _numAlignmentRead;
 	coretools::TCountDistributionVector<> _numAlignmentReadPerReadGroupPerChromosome;
 	size_t _numAlignmentsPassedQC;
 	bool _limitNumReads;
 	size_t _maxNumReadsToRead;

 	//current alignment
 	BamTools::BamAlignment _curBamAlignment;
 	size_t _curReadGroupID;
 	TCigar _curCigar;
 	genometools::TGenomePosition _curAlignmentPosition, _previousAlignmentPosition;
	bool _chrChanged;
	double _softClipFilterRatio = 1;

	//alignment filters
	bool _QCFiltersPassed;
 	TAlignmentList _blacklist;
 	bool _allowTooLongReads;
 	bool _keepAll; 	
	size_t _numNoReadGroup;
	std::vector<size_t> _numNotAligned;
	TBamFilterRange _readLengthFilter{500};
 	TBamFilterRange _mappedLengthFilter{500};
 	TBamFilterBool _duplicateFilter;
 	TBamFilterBool _softClippedRatioFilter;
 	TBamFilterBool _improperPairsFilter;
 	TBamFilterBool _unmappedFilter;
 	TBamFilterBool _failedQCFilter;
 	TBamFilterBool _secondaryFilter;
 	TBamFilterBool _supplementaryFilter;
 	TBamFilterBool _longerThanFragmentFilter;
 	TBamFilterBool _readGroupFilter;
 	TBamFilterBool _fwdStrandFilter;
 	TBamFilterBool _revStrandFilter;
 	TBamFilterBool _firstMateFilter;
 	TBamFilterBool _secondMateFilter;
 	TBamFilterBool _blacklistFilter;
 	TBamFilterRange _mappingQualityFilter;
 	TBamFilterRange _fragmentLengthFilter;
 	TBamFilter _externalFilter;

	void _fillSamHeader(TSamHeader & SamHeader);
	void _fillChromosomes(genometools::TChromosomes & chromosomes);
	void _fillReadGroups(TReadGroups & readGroups);
	bool _readNextAlignmentFromFile();
 	void _applyFilters();
	void _writeFilteringStats(std::string_view outputName) const;

	//output filtered reads
	coretools::TOutputFile _bamLog;

	//report progress
	mutable coretools::TTimer _timer;
	mutable size_t _progressFrequency;
	mutable size_t _lastProgressPrinted;

	std::string _millionReadsRead() const { return coretools::str::toStringWithPrecision((double) _numAlignmentRead / 1000000.0, 1); };
	void _openForWriting(BamTools::BamWriter & bamWriter, const std::string filename);

public:
	TBamFile();
	TBamFile(std::string_view Filename) : TBamFile() {open(Filename); setLimits();}

	//access header info READ ONLY
	const genometools::TChromosomes& chromosomes() const{ return _chromosomes; };
	const TReadGroups& readGroups() const { return _readGroups; };
	const TSamHeader samHeader() const{ return _samHeader; };

	//modify header info: know what you do!
	TReadGroups& readGroupsMutable(){ return _readGroups; };

	//filters
	void setFilters();
	void setLimits();
	void curFilterOut();
	void filterOut(const TAlignment & Alignment);
	void setExternalFilterReason(std::string_view reason);
	void openBamLog();
	void writeToBamLog(std::string_view alignmentName, bool isReverseStrand, std::string_view reason);

	//get filter status
 	bool duplicateFilterEnabled() const{ return _duplicateFilter.filters(); };
 	bool softClippedFilterEnabled() const{ return _softClippedRatioFilter.filters(); };
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
	void open(std::string_view Filename);
	bool isOpen() const{ return _bamReader.IsOpen(); };
	void close();
	bool readNextAlignment();
	bool readNextAlignmentThatPassesFilters();
	bool readNextAlignment(TAlignment & alignment);
	bool readNextAlignmentThatPassesFilters(TAlignment & alignment);
	void fill(TAlignment & alignment) const;

	bool jump(const genometools::TGenomePosition Position);
	void rewind();

	//writing
	void writeCurAlignment(TOutputBamFile & out);

	//getters for cur alignment
	const std::string curName() const{ return _curBamAlignment.Name; };
	genometools::TGenomePosition curPosition() const { return _curAlignmentPosition; };
	size_t refID() const{ return _curChromosome->refID(); };
	const genometools::TChromosome& curChromosome() const{ return *_curChromosome; };
	const TCigar& curCIGAR() const{ return _curCigar; };
	size_t curReadGroupID() const{ return _curReadGroupID; };
	bool chrChanged() const{ return _chrChanged; };
	bool curPassedQC() const{ return _QCFiltersPassed; };
	size_t curReadLength() const{ return _curCigar.lengthRead(); };
	size_t curUsableSequence() const{ return _curCigar.lengthSequenced(); };
	size_t curFragmentLength() const;
	uint16_t curMappingQuality() const{ return _curBamAlignment.MapQuality; };
	bool curIsPaired() const{ return _curBamAlignment.IsPaired(); };
	bool curIsProperPair() const{ return _curBamAlignment.IsProperPair(); };
	bool curIsReverseStrand() const{ return _curBamAlignment.IsReverseStrand(); };
	bool curIsDuplicate() const{ return _curBamAlignment.IsDuplicate(); };
	bool curIsMapped() const{ return _curBamAlignment.IsMapped(); };
	bool curIsFailedQC() const { return _curBamAlignment.IsFailedQC(); };
	bool curIsSecondary() const{ return !_curBamAlignment.IsPrimaryAlignment(); };
	bool curIsSupplementary() const{ return _curBamAlignment.IsSupplementary(); };
	bool curIsLongerThanFragment() const {return _curBamAlignment.IsProperPair() && _curBamAlignment.InsertSize < static_cast<int>(_curCigar.lengthAligned()); };
	bool curIsFirstMate() const{ return _curBamAlignment.IsFirstMate(); };
	bool curIsSecondMate() const{ return _curBamAlignment.IsSecondMate(); };
	std::string curQuerySequence(const size_t start, const size_t length) const;

	//modify cur alignment
	void curSetNewReadGroup(const size_t id);
	void curAddSamField();
	void curAddSamField(const std::string& tag, const std::string& value);
	void curAddSamField(const std::string& tag, float value);

	//other getters
	std::string filename() const{ return _filename; };
	size_t maxReadLength(){ return _readLengthFilter.range().max(); };
	size_t numAlignmentsRead(){ return _numAlignmentRead; };
	coretools::TCountDistributionVector<> numAlignmentReadPerReadGroupPerChromosome() { return _numAlignmentReadPerReadGroupPerChromosome; };
	double positionInFile() const { return (double) _bamReader.Tell() / (double) _fileSize; };
	size_t numReadGroups() const{ return _readGroups.size(); };
	TBamFilterBool getDuplicateFilter(){return _duplicateFilter; };

	//progress reporting
	void printSummaryNoEndIndent(std::string_view outputName) const;
	void printSummary(std::string_view outputName) const;
	void startProgressReporting(size_t Frequency=1000000) const;
	void printProgress() const;
	void printEndWithSummary(std::string_view outputName) const;
	void printEndNoEndIndent() const;

	//comparisons
	bool operator==(const genometools::TGenomePosition & Position) const{ return _curAlignmentPosition == Position; };
	bool operator<(const genometools::TGenomePosition & Position) const{ return _curAlignmentPosition < Position; };
	bool operator>(const genometools::TGenomePosition & Position) const{ return _curAlignmentPosition > Position; };
	bool operator<(const genometools::TGenomeWindow & Window) const{ return _curAlignmentPosition < Window; };
	bool operator>(const genometools::TGenomeWindow & Window) const{ return _curAlignmentPosition > Window; };
	bool operator<(const genometools::TChromosome & Chromosome) const{ return _curAlignmentPosition < Chromosome.start(); };
	bool operator>(const genometools::TChromosome & Chromosome) const{ return _curAlignmentPosition > Chromosome.end(); };
};

//------------------------------------------------
// TQualityAdjusterForWriting
// Manages the printing of quality scores when writing BAM files
//------------------------------------------------
class TQualityAdjusterForWriting{
private:
	bool _initialized;
	bool _adjust;
	bool _binIllumina;
	bool _limitRange;
	genometools::BaseQuality _minQual {genometools::BaseQuality::min()};
	genometools::BaseQuality _maxQual {genometools::BaseQuality::max()};

	char _adjustOneQuality(genometools::BaseQuality qual) const;

public:
	TQualityAdjusterForWriting();

	void initialize();
	bool adjusts() const { return _adjust; };
	void binQualitiesIllumina();
	void limitRange(const genometools::BaseQuality & min, const genometools::BaseQuality & max);
	void limitRange(const coretools::TNumericRange<uint8_t> & Range);
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
 	TOutputBamFile(const TQualityAdjusterForWriting & QualityAdjuster);
 	TOutputBamFile(const std::string Filename, const TBamFile & Original);
 	TOutputBamFile(const std::string Filename, const TSamHeader & Header, const genometools::TChromosomes & Chromosomes, const TReadGroups & ReadGroups);
 	TOutputBamFile(const std::string Filename, const TSamHeader & Header, const genometools::TChromosomes & Chromosomes, const TReadGroups & ReadGroups, const TQualityAdjusterForWriting & QualityAdjuster);

 	~TOutputBamFile();

 	TOutputBamFile(const TOutputBamFile&) = default;
 	TOutputBamFile(TOutputBamFile&&) noexcept = default;
 	TOutputBamFile& operator=(const TOutputBamFile&) = default;
 	TOutputBamFile& operator=(TOutputBamFile&&) noexcept = default;

 	void open(const std::string Filename, const TSamHeader & Header, const genometools::TChromosomes & Chromosomes, const TReadGroups & ReadGroups);
	void open(const std::string Filename, const TBamFile & Original);
	void setQualityAdjusterForWriting();
	void setQualityAdjusterForWriting(const TQualityAdjusterForWriting & QualityAdjuster);
	bool isOpen() const{ return _openForWriting; };
	void close();
	void closeNoIndex();
	void writeAlignment(const TAlignment & alignment);
	void writeAlignmentLater(const TAlignment & alignment);
};

}; //end namespace


#endif /* BAM_TBAMFILE_H_ */
