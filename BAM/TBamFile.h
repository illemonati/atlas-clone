/*
 * TBamFile.h
 *
 *  Created on: May 18, 2020
 *      Author: wegmannd
 */

#ifndef BAM_TBAMFILE_H_
#define BAM_TBAMFILE_H_

#include "BamData.h"
#include "bamtools/api/BamReader.h"
#include "bamtools/api/BamWriter.h"
#include "TChromosomes.h"
#include "TReadGroups.h"
#include "TAlignment.h"
#include "TBamFilter.h"

namespace BAM{

//TODO: think about splitting into a reader and writer

//-----------------------------------------------------
//TBamFile
//-----------------------------------------------------
class TBamFile{
private:
	//BAM file
	std::string _filename;
	BamTools::BamReader _bamReader;
	BamTools::BamRegion _bamRegion;
 	BamTools::SamHeader _bamHeader;
 	int64_t _fileSize;

 	//counters
 	uint64_t _numAlignmentRead;
 	uint64_t _numAlignmentsPassedQC;
 	bool _limitNumReads;
 	uint64_t _maxNumReadsToRead;

 	//current alignment
 	BamTools::BamAlignment _curBamAlignment;
 	uint16_t _curReadGroupID;
 	BAM::TCigar _curCigar;
 	uint32_t _previousAlignmentPos;
 	int _previousAlignmentChr; //negative at beginning to trigger chr change on first read
 	bool _chrChanged;

 	//writing
 	BamTools::BamWriter bamWriter;
 	bool _openForWriting;

	//alignment filters
 	bool _QCFiltersPassed;
 	uint16_t _maxReadLength;
 	TAlignmentBlacklist _blacklist;
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
 	TBamFileFilterRange _mappingQualityFilter;
 	TBamFileFilterRange _fragmentLengthfilter;
 	TBamFileFilter _externalFilter;

	void _limitChromosomes(TParameters & params, TLog* logfile);
	void _limitReadGroups(TParameters & params, TLog* logfile);
	void _setFilters(TParameters & params, TLog* logfile);
	void _fillChromosomes(TChromosomes & chromosomes);
	void _fillReadGroups(TReadGroups & readGroups);
 	void _applyFilters();

	//output filtered reads
	bool _updateLog;
	TBamFileLog* _bamLog;

	//report progress
	TLog* _logfile;
	TTimer _timer;
	uint32_t _progressFrequency;
	uint64_t _lastProgressPrinted;

	std::string _millionReadsRead(){ return to_string_with_precision((double) _numAlignmentRead / 1000000.0, 1); };

public:
	TChromosomes chromosomes;
	TReadGroups readGroups;

	TBamFile();

	//filters
	void setFiltersAndLimits(TParameters & params, TLog* logfile);
	void limitReadLength(const int MaxReadLength);
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
 	bool fragmentLengthfilterEnabled() const{ return _fragmentLengthfilter.filters(); };
 	bool externalFilterEnabled() const{ return _externalFilter.filters(); };

	//reading
	void open(const std::string filename, const bool indexNotRequired);
	bool readNextAlignment(); //TODO: make private
	bool readNextAlignmentThatPassesFilters(); //TODO: make private
	bool readNextAlignment(TAlignment & alignment);
	bool readNextAlignmentThatPassesFilters(TAlignment & alignment);
	void fill(TAlignment & alignment);

	bool jump(int chr, int position);
	void rewind();

	//writing
	void openOutput(std::string filename);
	void writeCurAlignment();
	void writeAlignment(TAlignment & alignment, const TGenotypeMap & genoMap, const TQualityMap & qualityMap);

	//getters for cur alignment
	const std::string curName() const{ return _curBamAlignment.Name; };
	uint32_t curPosition() const{ return _curBamAlignment.Position; };
	const TChromosome& curChromosome() const{ return chromosomes.curChromosome(); };
	const TCigar& curCIGAR() const{ return _curCigar; };
	uint16_t curReadGroupID() const{ return _curReadGroupID; };
	bool chrChanged() const{ return _chrChanged; };
	bool curPassedQC() const{ return _QCFiltersPassed; };
	uint16_t curReadLength() const{ return _curCigar.lengthRead(); };
	uint16_t curUsableSequence() const{ return _curCigar.lengthUsableSequence(); };
	uint16_t curFragmentLength() const;
	uint16_t curUsableAlignedLength(TQualityFilter & qualFilter) const;
	uint16_t curMappingQuality() const{ return _curBamAlignment.MapQuality; };
	bool curIsPaired() const{ return _curBamAlignment.IsPaired(); };
	bool curIsProperPair() const{ return _curBamAlignment.IsProperPair(); };
	bool curIsReverseStrand() const{ return _curBamAlignment.IsReverseStrand(); };

	//other getters
	std::string filename() const{ return _filename; };
	uint16_t maxReadLength(){ return _maxReadLength; };
	uint64_t numAlignmentsRead(){ return _numAlignmentRead; };
	double positionInFile(){ return (double) _bamReader.tell() / (double) _fileSize; };

	//progress reporting
	//TODO: try to make general and include in common utilities
	void printSummaryNoEndIndent(TLog* Logfile);
	void printSummary(TLog* Logfile);
	void startProgressReporting(uint32_t Frequency, TLog* Logfile);
	void printProgress();
	void printEndWithSummary();
	void printEndNoEndIndent();
};

}; //end namespace


#endif /* BAM_TBAMFILE_H_ */
