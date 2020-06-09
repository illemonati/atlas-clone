/*
 * TBamFile.h
 *
 *  Created on: May 18, 2020
 *      Author: wegmannd
 */

#ifndef BAM_TBAMFILE_H_
#define BAM_TBAMFILE_H_

#include "bamtools/api/BamReader.h"
#include "bamtools/api/BamWriter.h"

#include "TOutputBam.h"
#include "TBamFilter.h"

namespace BAM{

//-----------------------------------------------------
//TBamFile
//-----------------------------------------------------
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
 	uint32_t _previousAlignmentPos;
 	int _previousAlignmentChr; //negative at beginning to trigger chr change on first read
 	bool _chrChanged;

	//alignment filters
 	bool _QCFiltersPassed;
 	uint16_t _maxReadLength;
 	TAlignmentList _blacklist;
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

	void setFilters(TParameters & params, TLog* logfile);
	void _fillSamHeader(TSamHeader & SamHeader);
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
	void _openForWriting(BamTools::BamWriter & bamWriter, const std::string filename);

public:
	TBamFile();
	~TBamFile();

	//access header info
	const TChromosomes& chromosomes() const{ return _chromosomes; };
	const TReadGroups& readGroups() const { return _readGroups; };
	const TSamHeader& samHeader() const{ return _samHeader; };

	//filters
	void setLimits(TParameters & params, TLog* logfile);
	void setKeepAll();
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
	uint32_t curPosition() const{ return _curBamAlignment.Position; };
	const TChromosome& curChromosome() const{ return _chromosomes.curChromosome(); };
	const TCigar& curCIGAR() const{ return _curCigar; };
	uint16_t curReadGroupID() const{ return _curReadGroupID; };
	bool chrChanged() const{ return _chrChanged; };
	bool curPassedQC() const{ return _QCFiltersPassed; };
	uint16_t curReadLength() const{ return _curCigar.lengthRead(); };
	uint16_t curUsableSequence() const{ return _curCigar.lengthSequenced()(); };
	uint16_t curFragmentLength() const;
	uint16_t curUsableAlignedLength(TQualityFilter & qualFilter) const;
	uint16_t curMappingQuality() const{ return _curBamAlignment.MapQuality; };
	bool curIsPaired() const{ return _curBamAlignment.IsPaired(); };
	bool curIsProperPair() const{ return _curBamAlignment.IsProperPair(); };
	bool curIsReverseStrand() const{ return _curBamAlignment.IsReverseStrand(); };
	std::string curQuerySequence(const uint16_t start, const uint16_t length) const;

	//modify cur alignment
	void curSetNewReadGroup(const uint16_t id);
	void curFilterOut();
	void curAddSamField();
	void curAddSamField(const std::string tag, const std::string value);
	void curAddSamField(const std::string tag, const float value);

	//other getters
	bool isOpen(){ return _open; };
	std::string filename() const{ return _filename; };
	uint16_t maxReadLength(){ return _maxReadLength; };
	uint64_t numAlignmentsRead(){ return _numAlignmentRead; };
	double positionInFile(){ return (double) _bamReader.tell() / (double) _fileSize; };

	//progress reporting
	//TODO: try to make general and include in common utilities
	void printSummaryNoEndIndent();
	void printSummary();
	void startProgressReporting(uint32_t Frequency=1000000);
	void printProgress();
	void printEndWithSummary();
	void printEndNoEndIndent();
};


}; //end namespace


#endif /* BAM_TBAMFILE_H_ */
