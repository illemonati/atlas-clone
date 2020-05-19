/*
 * TBamFile.h
 *
 *  Created on: May 18, 2020
 *      Author: wegmannd
 */

#ifndef BAM_TBAMFILE_H_
#define BAM_TBAMFILE_H_

#include "bamtools/api/BamReader.h"
#include "TChromosomes.h"
#include "TReadGroups.h"
#include "TAlignment.h"
#include "TAlignmentBlacklist.h"

//-----------------------------------------------------
//TBamFile
//-----------------------------------------------------
class TBamFile{
private:
	//BAm file
	std::string _filename;
	BamTools::BamReader _bamReader;
	BamTools::BamRegion _bamRegion;
 	BamTools::SamHeader _bamHeader;
 	int64_t _fileSize;

 	//counters
 	uint64_t _numAlignmentRead;
 	uint64_t _numAlignmentsPassedQC;
 	uint32_t _previousAlignmentPos;
 	int _previousAlignmentChr; //negative at beginning to trigger chr change on first read
 	bool _chrChanged;
 	uint16_t _curReadGroupID;

	//filters
 	bool _QCFiltersPassed;
 	uint16_t _maxReadLength;
 	bool _keepAll;
 	bool _keepDuplicates;
 	bool _filterSoftClips;
 	bool _keepImproperPairs;
 	bool _keepUnmappedReads;
 	bool _keepFailedQC;
 	bool _keepSecondary;
 	bool _keepSupplementary;
 	bool _parse;

	bool applyQualityFilter;
	bool applyContextFilter;
	std::map<BaseContext,int> ignoreTheseContexts;
	size_t readUpToDepth, minDepth, maxDepth;
	bool applyMQFilter;
	int minMQ, maxMQ;
	int minFragmentLength, maxFragmentLength;
	bool _keepReadsLongerThanInsertSize;
	bool useStrand[2];
	bool useMate[2];

	//output filtered reads
	bool _updateBlacklist;
	TAlignmentBlacklist* blacklist;

	void _applyFilters();

public:
	TChromosomes chromosomes;
	TReadGroups readGroups;
	BamTools::BamAlignment bamAlignment;

	TBamFile();

	//filters
	void setAlignmentFiltersToKeepAll();
	void setAlignmentFilters(TParameters & params, TLog* logfile);
	void limitReadGroups(std::string readGroupList, TLog* logfile);

	//reading
	void open(std::string filename, uint32_t maxReadLength, bool indexNotRequired);
	bool readNextAlignment(); //TODO: make private
	bool readNextAlignmentThatPassesFilters(); //TODO: make private
	bool readNextAlignment(TAlignment & alignment);
	bool readNextAlignmentThatPassesFilters(TAlignment & alignment);

	bool jump(int chr, int position);
	void rewind();

	//getters
	std::string filename(){ return _filename; };
	uint32_t curPosition(){ return bamAlignment.Position; };
	uint16_t curReadGroupID(){ return _curReadGroupID; };
	bool chrChanged(){ return _chrChanged; };
	bool curIsLongerThanInsertSize();
	bool curPassedQC(){ return _QCFiltersPassed; };
	void fill(TAlignment & alignment); //TODO: make private

	uint64_t numAlignmentsRead(){ return _numAlignmentRead; };
	uint16_t maxReadLength(){ return _maxReadLength; };
};




#endif /* BAM_TBAMFILE_H_ */
