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
#include "TChromosomes.h"
#include "TReadGroups.h"
#include "TAlignment.h"
#include "TAlignmentBlacklist.h"
#include "BAMData.h"

//-----------------------------------------------------
//TBamFileFilter
//-----------------------------------------------------
class TBamFileFilter_base{
protected:
	bool _keep;
	uint64_t _counter;
	std::string _errorString;
	bool _updateBlacklist;
	TAlignmentBlacklist* _blacklist;

	void _filterOut(const std::string & alignmentName, const bool & isReverseStrand){
		++_counter;
		if(_updateBlacklist){
			_blacklist->addToBlacklist(alignmentName, isReverseStrand, _errorString);
		}
	};


public:
	TBamFileFilter_base(){
		_counter = 0;
		_keep = true;
		_updateBlacklist = false;
		_blacklist = nullptr;
	};

	void keep(){
		_keep = true;
	};

	void setBlacklist(const TAlignmentBlacklist* Blacklist){
		_blacklist = Blacklist;
		_updateBlacklist = true;

	};

};
class TBamFileFilter:public TBamFileFilter_base{
public:
	TBamFileFilter(){};

	void filter(const std::string ErrorString){
		_keep = false;
		_errorString = ErrorString;
	};

	bool pass(const bool state, const std::string & alignmentName, const bool & isReverseStrand){
		if(state && !_keep){
			_filterOut(alignmentName, isReverseStrand);
			return false;
		}
		return true;
	};
};

class TBamFileFilterRange:public TBamFileFilter{
private:
	uint32_t _min, _max;
public:
	TBamFileFilterRange(){
		_min = 0;
		_max = UINT32_MAX;
	};

	void filter(const uint32_t Min, const uint32_t Max, const std::string ErrorString){
		_keep = false;
		_min = Min;
		_max = Max;
		_errorString = ErrorString;
	};

	bool pass(const uint32_t value, const std::string & alignmentName, const bool & isReverseStrand){
		if(!_keep && (value < _min || value > _max)){
			_filterOut(alignmentName, isReverseStrand);
		}
		return true;
	};

};
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
 	bool _keepAll;

 	TBamFileFilter _duplicateFilter;
 	TBamFileFilter _softClippedFilter;
 	TBamFileFilter _improperPairsFilter;
 	TBamFileFilter _unmappedFilter;
 	TBamFileFilter _failedQCFilter;
 	TBamFileFilter _secondaryFilter;
 	TBamFileFilter _supplementaryFilter;
 	TBamFileFilter _longerThanFragmentFilter;
 	TBamFileFilter _readGroupFilter;
 	TBamFileFilter _fwdStrandFilter;
 	TBamFileFilter _revStrandFilter;
 	TBamFileFilter _firstMateFilter;
 	TBamFileFilter _secondMateFilter;
 	TBamFileFilterRange _mappingQualityFilter;
 	TBamFileFilterRange _fragmentLengthfilter;

 	//base filters -> to parser?
 	std::map<BaseContext,int> ignoreTheseContexts;



	//output filtered reads
	bool _updateBlacklist;
	TAlignmentBlacklist* _blacklist;

	void _applyFilters();

public:
	TChromosomes chromosomes;
	TReadGroups readGroups;


	TBamFile();

	//filters
	void setAlignmentFiltersToKeepAll();
	void setAlignmentFilters(TParameters & params, TLog* logfile);
	void limitReadGroups(std::string readGroupList, TLog* logfile);

	//reading
	void open(const std::string filename, const uint32_t maxReadLength, const bool indexNotRequired);
	bool readNextAlignment(); //TODO: make private
	bool readNextAlignmentThatPassesFilters(); //TODO: make private
	bool readNextAlignment(TAlignment & alignment);
	bool readNextAlignmentThatPassesFilters(TAlignment & alignment);

	bool jump(int chr, int position);
	void rewind();

	//writing
	void openOutput(std::string filename);
	void writeCurAlignment();
	void writeAlignment(TAlignment & alignment, const TGenotypeMap & genoMap, const TQualityMap & qualityMap);

	//getters
	std::string filename(){ return _filename; };
	uint32_t curPosition(){ return _curBamAlignment.Position; };
	uint16_t curReadGroupID(){ return _curReadGroupID; };
	bool chrChanged(){ return _chrChanged; };
	bool curIsLongerThanInsertSize(){
		return _curBamAlignment.IsPaired() && abs(_curBamAlignment.InsertSize) < _curBamAlignment.AlignedBases.length();
	};
	bool curPassedQC(){ return _QCFiltersPassed; };
	void fill(TAlignment & alignment);

	uint64_t numAlignmentsRead(){ return _numAlignmentRead; };
	uint16_t maxReadLength(){ return _maxReadLength; };
};




#endif /* BAM_TBAMFILE_H_ */
