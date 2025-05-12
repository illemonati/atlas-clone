/*
 * TBamFile.h
 *
 *  Created on: May 18, 2020
 *      Author: wegmannd
 */

#ifndef BAM_TBAMFILE_H_
#define BAM_TBAMFILE_H_

#include <string>
#include <vector>

#include "TBamFilters.h"
#include "api/BamAlignment.h"
#include "api/BamReader.h"
#include "api/SamHeader.h"

#include "coretools/TTimer.h"
#include "coretools/Types/probability.h"
#include "genometools/GenomePositions/TChromosomes.h"
#include "genometools/GenomePositions/TGenomePosition.h"

#include "TBamFilter.h"
#include "TCigar.h"
#include "TReadGroups.h"
#include "TSamHeader.h"

namespace BamTools {class BamWriter;}

namespace BAM{
class TAlignment;

//-----------------------------------------------------
//TBamFile
//-----------------------------------------------------
class TOutputBamFile; //forward declaration

class TBamFile{
private:
	constexpr static size_t _step = 100;
	constexpr static size_t _nope = -1;

	//BAM file
	std::string _filename;
	size_t _ID;
	BamTools::BamReader _bamReader;
	BamTools::SamHeader _bamHeader;
	size_t _fileSize = 0;

	//header
 	genometools::TChromosomes _chromosomes;
 	std::vector<genometools::TChromosome>::iterator _curChromosome;
 	TReadGroups _readGroups;
 	TSamHeader _samHeader;

 	//counters
 	coretools::TCountDistributionVector<> _numAlignmentReadPerReadGroupPerChromosome;
	size_t _numAlignmentRead        = 0;
	size_t _numAlignmentsPassedQC   = 0;
	size_t _maxNumReadsToRead       = -1;
	size_t _numIdentifiedDuplicates = 0;

	// downsample
	coretools::Probability _downProb;
	size_t _numDownsampled = 0;

	// duplicates
	bool _resetDuplicates = false;
	size_t _maxDupOverlap = _nope;
	struct TOld {
		size_t length{0};
		genometools::TGenomePosition position{_nope, _nope};
	};
	std::vector<TOld> _old;

	//current alignment
 	BamTools::BamAlignment _curBamAlignment;
 	TCigar _curCigar;
 	genometools::TGenomePosition _curAlignmentPosition, _previousAlignmentPosition;
	size_t _curReadGroupID      = 0;
	bool _chrChanged            = false;

	//alignment filters
	std::vector<size_t> _numNotAligned;
	bool _QCFiltersPassed  = false;
	size_t _numNoReadGroup = 0;

	TBamFilters _filters;

	//report progress
	mutable coretools::TTimer _timer;
	mutable size_t _progressFrequency   = 1000000;
	mutable size_t _lastProgressPrinted = 0;

	void _fillSamHeader();
	void _fillChromosomes();
	void _fillReadGroups();
	bool _readNextAlignmentFromFile();
	bool _applyFilters();
	bool _identifyDuplicate();
	void _writeFilteringStats(std::string_view outputName) const;

public:
	TBamFile(std::string_view Filename, size_t ID, bool EnableFilters = true);

	//access header info READ ONLY
	const genometools::TChromosomes& chromosomes() const noexcept { return _chromosomes; };
	const TReadGroups& readGroups() const noexcept { return _readGroups; };
	const TSamHeader& samHeader() const noexcept { return _samHeader; };

	//modify header info: know what you do!
	TReadGroups& readGroupsMutable() noexcept { return _readGroups; };

	//filters
	void setFilters(const TBamFilters& Filters);
	void curFilterOut();
	void filterOut(const TAlignment & Alignment);
	void setExternalFilterReason(std::string_view reason);

	//get filter status
	const TBamFilters& filters() const noexcept {return _filters;}
	const TBamFilter& filter(FilterType t) const noexcept {return _filters[t];}

	//reading
	bool readNextAlignment();
	bool readNextAlignmentThatPassesFilters();

	void fill(TAlignment & alignment) const;

	bool jump(const genometools::TGenomePosition Position);

	//writing
	void writeCurAlignment(TOutputBamFile & out);

	//getters for cur alignment
	const std::string &curName() const { return _curBamAlignment.Name; };
	genometools::TGenomePosition curPosition() const { return _curAlignmentPosition; };
	size_t refID() const noexcept { return _curChromosome->refID(); };
	const genometools::TChromosome &curChromosome() const noexcept { return *_curChromosome; };
	const TCigar &curCIGAR() const noexcept { return _curCigar; };
	constexpr size_t curReadGroupID() const noexcept { return _curReadGroupID; };
	constexpr bool chrChanged() const noexcept { return _chrChanged; };
	constexpr bool curPassedQC() const noexcept { return _QCFiltersPassed; };
	size_t curFragmentLength() const;
	uint16_t curMappingQuality() const noexcept { return _curBamAlignment.MapQuality; };
	bool curIsPaired() const { return _curBamAlignment.IsPaired(); };
	bool curIsProperPair() const { return _curBamAlignment.IsProperPair(); };
	bool curIsReverseStrand() const { return _curBamAlignment.IsReverseStrand(); };
	bool curIsDuplicate() const { return _curBamAlignment.IsDuplicate(); };
	bool curIsMapped() const { return _curBamAlignment.IsMapped(); };
	bool curIsFailedQC() const { return _curBamAlignment.IsFailedQC(); };
	bool curIsSecondary() const { return !_curBamAlignment.IsPrimaryAlignment(); };
	bool curIsSupplementary() const { return _curBamAlignment.IsSupplementary(); };
	bool curIsLongerThanFragment() const {
		return _curBamAlignment.IsProperPair() &&
			   _curBamAlignment.InsertSize < static_cast<int>(_curCigar.lengthAligned());
	};
	bool curIsFirstMate() const { return _curBamAlignment.IsFirstMate(); };
	bool curIsSecondMate() const { return _curBamAlignment.IsSecondMate(); };
	std::string curQuerySequence(const size_t start, const size_t length) const;

	//modify cur alignment
	void curSetNewReadGroup(size_t id);
	void curAddSamField(const std::string& tag, float value);

	//other getters
	const std::string& filename() const noexcept { return _filename; };
	size_t maxReadLength() const { return _filters.range(FilterType::ReadLength).max(); };
	const coretools::TCountDistributionVector<>& numAlignmentReadPerReadGroupPerChromosome() const noexcept { return _numAlignmentReadPerReadGroupPerChromosome; };
	size_t numReadGroups() const noexcept { return _readGroups.size(); };

	//progress reporting
	void printSummary(std::string_view outputName) const;
	void startProgressReporting(bool indent = true) const;
	void printProgress() const;
	void printEndWithSummary(std::string_view outputName, bool indent = true) const;
};

}; //end namespace


#endif /* BAM_TBAMFILE_H_ */
