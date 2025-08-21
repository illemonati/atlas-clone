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
#include "TRead.h"
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
	genometools::TGenomePosition _previousAlignmentPosition;

	//current read
	TRead _read;

	//alignment filters
	std::vector<size_t> _numNotAligned;
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
	bool jumpToEnd();

	size_t refID() const noexcept { return curPosition().refID(); }
	const genometools::TChromosome &curChromosome() const noexcept { return _chromosomes[refID()]; }
	constexpr bool chrChanged() const noexcept {
		return refID() != _previousAlignmentPosition.refID();
	}

	//writing
	void writeCurAlignment(TOutputBamFile & out);

	//getters for cur alignment
	const TRead& curRead() const {return _read;}
	genometools::TGenomePosition curPosition() const { return _read.position; }

	//modify cur alignment
	void setCurReadGroup(size_t id);
	void addCurSamField(const std::string& tag, float value);

	//other getters
	const std::string& filename() const noexcept { return _filename; }
	size_t maxReadLength() const { return _filters.range(FilterType::ReadLength).max(); }
	const coretools::TCountDistributionVector<>& numAlignmentReadPerReadGroupPerChromosome() const noexcept { return _numAlignmentReadPerReadGroupPerChromosome; }
	size_t numReadGroups() const noexcept { return _readGroups.size(); }
	double averageDepth();


	//progress reporting
	void printSummary(std::string_view outputName) const;
	void startProgressReporting(bool indent = true) const;
	void printProgress() const;
	void printEndWithSummary(std::string_view outputName, bool indent = true) const;
};

} // namespace BAM

#endif /* BAM_TBAMFILE_H_ */
