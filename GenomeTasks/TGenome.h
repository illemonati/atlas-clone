/*
 * TGenome.h
 */

#ifndef TGENOME_H_
#define TGENOME_H_

#include <string>

#include "coretools/Main/TParameters.h"
#include "genometools/TFastaReader.h"

#include "TAlignment.h"
#include "TBamFile.h"
#include "TBaseFilter.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "TReadGroupInfo.h"

namespace GenomeTasks {

//---------------------------------------------------------------
// TGenome_basic
// A base class without filters and genotype likelihoods
//---------------------------------------------------------------
class TGenome {
private:
	BAM::TBamFile _bamFile;
	std::string _outputName;
	BAM::RGInfo::TReadGroupInfo _rgInfo;

public:
	TGenome();
	~TGenome();

	const BAM::TBamFile& bamFile() const noexcept {return _bamFile;}
	BAM::TBamFile& bamFile() noexcept {return _bamFile;}

	const BAM::RGInfo::TReadGroupInfo& rgInfo() const noexcept {return _rgInfo;}
	BAM::RGInfo::TReadGroupInfo& rgInfo() noexcept {return _rgInfo;}

	const std::string& outputName() const noexcept {return _outputName;}
};

template<typename AlignmentHandler>
class TFiltered {
	AlignmentHandler _handleAlignment;
public:
	TFiltered(TGenome& genome) {
		genome.bamFile().setFilters(BAM::TBamFilters{true});
	}
	void operator()() {
		_handleAlignment();
	}
};

template<typename AlignmentHandler, bool RefRequired>
class TParsed {
	TGenome *_genome;
	AlignmentHandler _handleAlignment;

	bool _trimReads;
	int _trim3;
	int _trim5;

	// filters
	TQualityFilter _qualityFilter;
	TContextFilter _contextFilter;

	GenotypeLikelihoods::TGenotypeLikelihoodCalculator _genotypeLikelihoodCalculator;
	genometools::TFastaReader _reference;

	void _parse(BAM::TAlignment& alignment) {
		alignment.parse(_genotypeLikelihoodCalculator.sequencingErrorModels());

		if (_trimReads) { alignment.trimRead(_trim3, _trim5); }
		alignment.filter(_qualityFilter); // always on
		if (_contextFilter) alignment.filter(_contextFilter);
		if (_reference.isOpen()) { alignment.addReference(_reference); }
	}

public:
	TParsed(TGenome &genome) : _genome(&genome), _genotypeLikelihoodCalculator(genome.rgInfo()) {
		using coretools::instances::logfile;
		using coretools::instances::parameters;
		genome.bamFile().setFilters(BAM::TBamFilters{true});

		if (parameters().exists("trim3") || parameters().exists("trim5")) {
			_trim3 = parameters().template get<int>("trim3", 0);
			if (_trim3 < 0) UERROR("trimming distance trim3 must be >= 0!");
			_trim5 = parameters().template get<int>("trim5", 0);
			if (_trim5 < 0) UERROR("trimming distance trim5 must be >= 0!");
			if (_trim3 > 0 || _trim5 > 0) {
				coretools::instances::logfile().list(
					"Will trim first ", _trim3, " and ", _trim5,
					" bases from the 3' and 5' end, respectively. (parameters 'trim3', 'trim5')");
			}
			_trimReads = true;
		} else {
			_trim3     = 0;
			_trim5     = 0;
			_trimReads = false;
		}

		if (parameters().exists("fasta")) {
			std::string fastaFile = parameters().template get<std::string>("fasta");
			logfile().list("Reading reference sequence from '" + fastaFile + "'. (parameter fasta)");
			_reference.open(fastaFile);
		} else if (RefRequired) {
			UERROR("No reference provided! (Use parameter fasta to provide a reference)");
		}
	}

	void operator()() {
		BAM::TAlignment alignment;
		_genome->bamFile().fill(alignment);
		_parse(alignment);
		_handleAlignment(alignment);
	}
};

template<typename AlignmentHandler> class TReadTraverser {
	TGenome _genome;

	// Functors
	AlignmentHandler _handleAlignment;

	void _traverseBAMPassedQC() {
		// parse through bam file
		_genome.bamFile().startProgressReporting();
		while (_genome.bamFile().readNextAlignmentThatPassesFilters()) {
			_handleAlignment();
			_genome.bamFile().printProgress();
		}
		// report
		_genome.bamFile().printEndWithSummary(_genome.outputName());
	}

public:
	TReadTraverser() :_handleAlignment(_genome) {}
};

}; // namespace GenomeTasks

#endif /* GENOME_H_ */
