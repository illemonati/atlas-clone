/*
 * TGenome.h
 *
 *  Created on: Feb 19, 2015
 *      Author: wegmannd
 */

#ifndef TGENOME_H_
#define TGENOME_H_

#include <memory>
#include <set>
#include <stdint.h>
#include <string>
#include <vector>

#include "genometools/BED/TBed.h"
#include "TAlignment.h"
#include "TBamFile.h"
#include "genometools/GenomePositions/TChromosomes.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "TGenotypeData.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "coretools/Math/TNumericRange.h"
#include "coretools/TTimer.h"
#include "TWindow.h"
#include "genometools/TFastaReader.h"

namespace GenotypeLikelihoods {
class TSiteSubset;
}
namespace coretools {
class TLog;
}
namespace coretools {
class TParameters;
}
namespace coretools {
class TRandomGenerator;
}
namespace coretools {
class TSubsamplePicker;
}

namespace GenomeTasks {

//---------------------------------------------------------------
// TGenome_basic
// A base class without filters and genotype likelihoods
//---------------------------------------------------------------
class TGenome_basic {
protected:
	BAM::TBamFile _bamFile;
	std::string _outputName;

	void _openBamForWriting(const std::string &Filename, BAM::TOutputBamFile &OutBam);

public:
	TGenome_basic();
	virtual ~TGenome_basic() = default;
};

//---------------------------------------------------------------
// TGenome_filtered
// A base class without recalibration but BAM filters enabled
//---------------------------------------------------------------
class TGenome_filtered : public TGenome_basic {
protected:
	void _traverseBAMPassedQC();
	virtual void _handleAlignment() = 0;

public:
	TGenome_filtered();
};

//---------------------------------------------------------------
// TGenome_parsed
// A base class with BAM filters and a parsed, recalibrated alignment
//---------------------------------------------------------------
class TGenome_parsed : public TGenome_basic {
protected:
	BAM::TAlignment _alignment;
	GenotypeLikelihoods::TGenotypeLikelihoodCalculator _genotypeLikelihoodCalculator;
	GenotypeLikelihoods::TGenotypeLikelihoods _genoLik;

	// reference
	genometools::TFastaReader _reference;
	void _openReference(bool required = false);

	// read trimming
	bool _trimReads;
	int _trimmingLength3Prime;
	int _trimmingLength5Prime;

	// filters
	BAM::TQualityFilter _qualityFilter;
	BAM::TContextFilter _contextFilter;

	// functions for initialization
	void _setReadTrimming();

	void _parseAlignment(BAM::TAlignment &alignment);
	void _traverseBAMPassedQC();
	virtual void _handleAlignment() = 0;

public:
	TGenome_parsed();
};

//---------------------------------------------------------------
// TGenome_windows
// A base class to traverse a BAM file in windows
//---------------------------------------------------------------
class TGenome_windows : public TGenome_parsed {
protected:
	const genometools::TChromosomes &_chromosomes;
	std::vector<genometools::TChromosome>::const_iterator _curChromosome;

	// window params
	size_t _windowSize;
	size_t _numWindowsOnChr;
	size_t _windowNumber;
	bool _chrChangedWindow;

	// predefined windows
	genometools::TGenomeWindowList _predefinedWindows;
	std::multiset<genometools::TGenomeWindow>::iterator _curPredefinedWindow;

	// window limits
	size_t _limitWindows;
	size_t _skipWindows;

	// window filters
	double _maxMissing;
	double _maxRefN;

	// mask
	bool _doMasking, _considerRegions;
	genometools::TBed _mask;

	// sites
	std::unique_ptr<GenotypeLikelihoods::TSiteSubset> _subset;

	// site filters
	bool _applyDepthFilter;
	size_t _readUpToDepth;
	coretools::TNumericRange<size_t> _depthFilter;
	bool _filterCpG;
	size_t _downsampleDepth;
	std::unique_ptr<coretools::TSubsamplePicker> subsamplePicker;

	// tmp variables
	BAM::TAlignment _curAlignment;
	bool _hasWindowIndent;
	coretools::TTimer _windowTimer;

	// contructor functions
	void _setWindowParameters();
	void _setParsingLimits();
	void _setWindowFilters();
	void _setSiteFilters();
	void _setMasks();
	void _openSiteSubset(const std::string &filename);

	// functions to traverse BAM in windows
	GenotypeLikelihoods::TWindow _window;
	void _jumpToEnd();
	void _setCountersBeginningOfChromosome();
	bool _incrementWindow(GenotypeLikelihoods::TWindow_base &window);
	bool _moveToNextWindow(GenotypeLikelihoods::TWindow_base &window);
	bool _incrementPredefinedWindow();
	bool _moveToNextPredefinedWindow(GenotypeLikelihoods::TWindow_base &window);

	bool _moveWindow(GenotypeLikelihoods::TWindow_base &window);
	void _readAlignmentsIntoWindow(GenotypeLikelihoods::TWindow &window);
	void _applyWindowFilters(GenotypeLikelihoods::TWindow_base &window);
	bool _readAndParseAlignment();
	bool _readDataInNextWindow(GenotypeLikelihoods::TWindow &window);

	void _traverseBAMWindows();
	virtual void _handleWindow() = 0;

public:
	TGenome_windows();
};

}; // namespace GenomeTasks

#endif /* GENOME_H_ */
