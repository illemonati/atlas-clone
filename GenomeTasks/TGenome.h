/*
 * loci.h
 *
 *  Created on: Feb 19, 2015
 *      Author: wegmannd
 */

#ifndef GENOME_H_
#define GENOME_H_


#include <typeinfo>
#include <map>
#include <algorithm>

#include "TGenotypeData.h"
#include "TWindow.h"
#include "gzstream.h"
#include "TLog.h"
#include "BED/TBed.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "TGenotypePrior.h"
#include "TBamFile.h"

namespace GenomeTasks{

//---------------------------------------------------------------
// TGenome_basic
// A base class without filters and genotype likelihoods
//---------------------------------------------------------------
class TGenome_basic{
protected:
	coretools::TLog* _logfile;
	coretools::TParameters* _params;
	BAM::TBamFile _bamFile;
	coretools::TRandomGenerator* _randomGenerator;
	std::string _outputName;

	virtual void _openBamForWriting(const std::string filename, BAM::TOutputBamFile & outBam);
	void _initialize(coretools::TParameters & Params, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);

public:
    TGenome_basic();
    TGenome_basic(coretools::TParameters & Params, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
    virtual ~TGenome_basic(){};
};

//---------------------------------------------------------------
// TGenome_filtered
// A base class without recalibration but BAM filters enabled
//---------------------------------------------------------------
class TGenome_filtered:public TGenome_basic{
protected:
	virtual void _traverseBAMPassedQC();
	virtual void _handleAlignment(){ throw "_handleAlignment() not implemented for base class TGenome_filtered!"; };

public:
    TGenome_filtered();
    TGenome_filtered(coretools::TParameters & Params, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
};


//---------------------------------------------------------------
// TGenome_parsed
// A base class with BAM filters and a parsed, recalibrated alignment
//---------------------------------------------------------------
class TGenome_parsed:public TGenome_filtered{
protected:
	BAM::TAlignment _alignment;
	GenotypeLikelihoods::TGenotypeLikelihoodCalculator _genotypeLikelihoodCalculator;
	GenotypeLikelihoods::TGenotypeLikelihoods _genoLik;

	//reference
	BAM::TFastaBuffer _reference;
	void _openReference(bool required = false);

	//read trimming
	bool _trimReads;
	int _trimmingLength3Prime;
	int _trimmingLength5Prime;

	//filters
	BAM::TQualityFilter _qualityFilter;
	BAM::TContextFilter _contextFilter;

	//functions for initialization
	void _setReadTrimming(coretools::TParameters & params);

	void _parseAlignment(BAM::TAlignment & alignment);
	void _traverseBAMPassedQC();
public:
	TGenome_parsed(coretools::TParameters & Params, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
	virtual ~TGenome_parsed(){};
};

//---------------------------------------------------------------
// TGenome_windows
// A base class to traverse a BAM file in windows
//---------------------------------------------------------------
class TGenome_windows:public TGenome_parsed{
protected:
	const genometools::TChromosomes& _chromosomes;
	std::vector<genometools::TChromosome>::const_iterator _curChromosome;

	//window params
	uint32_t _windowSize;
	uint32_t _numWindowsOnChr;
	uint32_t _windowNumber;
	bool _chrChangedWindow;

	//predefined windows
    genometools::TGenomeWindowList _predefinedWindows;
	std::multiset<genometools::TGenomeWindow>::iterator _curPredefinedWindow;

	//window limits
	long _limitWindows;
	uint32_t _skipWindows;

	//window filters
	double _maxMissing;
	double _maxRefN;

	//mask
	bool _doMasking, _considerRegions;
    genometools::TBed _mask;

	//sites
	std::unique_ptr<GenotypeLikelihoods::TSiteSubset> _subset;

	//site filters
	bool _applyDepthFilter;
	uint32_t _readUpToDepth;
	coretools::TNumericRange<uint32_t> _depthFilter;
	bool _filterCpG;
	uint32_t _downsampleDepth;
	std::unique_ptr<coretools::TSubsamplePicker> subsamplePicker;

	//tmp variables
	BAM::TAlignment* _curAlignment;
	bool _hasWindowIndent;
	coretools::TTimer _windowTimer;

	//contructor functions
	void _setWindowParameters(coretools::TParameters & params);
	void _setParsingLimits(coretools::TParameters & params);
	void _setWindowFilters(coretools::TParameters & params);
	void _setSiteFilters(coretools::TParameters & params);
	void _setMasks(coretools::TParameters & params);
	void _openSiteSubset(const std::string filename);

	//functions to traverse BAM in windows
	GenotypeLikelihoods::TWindow _window;
	void _jumpToEnd();
	void _setCountersBeginningOfChromosome();
	bool _incrementWindow(GenotypeLikelihoods::TWindow_base & window);
	bool _moveToNextWindow(GenotypeLikelihoods::TWindow_base & window);
	bool _incrementPredefinedWindow();
	bool _moveToNextPredefinedWindow(GenotypeLikelihoods::TWindow_base & window);

	bool _moveWindow(GenotypeLikelihoods::TWindow_base & window);
	void _readAlignmentsIntoWindow(GenotypeLikelihoods::TWindow & window);
	void _applyWindowFilters(GenotypeLikelihoods::TWindow_base & window);
    bool _readAndParseAlignment(BAM::TAlignment & alignment);
    bool _readDataInNextWindow(GenotypeLikelihoods::TWindow & window);

	void _traverseBAMWindows();
	virtual void _handleWindow(){ throw "_handleWindow() not implemented for base class TGenome_windows!"; };

public:
	TGenome_windows(coretools::TParameters & Params, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
	~TGenome_windows();
};

}; //end namespace

#endif /* GENOME_H_ */
