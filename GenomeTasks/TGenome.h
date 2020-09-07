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
#include "TBed.h"
#include "TQualityMap.h"
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
	TLog* _logfile;
	TParameters* _params;
	BAM::TBamFile _bamFile;
	TRandomGenerator* _randomGenerator;
	std::string _outputName;

	GenotypeLikelihoods::TGenotypeMap _genoMap;
	BAM::TQualityMap _qualMap;

	virtual void _openBamForWriting(const std::string filename, BAM::TOutputBamFile & outBam);
	void _initialize(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);

public:
    TGenome_basic();
    TGenome_basic(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
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
    TGenome_filtered(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
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
	bool _applyQualityFilter;
	BAM::TQualityFilter _qualityFilter;
	bool _applyContextFilter;
	std::map<GenotypeLikelihoods::BaseContext,int> _ignoreTheseContexts;

	//functions for initialization
	void _setReadTrimming(TParameters & params);
	void _setQualityFilter(TParameters & params);
	void _setContextFilter(TParameters & params);

	void _parseAlignment(BAM::TAlignment & alignment);
	void _traverseBAMPassedQC();
public:
	TGenome_parsed(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	virtual ~TGenome_parsed(){};
};

//---------------------------------------------------------------
// TGenome_windows
// A base class to traverse a BAM file in windows
//---------------------------------------------------------------
class TGenome_windows:public TGenome_parsed{
protected:
	const BAM::TChromosomes& _chromosomes;
	std::vector<BAM::TChromosome>::const_iterator _curChromosome;

	//window params
	uint32_t _windowSize;
	uint32_t _numWindowsOnChr;
	uint32_t _windowNumber;
	bool _chrChangedWindow;

	//predefined windows
	BAM::TGenomeWindowList _predefinedWindows;
	std::multiset<BAM::TGenomeWindow>::iterator _curPredefinedWindow;

	//window limits
	long _limitWindows;
	uint32_t _skipWindows;

	//window filters
	double _maxMissing;
	double _maxRefN;

	//mask
	bool _doMasking, _considerRegions;
	BAM::TBed _mask;

	//sites
	std::unique_ptr<GenotypeLikelihoods::TSiteSubset> _subset;

	//site filters
	bool _applyDepthFilter;
	uint32_t _readUpToDepth, _minDepth, _maxDepth;
	bool _filterCpG;
	uint32_t _downsampleDepth;
	std::unique_ptr<TSubsamplePicker> subsamplePicker;

	//tmp variables
	BAM::TAlignment* _curAlignment;
	bool _hasWindowIndent;
	TTimer _windowTimer;

	//contructor functions
	void _setWindowParameters(TParameters & params);
	void _setParsingLimits(TParameters & params);
	void _setWindowFilters(TParameters & params);
	void _setSiteFilters(TParameters & params);
	void _setMasks(TParameters & params);
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
	TGenome_windows(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	~TGenome_windows();
};

}; //end namespace

#endif /* GENOME_H_ */
