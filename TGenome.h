/*
 * loci.h
 *
 *  Created on: Feb 19, 2015
 *      Author: wegmannd
 */

#ifndef LOCI_H_
#define LOCI_H_

#include <TGenotypeData.h>

#include "TWindow.h"
#include "gzstream.h"
#include "bamtools/api/BamWriter.h"
#include "TLog.h"
#include "TBed.h"
#include "TAlignmentParser.h"
#include "TQualityMap.h"
#include "TReadList.h"
#include <typeinfo>
#include <map>
#include <algorithm>
#include "counters.h"

#include "TRecalibration.h"
#include "TAllelicDepthCounts.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "TGenotypePrior.h"
#include "TBamFilter.h"

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

	TGenotypeMap _genoMap;
	TQualityMap _qualMap;


	virtual void _openBamForWriting(const std::string filename, BAM::TOutputBamFile & outBam);

public:
	TGenome_basic(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	virtual ~TGenome_basic(){};

	void mergeReadGroups(TParameters & params);
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

	//reference
	BAM::TFastaBuffer _reference;
	void _openReference(bool required = false);

	//read trimming
	bool _trimReads;
	int _trimmingLength3Prime;
	int _trimmingLength5Prime;

	//filters
	bool _applyQualityFilter;
	TQualityFilter _qualityFilter;
	bool _applyContextFilter;
	std::map<BaseContext,int> _ignoreTheseContexts;

	//functions for initialization
	void _setReadTrimming(TParameters & params);
	void _setQualityFilter(TParameters & params);
	void _setContextFilter(TParameters & params);
	void _setQualityRangeForPrinting(TParameters & params);

	void _parseAlignment(BAM::TAlignment & alignment);
	virtual void _traverseBAMPassedQC();
public:
	TGenome_parsed(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	virtual ~TGenome_parsed(){};

	void recalibrateBamFile(TParameters & params);
	void binQualityScores(TParameters & params);

};

//---------------------------------------------------------------
// TGenome_windows
// A base class to traverse a BAM file in windows
//---------------------------------------------------------------
class TGenome_windows:public TGenome_parsed{
protected:
	BAM::TAlignment* oldAlignment;
	bool oldAlignmentInitialized;
	bool oldAlignmentMustBeConsidered;

	//counters
	bool hasWindowIndent;

	//window params
	bool windowsPredefined;
	TBed* predefinedWindows;
	unsigned int windowSize;
	unsigned int numWindowsOnChr;
	unsigned int windowNumber;

	//window limits
	long limitWindows;
	int skipWindows;

	//window filters
	double maxMissing;
	double maxRefN;

	//masks
	bool doMasking, considerRegions;
	BAM::TBedReader* mask;

	//filters
	bool applyDepthFilter;
	uint32_t readUpToDepth, minDepth, maxDepth;

	//contructor functions
	void _setWindowParameters(TParameters & params);
	void _setParsingLimits(TParameters & params);
	void _setWindowFilters(TParameters & params);
	void _setSiteFilters(TParameters & params);
	void _setMasks(TParameters & params);

	//functions to traverse BAM in windows
	TWindow window;
	void _jumpToEnd();
	void _restartChromosomes(TWindow_base & window);
	void _moveChromosome(TWindow_base & window);
	bool _moveToNextWindowOnChr(TWindow_base & window);
	bool _moveToNextPredefinedWindow(TWindow_base & window);
	bool _moveWindow(TWindow_base & window);
	bool _readAlignment();
	bool _applyFilters();
	bool _fillAlignment(BAM::TAlignment & alignment);
	void _readAlignmentsIntoWindow(TWindow & window);
	void _applyWindowFilters(TWindow_base & window);
	void _readAlignmentsIntoWindow(TWindow & window);

	void traverseBAM();

public:
	TGenome_windows(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);

};

#endif /* LOCI_H_ */
