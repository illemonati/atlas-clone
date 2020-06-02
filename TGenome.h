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
#include "TAlignmentMerger.h"
#include "TAllelicDepthCounts.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "TGenotypePrior.h"

//---------------------------------------------------------------
// TGenome_basic
// A base class without filters and genotype likelihoods
//---------------------------------------------------------------
class TGenome_basic{
protected:
	TLog* _logfile;
	BAM::TBamFile _bamFile;
	TRandomGenerator* _randomGenerator;
	std::string _outputName;

	TGenotypeMap _genoMap;
	TQualityMap _qualMap;

public:
	TGenome_basic(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);

	void mergeReadGroups(TParameters & params);
};

//---------------------------------------------------------------
// TGenome_filtered
// A base class without recalibration but BAM filters enabled
//---------------------------------------------------------------
class TGenome_filtered:public TGenome_basic{
protected:
public:
	TGenome_filtered(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
};


//---------------------------------------------------------------
// TGenome_recalibrated
// A base class with BAM filters and a parsed, recalibrated alignment
//---------------------------------------------------------------
class TGenome_recalibrated:public TGenome_filtered{
protected:
	BAM::TAlignment _alignment;
	TGenotypeLikelihoodCalculator _genotypeLikelihoodCalculator;

	void _traverseBAMPassedQC();
	virtual void _handleAlignment(){ throw "_handleAlignment() not implemented for base class TGenome_recalibrated!"; };

public:
	TGenome_recalibrated(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	virtual ~TGenome_recalibrated(){};

	void printQualityDistribution(TParameters & params);
	void printQualityTransformation(TParameters & params);

	void recalibrateBamFile(TParameters & params);
	void binQualityScores(TParameters & params);

};




#endif /* LOCI_H_ */
