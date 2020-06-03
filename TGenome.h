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
// TGenome_recalibrated
// A base class with BAM filters and a parsed, recalibrated alignment
//---------------------------------------------------------------
class TGenome_recalibrated:public TGenome_filtered{
protected:
	BAM::TAlignment _alignment;
	GenotypeLikelihoods::TGenotypeLikelihoodCalculator _genotypeLikelihoodCalculator;

	virtual void _traverseBAMPassedQC();
public:
	TGenome_recalibrated(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	virtual ~TGenome_recalibrated(){};

	void recalibrateBamFile(TParameters & params);
	void binQualityScores(TParameters & params);

};




#endif /* LOCI_H_ */
