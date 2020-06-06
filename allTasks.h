/*
 * TTask.h
 *
 *  Created on: Mar 31, 2019
 *      Author: phaentu
 */

#ifndef ALLTASKS_H_
#define ALLTASKS_H_

#include "commonutilities/TTask.h"

//---------------------------------------------------------------------------
// all includes necessary for the application
//---------------------------------------------------------------------------
#include "TDistanceEstimator.h"
#include "Simulations/TSimulator.h"
#include "TMajorMinor.h"
#include "TAlleleCountEstimator.h"
#include "TAlleleFrequencyEstimator.h"
#include "TInbreedingEstimator.h"
#include "TPolymorhicWindowIdentifier.h"
#include "TBamDiagnoser.h"
#include "TBamFilter.h"
#include "VCF/TVcfConverter.h"
#include "VCF/TVcfDiagnostics.h"
#include "VCF/TVCFCompare.h"


//---------------------------------------------------------------------------
// Missing
//---------------------------------------------------------------------------

/*

class TTask_splitRGbyLength:public TTask_atlas{ -> splitmerger
public:
	TTask_splitRGbyLength(){ _explanation = "Splitting single end read groups in a BAM file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenomeWindows genome(logfile, parameters, _randomGenerator);
		genome.splitSingleEndReadGroups(parameters);
	};
};

class TTask_mateInfo:public TTask_atlas{ -> pileup
public:
	TTask_mateInfo(){ _explanation = "Writing read information per site"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenomeWindows genome(logfile, parameters, _randomGenerator);
		genome.printMateInformationPerSite(parameters);
	};
};

class TTask_contextStats:public TTask_atlas{
public:
	TTask_contextStats(){ _explanation = "Writing context statistics to file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenomeWindows genome(logfile, parameters, _randomGenerator);
		genome.contextStats(parameters);
	};
};

class TTask_writeDepthPerSite:public TTask_atlas{ -> pileup
public:
	TTask_writeDepthPerSite(){ _explanation = "Writing sequencing depth for each site"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenomeWindows genome(logfile, parameters, _randomGenerator);
		genome.writeDepthPerSite(parameters);
	};
};


class TTask_BQSR:public TTask_atlas{ -> suppress?
public:
	TTask_BQSR(){
		_explanation = "Estimating BQSR error re-calibration parameters";
		_citations.push_back("Hofmanova et al. (2016) PNAS");
	};

	void run(TParameters & parameters, TLog* logfile){
		TGenomeWindows genome(logfile, parameters, _randomGenerator);
		genome.BQSR(parameters);
	};
};

class TTask_binQualityScores:public TTask_atlas{ -> a feature of all BAM writing functions, e.g. BamFilter
public:
	TTask_binQualityScores(){ _explanation = "Binning quality scores"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenomeWindows genome(logfile, parameters, _randomGenerator);
		genome.binQualityScores(parameters);
	};
};


class TTask_recalBAM:public TTask_atlas{ -> feature of all BAM writing functions that parse, e.g. BamFilter
public:
	TTask_recalBAM(){ _explanation = "Recalibrating quality scores in a BAM file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenomeWindows genome(logfile, parameters, _randomGenerator);
		genome.recalibrateBamFile(parameters);
	};
};

class TTask_testBED:public TTask_atlas{ ->should be a proper integration test
public:
	TTask_testBED(){ _explanation = "Testing BED files"; };

	void run(TParameters & parameters, TLog* logfile){
		TBed bed;
		bed.test();
	};
};


*/

//---------------------------------------------------------------------------
// BAM file tools
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
// Depth
//---------------------------------------------------------------------------




//---------------------------------------------------------------------------
// PMD
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
// Recalibration
//---------------------------------------------------------------------------



//---------------------------------------------------------------------------
// Caller
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Diversity (theta)
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Population tools
//---------------------------------------------------------------------------




//---------------------------------------------------------------------------
// Simulations
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// VCF tools
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
// Debugging
//---------------------------------------------------------------------------



#endif /* ALLTASKS_H_ */
