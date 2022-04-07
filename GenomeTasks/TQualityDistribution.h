/*
 * TQualityDistribution.h
 *
 *  Created on: Jun 2, 2020
 *      Author: phaentu
 */

#ifndef TQUALITYDISTRIBUTION_H_
#define TQUALITYDISTRIBUTION_H_

#include <string>
#include <vector>

#include "TGenome.h"
#include "TLog.h"
#include "TParameters.h"
#include "TRandomGenerator.h"
#include "TSequencingErrorModels.h"
#include "TTask.h"
#include "counters.h"

namespace GenomeTasks{

//-----------------------------------
// TQualityDistribution
//-----------------------------------
class TQualityDistribution:public TGenome_parsed{
private:
	coretools::TCountDistributionVector _qualDist;

	void _handleAlignment();

public:
	TQualityDistribution(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
	void compileQualityDistribution();
};

//-----------------------------------
// TQualityTransformation
//-----------------------------------
class TQualityTransformation:public TGenome_parsed{
private:
	std::vector<coretools::TCountDistributionVector> _transformations;
	bool _compareToOtherSeqErrors;
	std::string _label1, _label2;
	GenotypeLikelihoods::SequencingError::TModels _otherSeqErrors;

	void _handleAlignment();

public:
	TQualityTransformation(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
	void compileQualityTransformation();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_qualityDist:public coretools::TTask{
public:
	TTask_qualityDist(){ _explanation = "Printing Quality Distribution"; };

	void run(){
		using namespace coretools::instances;
		TQualityDistribution qualDist(parameters(), &logfile(), &randomGenerator());
		qualDist.compileQualityDistribution();
	};
};

class TTask_qualityTransformation:public coretools::TTask{
public:
	TTask_qualityTransformation(){ _explanation = "Printing Quality Transformation"; };

	void run(){
		using namespace coretools::instances;
		TQualityTransformation transformer(parameters(), &logfile(), &randomGenerator());
		transformer.compileQualityTransformation();
	};
};


}; // end namespace


#endif /* TQUALITYDISTRIBUTION_H_ */
