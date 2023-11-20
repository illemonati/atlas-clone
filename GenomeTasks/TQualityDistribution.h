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

#include "TGenome_OLD.h"
#include "SequencingError/TModels.h"
#include "coretools/Main/TTask.h"
#include "coretools/Math/counters.h"

namespace GenomeTasks{

//-----------------------------------
// TQualityDistribution
//-----------------------------------
class TQualityDistribution:public old::TGenome_parsed{
private:
	coretools::TCountDistributionVector<> _qualDist;

	void _handleAlignment();

public:
	void compileQualityDistribution();
};

//-----------------------------------
// TQualityTransformation
//-----------------------------------
class TQualityTransformation:public old::TGenome_parsed{
private:
	std::vector<coretools::TCountDistributionVector<>> _transformations;
	bool _compareToOtherSeqErrors;
	std::string _label1, _label2;
	GenotypeLikelihoods::SequencingError::TModels _otherSeqErrors;

	void _handleAlignment();

public:
	TQualityTransformation();
	void run();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_qualityDist:public coretools::TTask{
public:
	TTask_qualityDist(){ _explanation = "Printing Quality Distribution"; };

	void run(){
		TQualityDistribution qualDist;
		qualDist.compileQualityDistribution();
	};
};

}; // end namespace


#endif /* TQUALITYDISTRIBUTION_H_ */
