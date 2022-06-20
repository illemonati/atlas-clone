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
#include "TSequencingErrorModels.h"
#include "TTask.h"
#include "counters.h"

namespace GenomeTasks{

//-----------------------------------
// TQualityDistribution
//-----------------------------------
class TQualityDistribution:public TGenome_parsed{
private:
	coretools::TCountDistributionVector<> _qualDist;

	void _handleAlignment();

public:
	void compileQualityDistribution();
};

//-----------------------------------
// TQualityTransformation
//-----------------------------------
class TQualityTransformation:public TGenome_parsed{
private:
	std::vector<coretools::TCountDistributionVector<>> _transformations;
	bool _compareToOtherSeqErrors;
	std::string _label1, _label2;
	GenotypeLikelihoods::SequencingError::TModels _otherSeqErrors;

	void _handleAlignment();

public:
	TQualityTransformation();
	void compileQualityTransformation();
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

class TTask_qualityTransformation:public coretools::TTask{
public:
	TTask_qualityTransformation(){ _explanation = "Printing Quality Transformation"; };

	void run(){
		TQualityTransformation transformer;
		transformer.compileQualityTransformation();
	};
};


}; // end namespace


#endif /* TQUALITYDISTRIBUTION_H_ */
