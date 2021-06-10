/*
 * TQualityDistribution.h
 *
 *  Created on: Jun 2, 2020
 *      Author: phaentu
 */

#ifndef TQUALITYDISTRIBUTION_H_
#define TQUALITYDISTRIBUTION_H_

#include "TGenome.h"
#include "TTask.h"

namespace GenomeTasks{

//-----------------------------------
// TQualityDistribution
//-----------------------------------
class TQualityDistribution:public TGenome_parsed{
private:
	coretools::TCountDistributionVector _qualDist;

	void _handleAlignment();

public:
	TQualityDistribution(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
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
	GenotypeLikelihoods::TSequencingErrorModels _otherSeqErrors;

	void _handleAlignment();

public:
	TQualityTransformation(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void compileQualityTransformation();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_qualityDist:public coretools::TTask{
public:
	TTask_qualityDist(){ _explanation = "Printing Quality Distribution"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TQualityDistribution qualDist(Parameters, Logfile, _randomGenerator);
		qualDist.compileQualityDistribution();
	};
};

class TTask_qualityTransformation:public coretools::TTask{
public:
	TTask_qualityTransformation(){ _explanation = "Printing Quality Transformation"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TQualityTransformation transformer(Parameters, Logfile, _randomGenerator);
		transformer.compileQualityTransformation();
	};
};


}; // end namespace


#endif /* TQUALITYDISTRIBUTION_H_ */
