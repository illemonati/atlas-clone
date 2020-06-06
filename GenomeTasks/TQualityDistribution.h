/*
 * TQualityDistribution.h
 *
 *  Created on: Jun 2, 2020
 *      Author: phaentu
 */

#ifndef TQUALITYDISTRIBUTION_H_
#define TQUALITYDISTRIBUTION_H_

#include "TGenome.h"

namespace GenomeTasks{

//-----------------------------------
// TQualityDistribution
//-----------------------------------
class TQualityDistribution:public TGenome_parsed{
private:
	TCountDistributionVector _qualDist;

	void _handleAlignment();

public:
	TQualityDistribution(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void compileQualityDistribution();
};

//-----------------------------------
// TQualityTransformation
//-----------------------------------
class TQualityTransformation:public TGenome_parsed{
private:
	std::vector<TCountDistributionVector> _transformations;
	bool _compareToOtherSeqErrors;
	std::string _label1, _label2;
	GenotypeLikelihoods::TSequencingErrorModels _otherSeqErrors;

	void _handleAlignment();

public:
	TQualityTransformation(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void compileQualityTransformation();
};

}; // end namespace


#endif /* TQUALITYDISTRIBUTION_H_ */
