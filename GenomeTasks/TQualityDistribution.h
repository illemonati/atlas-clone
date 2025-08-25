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

#include "TAlignmentTraverser.h"
#include "SequencingError/TModels.h"

namespace GenomeTasks{

class TQualityTransformation final  {
private:
	TAlignmentTraverser _alnTraverser;
	std::vector<coretools::TCountDistributionVector<>> _transformations;
	bool _compareToOtherSeqErrors;
	std::string _label1, _label2;
	GenotypeLikelihoods::SequencingError::TModels _otherSeqErrors;

	void _traverseAlignments();

public:
	TQualityTransformation();
	void run();
};

}; // end namespace


#endif /* TQUALITYDISTRIBUTION_H_ */
