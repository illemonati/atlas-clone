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

#include "TBamTraverser.h"
#include "SequencingError/TModels.h"

namespace GenomeTasks{

//-----------------------------------
// TQualityDistribution
//-----------------------------------
class TQualityDistribution final : public TBamTraverser<true> {
private:
	coretools::TCountDistributionVector<> _qualDist;

	void _handleAlignment(BAM::TAlignment& alignment) override;

public:
	void compileQualityDistribution();
};

//-----------------------------------
// TQualityTransformation
//-----------------------------------
class TQualityTransformation final : public TBamTraverser<true> {
private:
	std::vector<coretools::TCountDistributionVector<>> _transformations;
	bool _compareToOtherSeqErrors;
	std::string _label1, _label2;
	GenotypeLikelihoods::SequencingError::TModels _otherSeqErrors;

	void _handleAlignment(BAM::TAlignment& alignment) override;

public:
	TQualityTransformation();
	void run();
};

}; // end namespace


#endif /* TQUALITYDISTRIBUTION_H_ */
