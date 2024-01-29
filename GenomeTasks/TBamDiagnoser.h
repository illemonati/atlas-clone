/*
 * TBamDiagnoser.h
 *
 *  Created on: May 30, 2020
 *      Author: phaentu
 */

#ifndef TBAMDIAGNOSER_H_
#define TBAMDIAGNOSER_H_

#include <string>
#include <vector>

#include "TBamTraverser.h"
#include "coretools/Math/counters.h"
#include "genometools/GenomePositions/TGenomePosition.h"

namespace GenomeTasks{

class TBamDiagnoser final : public TBamReadTraverser<ReadType::Filtered> {
private:
	TQualityFilter _qualFilter;
	std::vector<std::string> _readGroupNames;
	bool _chromStats = false;
	std::vector<genometools::TGenomePosition> _oldPositions;
	genometools::TGenomePosition _oldPosition;

    // distributions
    coretools::TCountDistributionVector<> _passedQC;
	// std::vector per readgroup, countdistributionvector per chromosome
    std::vector<coretools::TCountDistributionVector<>> _readLength;
    std::vector<coretools::TCountDistributionVector<>> _readDist;
    coretools::TCountDistributionVector<> _allReadDist;
    std::vector<coretools::TCountDistributionVector<>> _usableLength;
    std::vector<coretools::TCountDistributionVector<>> _softClippedLength;
   	std::vector<coretools::TCountDistributionVector<>> _mappingQuality;
    std::vector<coretools::TCountDistributionVector<>> _fragmentLength;

    void _handleAlignment() override;

public:
	TBamDiagnoser() {
		_genome.bamFile().readGroups().fillVectorWithNames(_readGroupNames);
	}

	void run();
};

}; // end namespace

#endif /* TBAMDIAGNOSER_H_ */
