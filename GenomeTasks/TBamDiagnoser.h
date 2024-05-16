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
#include "coretools/Files/TOutputFile.h"
#include "genometools/GenomePositions/TGenomePosition.h"

namespace GenomeTasks{

class TBamDiagnoser final : public TBamReadTraverser<ReadType::Filtered> {
private:
	struct TOld {
		std::string name = "";
		size_t length    = 0;
		bool isReversed  = false;

		static constexpr size_t BIG = -1;
		genometools::TGenomePosition position{BIG, BIG};
	};
	TQualityFilter _qualFilter;
	std::vector<std::string> _readGroupNames;
	bool _chromStats = false;
	bool _identifyDuplicates = false;
	std::vector<TOld> _old;
	genometools::TGenomePosition _oldPosition;
	coretools::TOutputFile _duplicateFile;
	std::vector<size_t> _startCounter;
	size_t _allStart = 1;

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
    std::vector<coretools::TCountDistributionVector<>> _readStart;
    coretools::TCountDistributionVector<> _allReadStart;

    void _handleAlignment() override;

public:
	TBamDiagnoser();
	void run();
};

}; // end namespace

#endif /* TBAMDIAGNOSER_H_ */
