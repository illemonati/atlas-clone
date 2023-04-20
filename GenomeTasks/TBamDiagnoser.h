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

#include "TBamFilter.h"
#include "TGenome.h"
#include "coretools/Main/TTask.h"
#include "coretools/Math/counters.h"

namespace GenomeTasks{

//-------------------------------------------
// TBamDiagnoser
// A class to get basic characteristics of a BAM file
//-------------------------------------------

class TBamDiagnoser:public TGenome_filtered{
private:
	BAM::TQualityFilter _qualFilter;
	std::vector<std::string> _readGroupNames;
	bool _chromStats = false;

    // distributions
    coretools::TCountDistributionVector<> _passedQC;
    std::vector<coretools::TCountDistributionVector<>> _readLength;
    std::vector<coretools::TCountDistributionVector<>> _usableLength;
    std::vector<coretools::TCountDistributionVector<>> _softClippedLength;
   	std::vector<coretools::TCountDistributionVector<>> _mappingQuality;
    std::vector<coretools::TCountDistributionVector<>> _fragmentLength;

	void _writeHistogram(const std::vector<coretools::TCountDistributionVector<>> & distVec, const std::string& header, const std::string& name);
	void _writeTable(const coretools::TCountDistributionVector<> & distVec, const std::string& header, const std::string &name);
    void _handleAlignment() override;

public:
	TBamDiagnoser();
	double meanOverAllReadGroups(const std::vector<coretools::TCountDistributionVector<>> &vec);
	double meanForChromosome(const std::vector<coretools::TCountDistributionVector<>> &vec, size_t chromRefID);
	size_t maxOverAllReadGroups(const std::vector<coretools::TCountDistributionVector<>> &vec);
	size_t maxForChromosome(const std::vector<coretools::TCountDistributionVector<>> &vec, size_t chromRefID);
	size_t countsOverAllReadGroups(const std::vector<coretools::TCountDistributionVector<>> &vec);
	size_t countsForChromosome(const std::vector<coretools::TCountDistributionVector<>> &vec, size_t chromRefID);
	size_t countsLargerZeroOverAllReadGroups(const std::vector<coretools::TCountDistributionVector<>> &vec);
	size_t countsLargerZeroForChromosome(const std::vector<coretools::TCountDistributionVector<>> &vec, size_t chromRefID);
	size_t sumOverAllReadGroups(const std::vector<coretools::TCountDistributionVector<>> &vec);
	size_t sumForChromosome(const std::vector<coretools::TCountDistributionVector<>> &vec, size_t chromRefID);

	void run();
};

}; // end namespace

#endif /* TBAMDIAGNOSER_H_ */
