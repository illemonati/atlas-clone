/*
 * TReadSimulator.h
 *
 *  Created on: Oct 5, 2017
 *      Author: vivian
 */

#ifndef TSIMULATORREAD_H_
#define TSIMULATORREAD_H_

#include <array>
#include <memory>
#include <set>
#include <stdint.h>
#include <string>
#include <vector>

#include "../BAM/TBamFile.h"
#include "../BAM/TReadGroupInfo.h"
#include "PMD/TModel.h"
#include "SequencingError/TModels.h"
#include "TAlignment.h"
#include "TCigar.h"
#include "TReadGroupInfo.h"
#include "TReadGroups.h"
#include "TSamFlags.h"
#include "coretools/Main/TLog.h"
#include "coretools/Math/TCategoricalDistribution.h"
#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"

namespace GenotypeLikelihoods { namespace SequencingError { class TReadGroupModels; } }
namespace GenotypeLikelihoods { namespace SequencingError { struct TModel; } }
namespace Simulations { class TSimulatorBamFile; }
namespace Simulations { class TSimulatorReference; }

namespace Simulations {

using genometools::Base;
using genometools::PhredIntProbability;
using genometools::TGenomePosition;
using coretools::probdist::TCategoricalDistribution;
using BAM::RGInfo::TReadGroupInfoEntry;

//-------------------------------------
// TSimulatorRead
// a pure abstract base-class
//-------------------------------------
class TReadSimulator{
protected:
	const BAM::TReadGroup *_readGroup;
	//const TReadGroupInfoEntry & _readGroupInfo;
	std::string _readNamePrefix;

	// required distributions
	TCategoricalDistribution<size_t> _fragmentLengthDistr;
	TCategoricalDistribution<PhredIntProbability> _qualityDist;
	TCategoricalDistribution<PhredIntProbability> _mappingQualityDist;

	// Additional info
	int _readXPos = 1;
	int _readYPos = 1;
	std::unique_ptr<TCategoricalDistribution<size_t>> _softClipDist5;
	std::unique_ptr<TCategoricalDistribution<size_t>> _softClipDist3;
	const GenotypeLikelihoods::PMD::TModel *_pmd;
	const GenotypeLikelihoods::SequencingError::RGModels _recal;

	// contamination
	double _contaminationRate = 0.;
	TSimulatorReference *_contaminationSource = nullptr;


	// alignment
	BAM::TSamFlags _flags;
	BAM::TAlignment _alignment;

	//initialization functions
	template <typename Distr>
	void _initDistribution(Distr & Dist, const TReadGroupInfoEntry & RGInfo, const BAM::RGInfo::InfoType & Info){
		coretools::instances::logfile().list(coretools::str::capitalizeFirst(BAM::RGInfo::infos[Info].description), ": ", RGInfo.getString(Info));
		Dist.set(RGInfo.getString(Info));
	};

	// general functions
	double _calcMeanReadLength(size_t maxLen) const;
	std::string _getNextReadName();
	void _simulateAlignmentDetails(const TGenomePosition & Position);
	bool _simulateContamination();
	void _addSoftclippedBases(std::vector<Base> & bases, const std::unique_ptr<TCategoricalDistribution<size_t>> & softClippedDist, BAM::TCigar & Cigar);
	void _simulateBasesQualities(BAM::TAlignment &alignment, const std::vector<Base> &haplotype, size_t fragmentLength,
								 size_t readLength, bool readIsContaminated);

public:
	TReadSimulator(const BAM::TReadGroup & ReadGroup, const TReadGroupInfoEntry & RGInfo, const GenotypeLikelihoods::PMD::TModel & Pmd, const GenotypeLikelihoods::SequencingError::RGModels& Recal);
	virtual ~TReadSimulator() = default;

	//setters
	void setPMD(GenotypeLikelihoods::PMD::TModel const *Pmd);
	void setContamination(double rate, TSimulatorReference *source);

	//simulate
	virtual void simulate(const TGenomePosition & Position, const std::vector<Base> & Haplotype, BAM::TOutputBamFile &BamFile) = 0;

	//getters
	std::string name() const { return _readGroup->name_ID; };
	[[nodiscard]] virtual double meanReadLength() const = 0;
	double maxFragmentLength() {
		return _fragmentLengthDistr.max();
	};
};

//-------------------------------
// TSimulatorSingleEndRead
//-------------------------------
class TReadSimulatorSingleEnd final : public TReadSimulator {
private:
	coretools::StrictlyPositive<size_t> _numCycles;

public:
	TReadSimulatorSingleEnd(const BAM::TReadGroup & ReadGroup, const TReadGroupInfoEntry & RGInfo, const GenotypeLikelihoods::PMD::TModel & Pmd, const GenotypeLikelihoods::SequencingError::RGModels& Recal);
	~TReadSimulatorSingleEnd() = default;

	void simulate(const TGenomePosition & Position, const std::vector<Base> & Haplotype, BAM::TOutputBamFile &BamFile) override;
	[[nodiscard]] double meanReadLength() const override;
};

//-------------------------------
// TSimulatorPairedEndReads
//-------------------------------
class TReadSimulatorPairedEnd final : public TReadSimulator {
private:
	BAM::TAlignment _secondMate;
	BAM::TSamFlags _mateFlags;
	std::array<coretools::StrictlyPositive<size_t>, 2> _numCycles;

public:
	TReadSimulatorPairedEnd(const BAM::TReadGroup & ReadGroup, const TReadGroupInfoEntry & RGInfo, const GenotypeLikelihoods::PMD::TModel & Pmd, const GenotypeLikelihoods::SequencingError::RGModels& Recal);
	~TReadSimulatorPairedEnd() = default;

	void simulate(const TGenomePosition & Position, const std::vector<genometools::Base> & Haplotype, BAM::TOutputBamFile &BamFile) override;
	[[nodiscard]] double meanReadLength() const override;
};

}; // namespace Simulations

#endif /* TSIMULATORREAD_H_ */
