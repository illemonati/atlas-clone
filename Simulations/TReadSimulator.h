/*
 * TSimulatorQuality.h
 *
 *  Created on: Oct 5, 2017
 *      Author: vivian
 */

#ifndef TSIMULATORREAD_H_
#define TSIMULATORREAD_H_

#include <stdint.h>
#include <array>
#include <memory>
#include <string>
#include <vector>
#include <set>
#include "../BAM/TReadGroupInfo.h"
#include "GenotypeTypes.h"
#include "TAlignment.h"
#include "TCigar.h"
#include "TReadGroups.h"
#include "TSamFlags.h"
#include "PhredProbabilityTypes.h"
#include "TCategoricalDistribution.h"
#include "TReadGroupInfo.h"
#include "TLog.h"

namespace GenotypeLikelihoods { class TPMDType; }
namespace GenotypeLikelihoods { namespace SequencingError { class TModel; } }
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
	const BAM::TReadGroup &_readGroup;
	//const TReadGroupInfoEntry & _readGroupInfo;
	std::string _readNamePrefix;

	// required distributions
	TCategoricalDistribution<uint16_t> _fragmentLengthDistr;
	TCategoricalDistribution<PhredIntProbability> _qualityDist;
	TCategoricalDistribution<PhredIntProbability> _mappingQualityDist;

	// Additional info
	int _readXPos = 1;
	int _readYPos = 1;
	std::unique_ptr<TCategoricalDistribution<uint16_t>> _softClipDist5;
	std::unique_ptr<TCategoricalDistribution<uint16_t>> _softClipDist3;
	GenotypeLikelihoods::TPMDType const *_pmd = nullptr;
	std::array<GenotypeLikelihoods::SequencingError::TModel const *, 2> _recal;
	bool _hasRecal;

	// contamination
	double _contaminationRate = 0.;
	TSimulatorReference *_contaminationSource = nullptr;

	// alignment
	BAM::TSamFlags _flags;
	BAM::TAlignment _alignment;

	//initialization functions
	template <typename Distr>
	void _initDistribution(Distr & Dist, const TReadGroupInfoEntry & RGInfo, const BAM::RGInfo::InfoType & Info){
		coretools::instances::logfile().list(BAM::RGInfo::infos[Info].description, ": ", RGInfo[Info]);
		Dist.set(RGInfo[Info]);
	};

	// general functions
	double _calcMeanReadLength(const uint16_t maxLen) const;
	std::string _getNextReadName();
	void _simulateAlignmentDetails(const TGenomePosition & Position);
	bool _simulateContamination();
	void _addSoftclippedBases(std::vector<Base> & bases, const std::unique_ptr<TCategoricalDistribution<uint16_t>> & softClippedDist, BAM::TCigar & Cigar);
	void _simulateBasesQualities(BAM::TAlignment &alignment,
								 const std::vector<Base>& haplotype,
								 const uint16_t fragmentLength,
								 const uint16_t readLength,
								 bool readIsContaminated);

public:
	TReadSimulator(const BAM::TReadGroup & ReadGroup, const TReadGroupInfoEntry & RGInfo);
	virtual ~TReadSimulator() = default;

	//setters
	void setPMD(GenotypeLikelihoods::TPMDType const *Pmd);
	void setRecal(GenotypeLikelihoods::SequencingError::TModel const *Recal1, GenotypeLikelihoods::SequencingError::TModel const *Recal2);
	void setContamination(double rate, TSimulatorReference *source);

	//simulate
	virtual void simulate(const TGenomePosition & Position, const std::vector<Base>& Haplotype, TSimulatorBamFile &BamFile) = 0;
	virtual void writeUnwrittenAlignments(const genometools::TGenomePosition &, TSimulatorBamFile &){};

	//getters
	std::string name() const { return _readGroup.name_ID; };
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
	coretools::StrictlyPositive<uint16_t> _numCycles;

public:
	TReadSimulatorSingleEnd(const BAM::TReadGroup & ReadGroup, const TReadGroupInfoEntry & RGInfo);
	~TReadSimulatorSingleEnd() = default;

	void simulate(const TGenomePosition & Position, const std::vector<Base>& Haplotype, TSimulatorBamFile &BamFile) override;
	[[nodiscard]] double meanReadLength() const override;
};

//-------------------------------
// TSimulatorPairedEndReads
//-------------------------------
class TReadSimulatorPairedEnd final : public TReadSimulator {
private:
	std::array<coretools::StrictlyPositive<uint16_t>, 2> _numCycles;
	BAM::TSamFlags _mateFlags;

	std::multiset<BAM::TAlignment> _bamAlignmentSecondMates;

public:
	TReadSimulatorPairedEnd(const BAM::TReadGroup & ReadGroup, const TReadGroupInfoEntry & RGInfo);
	~TReadSimulatorPairedEnd() = default;

	void simulate(const TGenomePosition & Position, const std::vector<genometools::Base>& Haplotype, TSimulatorBamFile &BamFile) override;
	void writeUnwrittenAlignments(const genometools::TGenomePosition & Position, TSimulatorBamFile &BamFile) override;
	[[nodiscard]] double meanReadLength() const override;
};

}; // namespace Simulations

#endif /* TSIMULATORREAD_H_ */
