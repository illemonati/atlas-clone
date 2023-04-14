/*
 * TReadSimulator.h
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
#include "SequencingError/TModels.h"
#include "genometools/GenotypeTypes.h"
#include "TAlignment.h"
#include "TCigar.h"
#include "TReadGroups.h"
#include "TSamFlags.h"
#include "genometools/PhredProbabilityTypes.h"
#include "coretools/Math/TCategoricalDistribution.h"
#include "TReadGroupInfo.h"
#include "coretools/Main/TLog.h"
#include "../BAM/TBamFile.h"
#include "TSimulatedOutputFile.h"

namespace GenotypeLikelihoods { class TPMDType; }
namespace GenotypeLikelihoods { namespace SequencingError { class TReadGroupModels; } }
namespace GenotypeLikelihoods { namespace SequencingError { class TModel; } }
namespace Simulations { class TSimulatorBamFile; }
namespace Simulations { class TSimulatorReference; }
namespace Simulations { class TSimulatedOutputFile; }

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

	// contamination
	double _contaminationRate = 0.;
	TSimulatorReference *_contaminationSource = nullptr;

	GenotypeLikelihoods::SequencingError::TReadGroupModels _recalModels;

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
	void setContamination(double rate, TSimulatorReference *source);

	//simulate
	virtual void simulate(const TGenomePosition & Position, const std::vector<Base> & Haplotype, TSimulatedOutputFile & simulatedFile);
	virtual void simulate(const TGenomePosition & Position, const std::vector<genometools::Base> & Haplotype, BAM::TOutputBamFile & simulatedFile);

		/**
		 * 
		 * Used polymorphism to circumnavigate the problem with writeAlignmentLater in PairedEnd::simulate method implementation
		 * in this way, it is mandatory to use a BAM file in simulate in paired end but in the case of singleEnd it just needs TSimulatedOutputFile
		 * as it is the parent class of FASTQ and BAM.
		 * Override in child classes ensures the existance of said simulate methods
		 * 
		*/


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

	void SequenceSimulator(int sequenceLength);

	/***
	 * 
	 * Change BAM::TOutputBamFile to new TSimulatedOutputFile
	 * 
	*/

	void simulate(const TGenomePosition & Position, const std::vector<Base> & Haplotype, TSimulatedOutputFile & simulatedFile) override;

	[[nodiscard]] double meanReadLength() const override;

	void TReadSimulatorSingleEndFASTQ();
};

//-------------------------------
// TSimulatorPairedEndReads
//-------------------------------
class TReadSimulatorPairedEnd final : public TReadSimulator {
private:
	BAM::TAlignment _secondMate;
	BAM::TSamFlags _mateFlags;
	std::array<coretools::StrictlyPositive<uint16_t>, 2> _numCycles;

public:
	TReadSimulatorPairedEnd(const BAM::TReadGroup & ReadGroup, const TReadGroupInfoEntry & RGInfo);
	~TReadSimulatorPairedEnd() = default;

	void simulate(const TGenomePosition & Position, const std::vector<genometools::Base> & Haplotype, BAM::TOutputBamFile & simulatedFile) override;
	[[nodiscard]] double meanReadLength() const override;
};

}; // namespace Simulations

#endif /* TSIMULATORREAD_H_ */
