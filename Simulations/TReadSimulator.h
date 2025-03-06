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
#include <string>
#include <vector>

#include "coretools/Main/TLog.h"
#include "coretools/Math/TCategoricalDistribution.h"
#include "coretools/Types/probability.h"

#include "SequencingError/TModels.h"
#include "TAlignment.h"
#include "TReadGroupInfo.h"
#include "genometools/Genotypes/TwoBases.h"

namespace GenotypeLikelihoods::PMD {struct TModel;}

namespace BAM {
class TCigar;
class TOutputBamFile;
}

namespace Simulations { class TSimulatorReference; }

namespace Simulations {


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
	coretools::probdist::TCategoricalDistribution<size_t> _fragmentLengthDistr;
	coretools::probdist::TCategoricalDistribution<coretools::PhredInt> _qualityDist;
	coretools::probdist::TCategoricalDistribution<coretools::PhredInt> _mappingQualityDist;

	// Additional info
	int _readXPos = 1;
	int _readYPos = 1;
	std::unique_ptr<coretools::probdist::TCategoricalDistribution<size_t>> _softClipDist5;
	std::unique_ptr<coretools::probdist::TCategoricalDistribution<size_t>> _softClipDist3;
	const GenotypeLikelihoods::PMD::TModel *_pmd;
	const GenotypeLikelihoods::SequencingError::RGModels _recal;

	coretools::Probability _duplicationRate{0.};
	coretools::Probability _duplicationRateAmongSimulated{0.};
	coretools::Probability _baseN;
	coretools::Probability _refBias;

	// contamination
	double _contaminationRate = 0.;
	TSimulatorReference *_contaminationSource = nullptr;

	// alignment
	BAM::TAlignment _alignment;

	//initialization functions
	template <typename Distr>
	void _initDistribution(Distr & Dist, const BAM::RGInfo::TReadGroupInfoEntry & RGInfo, const BAM::RGInfo::InfoType & Info){
		coretools::instances::logfile().list(coretools::str::capitalizeFirst(BAM::RGInfo::infos[Info].description), ": ", RGInfo.getString(Info));
		Dist.set(RGInfo.getString(Info));
	};

	// general functions
	double _calcMeanReadLength(size_t MaxLen) const;
	std::string _getNextReadName();
	void _simulateAlignmentDetails(const genometools::TGenomePosition & Position);
	bool _simulateContamination();
	void _addSoftclippedBases(std::vector<genometools::Base> & Bases, const size_t &SoftClipLength, BAM::TCigar & Cigar);
	void _simulateBasesQualities(BAM::TAlignment& Alignment, const std::vector<genometools::TwoBase>& Haplotype,
								 bool FirstHaplo, size_t ReadLength, bool ReadIsContaminated);

	virtual void _simulate(const genometools::TGenomePosition &Position,
						   const std::vector<genometools::TwoBase> &Haplotype,
						   const TSimulatorReference &Reference)         = 0;
	virtual void _writeSimulatedAlignments(BAM::TOutputBamFile &BamFile) = 0;

public:
	TReadSimulator(const BAM::TReadGroup &ReadGroup, const BAM::RGInfo::TReadGroupInfoEntry &RGInfo,
				   const GenotypeLikelihoods::PMD::TModel &Pmd,
				   const GenotypeLikelihoods::SequencingError::RGModels &Recal);
	virtual ~TReadSimulator() = default;

	//setters
	void setPMD(GenotypeLikelihoods::PMD::TModel const *Pmd);
	void setContamination(double rate, TSimulatorReference *Source);

	//simulate
	size_t simulate(const genometools::TGenomePosition &Position, const std::vector<genometools::TwoBase> &Haplotype,
					const TSimulatorReference &Reference, BAM::TOutputBamFile &BamFile);

	// getters
	const std::string& name() const noexcept { return _readGroup->name_ID; }
	[[nodiscard]] virtual double meanReadLength() const = 0;
	double maxFragmentLength() {
		return _fragmentLengthDistr.max();
	}
};

//-------------------------------
// TSimulatorSingleEndRead
//-------------------------------
class TReadSimulatorSingleEnd final : public TReadSimulator {
private:
	size_t _numCycles;

	void _simulate(const genometools::TGenomePosition & Position, const std::vector<genometools::TwoBase> & Haplotype, const TSimulatorReference &Reference) override;
	void _writeSimulatedAlignments(BAM::TOutputBamFile & BamFile) override;

public:
	TReadSimulatorSingleEnd(const BAM::TReadGroup & ReadGroup, const BAM::RGInfo::TReadGroupInfoEntry & RGInfo, const GenotypeLikelihoods::PMD::TModel & Pmd, const GenotypeLikelihoods::SequencingError::RGModels& Recal);
	~TReadSimulatorSingleEnd() = default;

	[[nodiscard]] double meanReadLength() const override;
};

//-------------------------------
// TSimulatorPairedEndReads
//-------------------------------
class TReadSimulatorPairedEnd final : public TReadSimulator {
private:
	BAM::TAlignment _mate;
	std::array<size_t, 2> _numCycles;

	void _simulate(const genometools::TGenomePosition & Position, const std::vector<genometools::TwoBase> & Haplotype, const TSimulatorReference &Reference) override;
	void _writeSimulatedAlignments(BAM::TOutputBamFile & BamFile) override;

public:
	TReadSimulatorPairedEnd(const BAM::TReadGroup & ReadGroup, const BAM::RGInfo::TReadGroupInfoEntry & RGInfo, const GenotypeLikelihoods::PMD::TModel & Pmd, const GenotypeLikelihoods::SequencingError::RGModels& Recal);
	~TReadSimulatorPairedEnd() = default;

	[[nodiscard]] double meanReadLength() const override;
};

}; // namespace Simulations

#endif /* TSIMULATORREAD_H_ */
