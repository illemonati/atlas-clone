/*
 * TRecalibrationEMModel.h
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMMODEL_H_
#define TRECALIBRATIONEMMODEL_H_

#include <limits>
#include <stdint.h>
#include <string>

#include "TEpsilon.h"
#include "TGenotypeData.h"
#include "TReadGroupInfo.h"
#include "TRho.h"
#include "TSequencedBase.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Types/probability.h"
#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"

namespace BAM {
class TSequencedBase;
}

namespace GenotypeLikelihoods {
namespace SequencingError {

//--------------------------------------------------------------------
// TModel
// pure abstract base class
//--------------------------------------------------------------------
struct TModel {
	virtual ~TModel()                                                                                 = default;
	virtual bool recalibrates() const noexcept                                                        = 0;
	virtual coretools::Probability errorRate(const BAM::TSequencedBase &base) const noexcept          = 0;
	virtual genometools::PhredIntProbability phredInt(const BAM::TSequencedBase &base) const noexcept = 0;
	virtual TBaseLikelihoods baseLikelihoods(const BAM::TSequencedBase &base) const noexcept          = 0;
	virtual void simulate(BAM::TSequencedBase &base) const noexcept                                   = 0;
	virtual std::string epsilonDefinition() const noexcept                                            = 0;
	virtual std::string rhoDefinition() const noexcept                                                = 0;
	virtual BAM::RGInfo::TInfo info() const                                                           = 0;
};

//------------------------------------------------
// TModelNoRecal
//------------------------------------------------
class TModelNoRecal final : public TModel {
public:
	bool recalibrates() const noexcept override { return false; };

	coretools::Probability errorRate(const BAM::TSequencedBase &base) const noexcept override;
	genometools::PhredIntProbability phredInt(const BAM::TSequencedBase &base) const noexcept override;
	TBaseLikelihoods baseLikelihoods(const BAM::TSequencedBase &base) const noexcept override;
	void simulate(BAM::TSequencedBase &base) const noexcept override;

	std::string epsilonDefinition() const noexcept override { return "-"; };
	std::string rhoDefinition() const noexcept override { return "-"; };
	BAM::RGInfo::TInfo info() const override { return BAM::RGInfo::TInfo{}; };
};

//------------------------------------------------
// TModelRecal
//------------------------------------------------
class TModelRecal final : public TModel {
private:
	TRho _rho;
	TEpsilon _epsilon;
public:
	TModelRecal(std::string_view EpsilonDef, std::string_view RhoDef = "default")
		: _rho(RhoDef), _epsilon(EpsilonDef) {}
	TModelRecal(const BAM::RGInfo::TInfo &info) : _rho(info["rho"]), _epsilon(info) {}

	// Virtual
	bool recalibrates() const noexcept override { return true; };

	// get error rates
	coretools::Probability errorRate(const BAM::TSequencedBase &base) const noexcept override;
	genometools::PhredIntProbability phredInt(const BAM::TSequencedBase &base) const noexcept override;
	TBaseLikelihoods baseLikelihoods(const BAM::TSequencedBase &base) const noexcept override;
	virtual void simulate(BAM::TSequencedBase &base) const noexcept override;

	std::string epsilonDefinition() const noexcept override {return _epsilon.definition();};
	std::string rhoDefinition() const noexcept override { return _rho.definition(); };
	BAM::RGInfo::TInfo info() const override;

	// Model-Access
	TRho& rho() noexcept {return _rho;}
	TEpsilon& epsilon() noexcept {return _epsilon;}
};

} // namespace SequencingError
} // namespace GenotypeLikelihoods

#endif /* TRECALIBRATIONEMMODEL_H_ */
