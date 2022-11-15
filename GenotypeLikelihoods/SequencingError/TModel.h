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

#include "TReadGroupInfo.h"
#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"
#include "TGenotypeData.h"
#include "TSequencedBase.h"
#include "coretools/Containers/TStrongArray.h"
#include "TEpsilon.h"
#include "coretools/Types/probability.h"

namespace BAM {
class TSequencedBase;
}
namespace GenotypeLikelihoods {
namespace RecalEstimatorTools {
class TRecalDataTable;
}
} // namespace GenotypeLikelihoods

namespace GenotypeLikelihoods {
namespace SequencingError {

//--------------------------------------------------------------------
// TRho
//--------------------------------------------------------------------
class TRho {
private:
	coretools::TStrongArray<TBaseProbabilities, genometools::Base> _rho{
		{TBaseProbabilities::normalize({0., 1. / 3, 1. / 3, 1. / 3}),
		 TBaseProbabilities::normalize({1. / 3, 0., 1. / 3, 1. / 3}),
		 TBaseProbabilities::normalize({1. / 3, 1. / 3, 0., 1. / 3}),
		 TBaseProbabilities::normalize({1. / 3, 1. / 3, 1. / 3, 0.})}}; //[from k][to l]

	coretools::TStrongArray<coretools::TStrongArray<double, genometools::Base>, genometools::Base> _rhoSum{};
public:
	TRho() = default;
	TRho(std::string_view Def);

	const TBaseProbabilities& operator[](genometools::Base from) const noexcept {
		return _rho[from];
	}
	std::string definition() const noexcept;

	// functions used to estimate
	void add(genometools::Base l, coretools::Probability P_g_I_d, const TBaseProbabilities &P_bbar_I_d) noexcept;
	void estimate() noexcept;

	BAM::RGInfo::TInfo info() const ;
};


//--------------------------------------------------------------------
// TModel
// pure abstract base class
//--------------------------------------------------------------------
class TModel {
public:
	virtual ~TModel()                                                                                 = default;
	virtual bool estimatable() const noexcept                                                         = 0;
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
	bool estimatable() const noexcept override { return false; };
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
	TModelRecal(const BAM::RGInfo::TInfo &) : _epsilon("") {}

	bool estimatable() const noexcept override { return true; };
	bool recalibrates() const noexcept override { return true; };
	std::string epsilonDefinition() const noexcept override {return _epsilon.definition();};
	std::string rhoDefinition() const noexcept override { return _rho.definition(); };
	BAM::RGInfo::TInfo info() const override;

	// get error rates
	coretools::Probability errorRate(const BAM::TSequencedBase &base) const noexcept override;
	genometools::PhredIntProbability phredInt(const BAM::TSequencedBase &base) const noexcept override;
	TBaseLikelihoods baseLikelihoods(const BAM::TSequencedBase &base) const noexcept override;
	virtual void simulate(BAM::TSequencedBase &base) const noexcept override;

	// functions to estimate
	void checkOrInit(const RecalEstimatorTools::TRecalDataTable &DataTable) { _epsilon.checkOrInit(DataTable); };

	// functions to estimate rho
	void addToRho(const BAM::TSequencedBase &data, coretools::Probability P_g_I_d, const TBaseProbabilities &P_bbar_I_d) noexcept {_rho.add(data.base, P_g_I_d, P_bbar_I_d);}; 
	void estimateRho() noexcept {_rho.estimate();};

	// functions to estimate betas
	double Q() const noexcept {return _epsilon.Q();}
	template<bool updateJF, bool isInvariant>
	void addToEpsilon(const BAM::TSequencedBase &base, const TGenotypeLikelihoods &P_g_I_ds, const TGenotypeLikelihoods &P_bbar_I_gds) {
		_epsilon.addToEpsilon<updateJF, isInvariant>(base, P_g_I_ds, P_bbar_I_gds);
	}
	void solveJxF() {_epsilon.solveJxF();}
	void propose(double lambda) {_epsilon.propose(lambda);}
	bool acceptOrReject() {return _epsilon.acceptOrReject();}
	void adjust() {_epsilon.adjust();}
	double maxF() const noexcept {return _epsilon.maxF();}
};

} // namespace SequencingError
} // namespace GenotypeLikelihoods

#endif /* TRECALIBRATIONEMMODEL_H_ */
