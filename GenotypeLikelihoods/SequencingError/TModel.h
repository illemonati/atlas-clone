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
	TRho(const std::string &Def);
	TRho(const TRho &other) = default;
	TRho &operator=(const TRho &other) = default;


	const TBaseProbabilities& operator[](genometools::Base from) const noexcept {
		return _rho[from];
	}
	std::string getDefinition() const noexcept;

	// functions used to estimate
	void add(genometools::Base l, coretools::Probability P_g_I_d, const TBaseProbabilities &P_bbar_I_d) noexcept;
	void estimate() noexcept;
};


//--------------------------------------------------------------------
// TModel
// pure abstract base class
//--------------------------------------------------------------------
class TModel {
public:
	virtual ~TModel()                                                                                    = default;
	virtual bool estimatable() const noexcept                                                            = 0;
	virtual bool recalibrates() const noexcept                                                           = 0;
	virtual coretools::Probability getErrorRate(const BAM::TSequencedBase &base) const noexcept          = 0;
	virtual genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase &base) const noexcept = 0;
	virtual TBaseLikelihoods getBaseLikelihoods(const BAM::TSequencedBase &base) const noexcept          = 0;
	virtual void simulate(BAM::TSequencedBase &base) const noexcept                                      = 0;
	virtual std::string getEpsilonDefinition() const noexcept                                          = 0;
	virtual std::string getRhoDefinition() const noexcept                                                = 0;
};

//------------------------------------------------
// TModelNoRecal
//------------------------------------------------
class TModelNoRecal final : public TModel {
public:
	bool estimatable() const noexcept override { return false; };
	bool recalibrates() const noexcept override { return false; };

	coretools::Probability getErrorRate(const BAM::TSequencedBase &base) const noexcept override;
	genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase &base) const noexcept override;
	TBaseLikelihoods getBaseLikelihoods(const BAM::TSequencedBase &base) const noexcept override;
	void simulate(BAM::TSequencedBase &base) const noexcept override;

	std::string getEpsilonDefinition() const noexcept override { return "-"; };
	std::string getRhoDefinition() const noexcept override { return "-"; };
};

//------------------------------------------------
// TModelRecal
//------------------------------------------------
class TModelRecal final : public TModel {
private:
	TRho _rho;
	TEpsilon _epsilon;
public:
	TModelRecal(const std::string &EpsilonDef, const std::string& RhoDef="default");
	TModelRecal(const BAM::RGInfo::TInfo & info) : _epsilon(""){};

	bool estimatable() const noexcept override { return true; };
	bool recalibrates() const noexcept override { return true; };
	std::string getEpsilonDefinition() const noexcept override {return _epsilon.getDefinition();};
	std::string getRhoDefinition() const noexcept override { return _rho.getDefinition(); };

	// get error rates
	coretools::Probability getErrorRate(const BAM::TSequencedBase &base) const noexcept override;
	genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase &base) const noexcept override;
	TBaseLikelihoods getBaseLikelihoods(const BAM::TSequencedBase &base) const noexcept override;
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
