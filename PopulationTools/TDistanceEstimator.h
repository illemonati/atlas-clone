/*
 * TDistanceCalculator.h
 *
 *  Created on: Jul 25, 2017
 *      Author: phaentu
 */

#ifndef TDISTANCEESTIMATOR_H_
#define TDISTANCEESTIMATOR_H_

#include <array>
#include <string>
#include <vector>

#include "coretools/Containers/TStrongArray.h"

#include "coretools/Files/TLineWriter.h"
#include "genometools/Genotypes/Base.h"
#include "genometools/Genotypes/Containers.h"

#include "genometools/GLF/GLF.h"
#include "genometools/GLF/TGLFVector.h"

namespace gz { class ogzstream; }

namespace PopulationTools{

//------------------------------------------------
// DistancePhi
//------------------------------------------------

enum class DistancePhi : uint8_t {min=0, aa_aa=min, aa_ab, ab_aa, aa_bb, ab_ab, ab_ac, aa_bc, ab_cc, ab_cd, nn_nn, max = nn_nn};

constexpr DistancePhi distancePhi(genometools::Genotype a, genometools::Genotype b) {
	using DP                                          = DistancePhi;
	constexpr std::array<std::array<DP, 11>, 11> dphi = {{
		// AA,	AC,		AG,		AT,		CC,		CG,		CT,		GG,		GT,		TT,		NN
		{DP::aa_aa, DP::aa_ab, DP::aa_ab, DP::aa_ab, DP::aa_bb, DP::aa_bc, DP::aa_bc, DP::aa_bb, DP::aa_bc, DP::aa_bb,
	     DP::nn_nn}, // AA
		{DP::ab_aa, DP::ab_ab, DP::ab_ac, DP::ab_ac, DP::ab_aa, DP::ab_ac, DP::ab_ac, DP::ab_cc, DP::ab_cd, DP::ab_cc,
	     DP::nn_nn}, // AC
		{DP::ab_aa, DP::ab_ac, DP::ab_ab, DP::ab_ac, DP::ab_cc, DP::ab_ac, DP::ab_cd, DP::ab_aa, DP::ab_ac, DP::ab_cc,
	     DP::nn_nn}, // AG
		{DP::ab_aa, DP::ab_ac, DP::ab_ac, DP::ab_ab, DP::ab_cc, DP::ab_cd, DP::ab_ac, DP::ab_cc, DP::ab_ac, DP::ab_aa,
	     DP::nn_nn}, // AT
		{DP::aa_bb, DP::aa_ab, DP::aa_bc, DP::aa_bc, DP::aa_aa, DP::aa_ab, DP::aa_ab, DP::aa_bb, DP::aa_bc, DP::aa_bb,
	     DP::nn_nn}, // CC
		{DP::ab_cc, DP::ab_ac, DP::ab_ac, DP::ab_cd, DP::ab_aa, DP::ab_ab, DP::ab_ac, DP::ab_aa, DP::ab_ac, DP::ab_cc,
	     DP::nn_nn}, // CG
		{DP::ab_cc, DP::ab_ac, DP::ab_cd, DP::ab_ac, DP::ab_aa, DP::ab_ac, DP::ab_ab, DP::ab_cc, DP::ab_ac, DP::ab_aa,
	     DP::nn_nn}, // CT
		{DP::aa_bb, DP::aa_bc, DP::aa_ab, DP::aa_bc, DP::aa_bb, DP::aa_ab, DP::aa_bc, DP::aa_aa, DP::aa_ab, DP::aa_bb,
	     DP::nn_nn}, // GG
		{DP::ab_cc, DP::ab_cd, DP::ab_ac, DP::ab_ac, DP::ab_cc, DP::ab_ac, DP::ab_ac, DP::ab_aa, DP::ab_ab, DP::ab_aa,
	     DP::nn_nn}, // GT
		{DP::aa_bb, DP::aa_bc, DP::aa_bc, DP::aa_ab, DP::aa_bb, DP::aa_bc, DP::aa_ab, DP::aa_bb, DP::aa_ab, DP::aa_aa,
	     DP::nn_nn}, // TT
		{DP::nn_nn, DP::nn_nn, DP::nn_nn, DP::nn_nn, DP::nn_nn, DP::nn_nn, DP::nn_nn, DP::nn_nn, DP::nn_nn, DP::nn_nn,
	     DP::nn_nn} // NN
	}};
	return dphi[coretools::index(a)][coretools::index(b)];
}

inline std::string toString(DistancePhi dp) {
	std::array strs = {"aa_aa", "aa_ab", "ab_aa", "aa_bb", "ab_ab", "ab_ac", "aa_bc", "ab_cc", "ab_cd", "nn_nn"};
	return strs[coretools::index(dp)];
}

//-------------------------------------
// TDistanceData
//-------------------------------------

using TDistanceData = coretools::TStrongArray<double, DistancePhi>;

//----------------------------------------------------
//TDistanceClass
//----------------------------------------------------
class TDistance{
protected:
	TDistanceData _distanceWeight;

public:
	TDistance();
	virtual ~TDistance() = default;

	[[nodiscard]] virtual double calculateDistance(const TDistanceData & phi);
	[[nodiscard]] const TDistanceData weights(){ return _distanceWeight; };
};

class TDistanceProbMismatch:public TDistance{
public:
	TDistanceProbMismatch();
};

class TDistanceEuclidian:public TDistance{
public:
	double calculateDistance(const TDistanceData & phi);
};

class TDistanceUser:public TDistance{
public:
	TDistanceUser(std::vector<double> vec);
};

//--------------------------------------------
//TDistanceEstimate
//--------------------------------------------

typedef std::vector<genometools::TGLFLikelihoods> GenotypeQualityVector;

class TEMforDistanceEstimation{
private:
	coretools::TStrongArray<coretools::TStrongArray<coretools::TStrongArray<bool, genometools::Base>, genometools::Genotype>, genometools::Genotype> _genoToBaseMap();

	//settings
	int _maxIterations;
	double _epsilon;

	//tmp variables
	double _old_LL;
	TDistanceData _K; //normalizing constant

	std::array< std::array<double, 10>, 10 > _probGeno;
	std::array< std::array<double, 10>, 10 > _P_G;
	std::array< std::array<double, 10>, 10 > _P_G_one_site;

	std::unique_ptr<TDistance> _distanceObject;

	genometools::TBaseData _pi;
	TDistanceData _phi;
	double _LL;
	double _distance;

	void _guessPi(GenotypeQualityVector & genoQual1, GenotypeQualityVector & genoQual2);
	void _guessPhi(GenotypeQualityVector & genoQual1, GenotypeQualityVector & genoQual2);
	void _fill_K(genometools::TBaseData  & thesePi);
	void _fill_P_g_given_phi_pi(const TDistanceData & phi, genometools::TBaseData & pi);

public:

	TEMforDistanceEstimation();
	bool estimatePhiWithEM(GenotypeQualityVector & genoQual1, GenotypeQualityVector & genoQual2);
	double distance() const noexcept {return _distance;}
	const genometools::TBaseData& pi() const noexcept {return _pi;};
	const TDistanceData& phi() const noexcept {return _phi;}
};

//--------------------------------------------
//TDistanceEstimator
//--------------------------------------------
class TDistanceEstimator{
private:
	int maxterations;
	double _epsilon;
	std::string _outputName;

	//GLF files
	genometools::TGLFVector _GLFs;

	void _estimateDistanceGenomeWide(TEMforDistanceEstimation & EM_object);
	void _moveToNextCommonChr(genometools::TGLFReader & g1, genometools::TGLFReader & g2);
	void _advance(genometools::TGLFReader & g1, genometools::TGLFReader & g2);
	void _readCommonSites(GenotypeQualityVector & genoQual1, GenotypeQualityVector & genoQual2, genometools::TGLFReader & g1, genometools::TGLFReader & g2);
	void _estimateDistanceGenomeWide(TEMforDistanceEstimation & EM_object, genometools::TGLFReader & g1, genometools::TGLFReader & g2, coretools::TLineWriter & out);

	void _estimateDistanceInWindows(TEMforDistanceEstimation & EM_object, uint32_t windowLen);
	void _estimateDistanceInWindows(TEMforDistanceEstimation & EM_object, std::string filename, genometools::TGLFReader & g1, genometools::TGLFReader & g2, uint32_t windowLen);

	void _writeDistanceEstimates(coretools::TLineWriter & out, const genometools::TChromosome& Chr, genometools::TGenomeWindow& Window, uint32_t numsitesWithData, TEMforDistanceEstimation & EM_object);
	void _writeDistanceEstimates(coretools::TLineWriter & out, int numsitesWithData, TEMforDistanceEstimation & EM_object);
	void _writeDistanceEstimatesNoData(coretools::TLineWriter & out, const genometools::TChromosome& Chr, genometools::TGenomeWindow& Window);
	void _writeDistanceEstimatesNoData(coretools::TLineWriter & out);

public:
	TDistanceEstimator();
	void run();
};

}; //end namespace

#endif /* TDISTANCEESTIMATOR_H_ */
