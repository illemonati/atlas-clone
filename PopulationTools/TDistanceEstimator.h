/*
 * TDistanceCalculator.h
 *
 *  Created on: Jul 25, 2017
 *      Author: phaentu
 */

#ifndef TDISTANCEESTIMATOR_H_
#define TDISTANCEESTIMATOR_H_

#include <stdint.h>
#include <algorithm>
#include <array>
#include <string>
#include <vector>
#include "GenotypeTypes.h"
#include "TGLF.h"
#include "TGenotypeData.h"
#include "TParameters.h"
#include "TTask.h"
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
	return dphi[genometools::index(a)][genometools::index(b)];
}

inline std::string toString(DistancePhi dp) {
	std::array strs = {"aa_aa", "aa_ab", "ab_aa", "aa_bb", "ab_ab", "ab_ac", "aa_bc", "ab_cc", "ab_cd", "nn_nn"};
	return strs[genometools::index(dp)];
}

//-------------------------------------
// TDistanceData
//-------------------------------------
class TDistanceData : public GenotypeLikelihoods::TData_base<double, DistancePhi, genometools::index(DistancePhi::max)>{
private:
	using GenotypeLikelihoods::TData_base<double, DistancePhi, genometools::index(DistancePhi::max)>::_data;

public:
	TDistanceData() : GenotypeLikelihoods::TData_base<double, DistancePhi, genometools::index(DistancePhi::max)>(0.0) {};
	~TDistanceData() = default;

	double& operator()(const genometools::Genotype & g1, const genometools::Genotype & g2){
		return _data[genometools::index(distancePhi(g1, g2))];
	};
};

//--------------------------------------------
//TGenocombinationToBaseMap
//--------------------------------------------
class TGenocombinationToBaseMap{
public:
	std::array< std::array< std::array< bool, 4>, 10>, 10> genotypeCombinationHasBase{};

	TGenocombinationToBaseMap();
	~TGenocombinationToBaseMap() = default;

	bool operator()(const genometools::Genotype &g1, const genometools::Genotype &g2, const genometools::Base &b) {
		using genometools::index;
		return genotypeCombinationHasBase[index(g1)][index(g2)][index(b)];
	};
};

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

typedef std::vector<GLF::TGLFLikelihoods> GenotypeQualityVector;

class TEMforDistanceEstimation{
private:
	TGenocombinationToBaseMap genoToBaseMap;

	//settings
	int maxNumEMIterations;
	double epsilonForEM;

	//tmp variables
	double old_LL;
	TDistanceData K; //normalizing constant

	std::array< std::array<double, 10>, 10 > probGeno;
	std::array< std::array<double, 10>, 10 > P_G;
	std::array< std::array<double, 10>, 10 > P_G_one_site;

	TDistance* distanceObject;

	void guessPi(GenotypeQualityVector & genoQual1, GenotypeQualityVector & genoQual2);
	void guessPhi(GenotypeQualityVector & genoQual1, GenotypeQualityVector & genoQual2);
	void fill_K(GenotypeLikelihoods::TBaseData  & thesePi);
	void fill_P_g_given_phi_pi(const TDistanceData & phi, GenotypeLikelihoods::TBaseData & pi);

public:
	GenotypeLikelihoods::TBaseData pi;
	TDistanceData phi;
	double LL;
	double distance;

	TEMforDistanceEstimation();
	~TEMforDistanceEstimation(){
		delete distanceObject;
	};

	bool estimatePhiWithEM(GenotypeQualityVector & genoQual1, GenotypeQualityVector & genoQual2);
};

//--------------------------------------------
//TDistanceEstimator
//--------------------------------------------
class TDistanceEstimator{
private:
	int maxNumEMIterations;
	double epsilonForEM;
	std::string outputName;

	//GLF files
	int numGLFs;
	std::vector<std::string> GLFNames;
	GLF::TGlfReader* glfs;
	bool readersOpened;

	void openGLF();
	void closeGLF();


	void estimateDistanceGenomeWide(TEMforDistanceEstimation & EM_object);
	bool moveToNextCommonChr(GLF::TGlfReader & g1, GLF::TGlfReader & g2);
	bool advance(GLF::TGlfReader & g1, GLF::TGlfReader & g2);
	void readCommonSites(GenotypeQualityVector & genoQual1, GenotypeQualityVector & genoQual2, GLF::TGlfReader & g1, GLF::TGlfReader & g2);
	void estimateDistanceGenomeWide(TEMforDistanceEstimation & EM_object, GLF::TGlfReader & g1, GLF::TGlfReader & g2, gz::ogzstream & out);

	void estimateDistanceInWindows(TEMforDistanceEstimation & EM_object, uint32_t windowLen);
	void estimateDistanceInWindows(TEMforDistanceEstimation & EM_object, std::string filename, GLF::TGlfReader & g1, GLF::TGlfReader & g2, uint32_t windowLen);

	void writeDistanceEstimates(gz::ogzstream & out, std::string & chr, uint32_t windowStart, uint32_t windowEnd, uint32_t numsitesWithData, TEMforDistanceEstimation & EM_object);
	void writeDistanceEstimates(gz::ogzstream & out, int numsitesWithData, TEMforDistanceEstimation & EM_object);
	void writeDistanceEstimatesNoData(gz::ogzstream & out, std::string & chr, uint32_t windowStart, uint32_t windowEnd);
	void writeDistanceEstimatesNoData(gz::ogzstream & out);

public:
	TDistanceEstimator();
	~TDistanceEstimator(){
		closeGLF();
	};

	void printGLF();
	void estimateDistances();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_estimateDist:public coretools::TTask{
public:
	TTask_estimateDist(){ _explanation = "Estimating the genetic distance between individuals"; };

	void run(){
		using namespace coretools::instances;
		TDistanceEstimator distEst;
		distEst.estimateDistances();
	};
};

}; //end namespace

#endif /* TDISTANCEESTIMATOR_H_ */
