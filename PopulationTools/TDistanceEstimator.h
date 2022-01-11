/*
 * TDistanceCalculator.h
 *
 *  Created on: Jul 25, 2017
 *      Author: phaentu
 */

#ifndef TDISTANCEESTIMATOR_H_
#define TDISTANCEESTIMATOR_H_

#include "TParameters.h"
#include "TGLF.h"
#include <math.h>
#include "GenotypeTypes.h"
#include "strongTypes.h"

namespace PopulationTools{

using namespace GenotypeLikelihoods;
using genometools::Genotype;
using genometools::Base;
using genometools::HighPrecisionPhredIntProbability;

//------------------------------------------------
// DistancePhi
//------------------------------------------------
enum DistancePhiEnum : uint8_t {aa_aa=0, aa_ab, ab_aa, aa_bb, ab_ab, ab_ac, aa_bc, ab_cc, ab_cd, nn_nn};

class DistancePhi : public coretools::StrongTypes::StrongType<DistancePhiEnum, DistancePhi> {
private:
    static std::string _toString[nn_nn + 1];
	static constexpr std::array< std::array<DistancePhiEnum, 11>, 11> _genoToPhiMap {{
    		//AA,	AC,		AG,		AT,		CC,		CG,		CT,		GG,		GT,		TT,		NN
    		{aa_aa,	aa_ab,	aa_ab,	aa_ab,	aa_bb,	aa_bc,	aa_bc,	aa_bb,	aa_bc,	aa_bb,	nn_nn}, //AA
    		{ab_aa,	ab_ab,	ab_ac,	ab_ac,	ab_aa,	ab_ac,	ab_ac,	ab_cc,	ab_cd,	ab_cc,	nn_nn}, //AC
    		{ab_aa,	ab_ac,	ab_ab,	ab_ac,	ab_cc,	ab_ac,	ab_cd,	ab_aa,	ab_ac,	ab_cc,	nn_nn}, //AG
    		{ab_aa,	ab_ac,	ab_ac,	ab_ab,	ab_cc,	ab_cd,	ab_ac,	ab_cc,	ab_ac,	ab_aa,	nn_nn}, //AT
    		{aa_bb,	aa_ab,	aa_bc,	aa_bc,	aa_aa,	aa_ab,	aa_ab,	aa_bb,	aa_bc,	aa_bb,	nn_nn}, //CC
    		{ab_cc,	ab_ac,	ab_ac,	ab_cd,	ab_aa,	ab_ab,	ab_ac,	ab_aa,	ab_ac,	ab_cc,	nn_nn}, //CG
    		{ab_cc,	ab_ac,	ab_cd,	ab_ac,	ab_aa,	ab_ac,	ab_ab,	ab_cc,	ab_ac,	ab_aa,	nn_nn}, //CT
    		{aa_bb,	aa_bc,	aa_ab,	aa_bc,	aa_bb,	aa_ab,	aa_bc,	aa_aa,	aa_ab,	aa_bb,	nn_nn}, //GG
    		{ab_cc,	ab_cd,	ab_ac,	ab_ac,	ab_cc,	ab_ac,	ab_ac,	ab_aa,	ab_ab,	ab_aa,	nn_nn}, //GT
    		{aa_bb,	aa_bc,	aa_bc,	aa_ab,	aa_bb,	aa_bc,	aa_ab,	aa_bb,	aa_ab,	aa_aa,	nn_nn}, //TT
    		{nn_nn,	nn_nn,	nn_nn,	nn_nn,	nn_nn,	nn_nn,	nn_nn,	nn_nn,	nn_nn,	nn_nn,	nn_nn}  //NN
		}};

public:
    //constructors
    explicit constexpr DistancePhi() : StrongType(nn_nn) {};

    explicit constexpr DistancePhi(const DistancePhiEnum & Phi) : StrongType(Phi) {};

    explicit constexpr DistancePhi(const Genotype & First, const Genotype & Second){
        _value = _genoToPhiMap[First.get()][Second.get()];
    };

    //assignments
    constexpr void operator=(const DistancePhiEnum & Phi){
        _value = Phi;
    };

    void set(const Genotype & First, const Genotype & Second){
        _value = _genoToPhiMap[First.get()][Second.get()];
    };

    //convert
    [[nodiscard]] explicit operator std::string() const {
        return _toString[static_cast<uint8_t>( _value)];
    };

    //++operator, <operator and range info (to loop)
    [[nodiscard]] static constexpr DistancePhi min() { return DistancePhi(aa_aa); };
    [[nodiscard]] static constexpr DistancePhi max() { return DistancePhi(nn_nn); };

    constexpr DistancePhi& operator++(){
        if(_value == nn_nn){
            throw std::runtime_error("constexpr DistancePhi& operator++(): overflow!");
        };
        _value = static_cast<DistancePhiEnum>( static_cast<uint8_t>(_value) + 1);
        return *this;
    };

    constexpr bool operator<(const DistancePhi & other){
        return static_cast<uint8_t>(_value) < static_cast<uint8_t>(other.get());
    };
};

std::ostream& operator<<(std::ostream& os, const DistancePhi & Phi);

//-------------------------------------
// TDistanceData
//-------------------------------------
class TDistanceData : public TData_base<double, DistancePhi, DistancePhiEnum, nn_nn>{
private:
	using TData_base<double, DistancePhi, DistancePhiEnum, nn_nn>::_data;

public:
	TDistanceData() : TData_base<double, DistancePhi, DistancePhiEnum, nn_nn>(0.0) {};
	~TDistanceData() = default;

	double& operator()(const Genotype & g1, const Genotype & g2){
		return _data[DistancePhi(g1, g2).get()];
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

	bool operator()(const Genotype & g1,const  Genotype & g2, const Base & b){
		return genotypeCombinationHasBase[g1.get()][g2.get()][b.get()];
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
class TEMforDistanceEstimation{
private:
	coretools::TLog* logfile;
	TGenocombinationToBaseMap genoToBaseMap;

	//settings
	int maxNumEMIterations;
	double epsilonForEM;

	//tmp variables
	double old_LL;
	TDistanceData K; //normalizing constant

	double** probGeno;
	double** P_G;
	double** P_G_one_site;
//	double* distanceWeight; //weight for each phi class towards the distance.
	TDistance* distanceObject;

//	void calculateDistance();
	void guessPi(std::vector<HighPrecisionPhredIntProbability*> & genoQual1, std::vector<HighPrecisionPhredIntProbability*> & genoQual2);
	void guessPhi(std::vector<HighPrecisionPhredIntProbability*> & genoQual1, std::vector<HighPrecisionPhredIntProbability*> & genoQual2);
	void fill_K(GenotypeLikelihoods::TBaseData  & thesePi);
	void fill_P_g_given_phi_pi(const TDistanceData & phi, GenotypeLikelihoods::TBaseData & pi);

public:
	GenotypeLikelihoods::TBaseData pi;
	TDistanceData phi;
	double LL;
	double distance;

	TEMforDistanceEstimation(coretools::TLog* Logfile, coretools::TParameters & params);
	~TEMforDistanceEstimation(){
		for(int g1=0; g1<10; ++g1){
			delete[] probGeno[g1];
			delete[] P_G[g1];
			delete[] P_G_one_site[g1];
		}
		delete[] probGeno;
		delete[] P_G;
		delete[] P_G_one_site;
		delete distanceObject;
//		delete[] distanceWeight;
	};

	bool estimatePhiWithEM(std::vector<HighPrecisionPhredIntProbability*> & genoQual1, std::vector<HighPrecisionPhredIntProbability*> & genoQual2);
};

//--------------------------------------------
//TDistanceEstimator
//--------------------------------------------
class TDistanceEstimator{
private:
	coretools::TLog* logfile;
	int maxNumEMIterations;
	double epsilonForEM;
	std::string outputName;

	//GLF files
	int numGLFs;
	std::vector<std::string> GLFNames;
	GLF::TGlfReader* glfs;
	bool readersOpened;

	void openGLF(coretools::TParameters & params);
	void closeGLF();


	void estimateDistanceGenomeWide(TEMforDistanceEstimation & EM_object);
	bool moveToNextCommonChr(GLF::TGlfReader & g1, GLF::TGlfReader & g2);
	bool advance(GLF::TGlfReader & g1, GLF::TGlfReader & g2);
	void readCommonSites(std::vector<HighPrecisionPhredIntProbability*> & genoQual1, std::vector<HighPrecisionPhredIntProbability*> & genoQual2, GLF::TGlfReader & g1, GLF::TGlfReader & g2);
	void estimateDistanceGenomeWide(TEMforDistanceEstimation & EM_object, GLF::TGlfReader & g1, GLF::TGlfReader & g2, gz::ogzstream & out);

	void estimateDistanceInWindows(TEMforDistanceEstimation & EM_object, uint32_t windowLen);
	void estimateDistanceInWindows(TEMforDistanceEstimation & EM_object, std::string filename, GLF::TGlfReader & g1, GLF::TGlfReader & g2, uint32_t windowLen);

	void writeDistanceEstimates(gz::ogzstream & out, std::string & chr, uint32_t windowStart, uint32_t windowEnd, uint32_t numsitesWithData, TEMforDistanceEstimation & EM_object);
	void writeDistanceEstimates(gz::ogzstream & out, int numsitesWithData, TEMforDistanceEstimation & EM_object);
	void writeDistanceEstimatesNoData(gz::ogzstream & out, std::string & chr, uint32_t windowStart, uint32_t windowEnd);
	void writeDistanceEstimatesNoData(gz::ogzstream & out);

public:
	TDistanceEstimator(coretools::TLog* Logfile, coretools::TParameters & params);
	~TDistanceEstimator(){
		closeGLF();
	};

	void printGLF(coretools::TParameters & params);
	void estimateDistances(coretools::TParameters & params);
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_estimateDist:public coretools::TTask{
public:
	TTask_estimateDist(){ _explanation = "Estimating the genetic distance between individuals"; };

	void run(coretools::TParameters & Parameters, coretools::TLog* Logfile){
		TDistanceEstimator distEst(Logfile, Parameters);
		distEst.estimateDistances(Parameters);
	};
};

}; //end namespace

#endif /* TDISTANCEESTIMATOR_H_ */
