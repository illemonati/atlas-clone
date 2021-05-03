/*
 * TRecalibrationEM.h
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMESTIMATOR_H_
#define TRECALIBRATIONEMESTIMATOR_H_

#include "auxiliaryTools.h"
#include "TPostMortemDamage.h"
#include "TSequencingErrorModels.h"
#include "TSite.h"
#include "TGenotypeDistribution.h"
#include "TGenotypeLikelihoodCalculator.h"

namespace GenotypeLikelihoods{

namespace RecalEstimatorTools {

//--------------------------------------------------------------------
// TRecalibrationEMDataTable
// Object to store for which qualities and positions data is available.
//--------------------------------------------------------------------
class TDataTable{
private:
	uint64_t counts;
	bool initialized;

public:
	int maxQual, maxFragmentLength, maxMQ;
	uint16_t maxPos;
	unsigned int* qualities;
	unsigned int* fragmentLengths;
	unsigned int* MQ;

	TDataTable();
	TDataTable(const int MaxQual, const int MaxFragmentLength, const int MQ);
	~TDataTable();

	void initialize(const int MaxQual, const int MaxFragmentLength, const int MQ);
	void clear();
	void add(const BAM::TBase & base);
	void add(const TSite & site);
	size_t size() const;
	void fillVectorWithUsedQualities(std::vector<uint16_t> & Q) const;
	void fillVectorWithUsedFragmentLengths(std::vector<uint16_t> & lengths) const;
	void fillVectorWithUsedMQ(std::vector<uint16_t> & MQ) const;
};

class TDataTables{
private:
	uint16_t _numReadGroups;
	uint16_t _maxQual;
	TDataTable** _tables; //tables[readGroup][first/second]
	uint64_t _totalCounts;
	bool _initialized;

public:
	TDataTables();
	TDataTables(const int NumReadGroups, const int MaxQual, const int MaxFragmentLength, const int MQ);
	~TDataTables();

	void clear();
	void init(const int NumReadGroups, const int MaxQual, const int MaxFragmentLength, const int MaxMQ);
	void reset();
	void add(const BAM::TBase & base);
	void add(const TSite & site);
	void fillVectorWithUsedQualities(const int readGroupId, const bool isSecondMate, std::vector<uint16_t> & Q) const;

	uint64_t size() const;
	const TDataTable& table(const int readGroupId, const bool isSecondMate) const;
	const TDataTable& table(const int readGroupId, const int isSecondMate) const;
};

//------------------------------------------------
// Classes to keep track of models to estimate
//------------------------------------------------
class TModelStatusEntry{
private:
	bool _first;
	bool _second;

public:
	TModelStatusEntry(){
		_first = false;
		_second = false;
	};

	uint16_t size(){
		return _first + _second;
	};

	void set(const bool & IsSecondMate){
		if(!IsSecondMate){
			_first = true;
		} else {
			_second = true;
		}
	};

	std::string getString() const{
		if(_first && _second){
			return "(first and second mates)";
		} else if(_first){
			return "(first mate)";
		} else if(_second){
			return "(second mate)";
		} else {
			return "none";
		}
	}
};

enum ModelStatusTypes : uint8_t {copied=0, noData, littleData, dataButNoRecal};

class TModelStatus{
private:
	std::array<TModelStatusEntry, 4> _status;

public:
	TModelStatusEntry& operator[](const ModelStatusTypes & Type){
		return _status[Type];
	};
};

class TModelStati{
private:
	std::map<uint16_t, TModelStatus> modelStatus;

public:

	void add(const uint16_t & ReadGroupId){
		modelStatus.emplace(ReadGroupId);
	};

	TModelStatus& operator[](const uint16_t & ReadGroupId){
		return modelStatus[ReadGroupId];
	};

	uint16_t num(const ModelStatusTypes & Type){
		uint16_t num = 0;
		for(auto& m : modelStatus){
			m.second[Type].size();
		}
		return num;
	};

	void report(const ModelStatusTypes & Type, const std::string & Title, const BAM::TReadGroups & ReadGroups, TLog* Logfile){
		if(num(Type) > 0){
			Logfile->startIndent(Title);
			for(auto& m : modelStatus){
				if(m.second[Type].size() > 0){
					Logfile->list(ReadGroups.getName(m.first), " ", m.second[Type].getString());
				}
			}
			Logfile->endIndent();
		}
	}
};

//--------------------------------------------------------------------------------------
// TModelIndex
// Object to map read group ID and mate to an internal index used to store recal models
//--------------------------------------------------------------------------------------
class TModelIndex{
private:
	std::vector< std::array< std::shared_ptr<TSequencingErrorModelRecal> , 2> > _index;

public:
	TModelIndex(const BAM::TReadGroups& ReadGroups){
		_index.resize(ReadGroups.size());
	};
	~TModelIndex() = default;

	void set(const uint16_t & ReadGroupId, const bool & IsSecondMate, std::shared_ptr<TSequencingErrorModelRecal> & Model, const BAM::TReadGroupMap & ReadGroupMap){
		 _index[ReadGroupId][IsSecondMate] = Model;
	};

	bool inUse(const BAM::TBase & base) const{
		return _index[base.readGroupID][base.isSecondMate()];
	};

	std::shared_ptr<TSequencingErrorModelRecal>& operator()(const BAM::TBase & base) const{
		return _index[base.readGroupID][base.isSecondMate()];
	};
};

}; //end namespace RecalEstimatorTools

//--------------------------------------------------------------------------
// TSequencingErrorModelVectorForEstimation
// store pointers to models so that they can be estimated
//--------------------------------------------------------------------------

class TSequencingErrorModelVectorForEstimation{
private:
	//vector of pointers to models that require estimation
	std::vector< std::shared_ptr<TSequencingErrorModelRecal> > _models;
	RecalEstimatorTools::TModelIndex _modelIndex;

public:
	TSequencingErrorModelVectorForEstimation::TSequencingErrorModelVectorForEstimation(TSequencingErrorModels & SequencingErrorModels,
																					   const RecalEstimatorTools::TDataTables & DataTables,
																					   const BAM::TReadGroups & ReadGroups,
																					   const BAM::TReadGroupMap & ReadGroupMap,
																					   const uint32_t & MinRequiredObservations,
																					   TLog* Logfile);
	~TSequencingErrorModelVectorForEstimation() = default;


	//functions to estimate rho
	void prepareRhoEstimationFromEMWeights();
	void addBaseForRhoEstimation(BAM::TBase & base, const TBaseData & EMWeights);
	void estimateRho();

	//functions to estimate beta
	void setNewtonRaphsonParamsToZero();

	void addToFandJacobian(const BAM::TBase & base, const TBaseData & EM_weights_bbar_given_d);
	void addToQ(const BAM::TBase & base, const TBaseData & EM_weights_bbar_given_d);

	void setQToZero();

	double curQ();
	bool solveJxF();
	void proposeNewParameters(double lambda);
	void scaleParameters();
	unsigned int acceptProposedParametersBasedOnQ();
	void adjustParametersPostEstimation();
	double getSteepestGradient();
};

//--------------------------------------------------------------------
// TRecalibrationEMEstimator
//--------------------------------------------------------------------
class TRecalibrationEMEstimator:public TGenotypeLikelihoodCalculator{
private:
	std::vector<TSite> _sites;
	std::unique_ptr<TGenotypeDistribution> _genoDist;
	RecalEstimatorTools::TDataTables _dataTables;

	//variables for estimation
	bool equalBaseFrequencies;
	int _numEMIterations;
	double _minDeltaLL;
	int _NewtonRaphsonNumIterations;
	double _NewtonRaphsonMaxF;
	unsigned int _minRequiredObservations;
	std::string _recalFile; //file name in case a file with model is provided
	TSequencingErrorModelDefinition _modelDefitionForEstimation;

	void _initializeModels();
	void _runEM(std::string outputName, bool & writeTmpTables);
	void _fillRelevantBaseFrequencies(TBaseData & baseFreq, const Genotype genotype);

	//functions to estimate theta_epsilon (sequencing error rates)
	void _calculate_EMWeights_epsilon(std::vector<TBaseData> & EMWeights);
	double _calculate_Q_beta(const std::vector<TBaseData> & EM_weights_bbar_given_d);
	void _calculate_J_F_beta(const std::vector<TBaseData> & EM_weights_bbar_given_d);
	void _updateEM_theta_epsilon();

	//function to calculate LL (currently uses haploid model)
	double _calculateLL_fullModel();

public:
	TRecalibrationEMEstimator(TParameters & args, BAM::TReadGroups* ReadGroups, TLog* Logfile, BAM::TReadGroupMap* ReadGroupMap);

	//functions to add data
	void addSite(const TSite & site);
	size_t numSites();
	size_t numSitesDepthTwoOrMore();
	void compileDataTable(TDataTables & dataTable);
	size_t cumulativeDepth();

	//function to estimate
	void initializeFromFile(const std::string string);
	void performEstimation(std::string outputName, bool & writeTmpTables);
	void performEstimationKnownGenotypes(std::string outputName, bool & writeTmpTables);

	void writeCurrentEstimates(std::string filename);
	double calcLL();
	void calcLikelihoodSurface(std::string filename, int numMarginalGridPoints);
	void calcQSurface(std::string filename, int numMarginalGridPoints);
};

}; //end namespace

#endif /* TRECALIBRATIONEMESTIMATOR_H_ */
