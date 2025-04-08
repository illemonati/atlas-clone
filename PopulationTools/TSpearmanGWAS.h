/*
 * THardyWeinbergTest.h
 *
 *  Created on: Jul 31, 2020
 *      Author: wegmannd
 */

#ifndef POPULATIONTOOLS_TSPEARMANGWAS_H_
#define POPULATIONTOOLS_TSPEARMANGWAS_H_

#include <map>
#include <string>
#include <vector>
#include <armadillo>

#include "genometools/TSampleLikelihoods.h"
#include "genometools/VCF/TPopulation.h"
#include "genometools/VCF/TVcfFile.h"

namespace coretools { class TOutputFile; }

namespace PopulationTools{

namespace SpearmanGWASimpl{

//--------------------------------------------
// TSpearmanGWASBootstraps
//--------------------------------------------
class TSpearmanGWASBootstraps{
private:
	size_t _sampleSize{0};
	std::vector< std::vector<size_t> > _bootstraps;

public:
	TSpearmanGWASBootstraps() = default;
	TSpearmanGWASBootstraps(size_t SampleSize, size_t NumBootstraps);

	size_t sampleSize() const { return _sampleSize; }
	size_t getBootstrappedIndex(size_t Bootstrap, size_t Index) const { return _bootstraps[Bootstrap][Index]; }
	const std::vector<size_t>& get(size_t Bootstrap) const { return _bootstraps[Bootstrap]; } 
};

class TSpearmanGWASBootstrapLibrary{
private:
	size_t _numBootstraps;
	std::map<size_t, TSpearmanGWASBootstraps> _bootstraps;

public:
	TSpearmanGWASBootstrapLibrary(size_t NumBootstraps);

	size_t numBootstraps() const { return _numBootstraps; }
	const TSpearmanGWASBootstraps& get(size_t SampleSize);
};

//------------------------------------------------
//TSpearmanGWASPopulation
//------------------------------------------------
class TSpearmanGWASPopulation{
private:
	
	std::vector<double> _data;
	std::vector<double> _dataWithData; // vector only containing the data of individuals with genotypic information
	std::vector<double> _dosage; //mean posterior genotype
	std::vector<double> _ranksDosage;
	size_t _sampleSize;

	std::vector<double> _ranksData;	
	double _RSStotal;
	std::vector<double> _XT_X_inv_X;

	double _calcAbsRho(const double sumPairwiseProductDataDosage, const double productOfMeans, const double sqrtProductOfVariances) const;
	double _sumPairwiseProductBootstrap(const size_t b, const std::vector<double>& ranksDosage);
	
	void _regressSpearmanAndAddRSSBootstrap(double &RSS_null, double &RSS_model, const std::vector<size_t> &Bootstraps);

public:
	void clear();
	void addSample(double data);
	void finalizeDataReading();
		
	void updateDosage(genometools::TSampleLikelihoods<coretools::HPPhredInt> *GenotypeLikelihoods);
	size_t sampleSize(){ return _sampleSize; }
	void regressSpearmanAndAddRSS(double& RSS_null, double& RSS_model);
	void bootstrap(std::vector<double> &RSS_null, std::vector<double> &RSS_model,
	               TSpearmanGWASBootstrapLibrary &bootstrapLib);
};

} // end namesapce impl

//------------------------------------------------
//TSpearmanGWAS
//------------------------------------------------
class TSpearmanGWAS{
private:
	std::string _outname;

	//vcf-file
	std::string _vcfFilename;
    genometools::TVcfFileSingleLine _vcfFile;
	bool _limitLines;
	size_t _maxNumLines;
	size_t _numBootstraps;
	double _bootstrapPThreshold;

	//samples
    genometools::TPopulationSamples _samples;
	std::vector<SpearmanGWASimpl::TSpearmanGWASPopulation> _populations;

	void _openVCF();
	void _closeVCF();

	std::map<std::string, double> _readDataIntoMap(std::string_view Filename);
	double _calculateF(const double RSS_null, const double RSS_model, const double F_scale);

public:
	TSpearmanGWAS();
	
	void run();
};

} // namespace PopulationTools

#endif /* POPULATIONTOOLS_THARDYWEINBERGTEST_H_ */
