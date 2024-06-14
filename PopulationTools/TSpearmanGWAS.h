/*
 * THardyWeinbergTest.h
 *
 *  Created on: Jul 31, 2020
 *      Author: wegmannd
 */

#ifndef POPULATIONTOOLS_TSPEARMANGWAS_H_
#define POPULATIONTOOLS_TSPEARMANGWAS_H_

#include <array>
#include <map>
#include <string>
#include <vector>
#include <armadillo>

#include "genometools/GenotypeTypes.h"
#include "genometools/VCF/TPopulation.h"
#include "genometools/VCF/TVcfFile.h"
#include "genometools/VCF/TPopulationLikelihoods.h"

namespace coretools { class TOutputFile; }

namespace PopulationTools{

//------------------------------------------------
//TSpearmanGWASPopulation
//------------------------------------------------
class TSpearmanGWASPopulation{
private:
	
	std::vector<double> _data;
	std::vector<double> _dosage; //mean posterior genotype

	arma::vec _ranksData;	
	std::vector< std::vector<size_t> > _bootstraps;

	double _calcAbsRho(const double sumPairwiseProductDataDosage, const double productOfMeans, const double sqrtProductOfVariances) const;
	double _sumPairwiseProductBootstrap(const size_t b, const std::vector<double>& ranksDosage);
	
public:
	TSpearmanGWASPopulation() = default;
	TSpearmanGWASPopulation(size_t Size);	
	void resize(size_t Size);
	void clear();
	void addSample(double data);
	void finalizeDataReading();
	void prepareBootstraps(const size_t NumBootstraps);
		
	void updateDosage(genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> *GenotypeLikelihoods);
	//double calcSpearmanRhoAndBootstrap(std::vector<double> & bootstrapsRho);
	double regressSpearmanAndBootstrap(std::vector<double> & bootstrapsRho);
};


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

	//samples
    genometools::TPopulationSamples _samples;
	std::vector<TSpearmanGWASPopulation> _populations;

	//genotype data
	//THWPopulations _populations;

	void _openVCF();
	void _closeVCF();

	std::map<std::string, double> _readDataIntoMap(std::string_view Filename);

public:
	TSpearmanGWAS();
	void run();
};

} // namespace PopulationTools

#endif /* POPULATIONTOOLS_THARDYWEINBERGTEST_H_ */
