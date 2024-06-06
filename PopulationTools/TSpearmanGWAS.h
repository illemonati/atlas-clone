/*
 * THardyWeinbergTest.h
 *
 *  Created on: Jul 31, 2020
 *      Author: wegmannd
 */

#ifndef POPULATIONTOOLS_THARDYWEINBERGTEST_H_
#define POPULATIONTOOLS_THARDYWEINBERGTEST_H_

#include <array>
#include <map>
#include <string>
#include <vector>

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

	std::vector<double> _ranksData;	
	std::pair<double, double> _meanVarData;		

	std::vector< std::vector<size_t> > _bootstraps;
	std::vector<double> _bootstrapsRho;

	double _calcRho(const double sumPairwiseProductDataDosage, const double productOfMeans, const double sqrtProductOfVariances) const;
	double _sumPairwiseProductBootstrap(const size_t b);
	
public:
	TSpearmanGWASPopulation() = default;
	TSpearmanGWASPopulation(size_t Size);	
	void resize(size_t Size);
	void clear();
	void addSample(double data);
	void finalizeDataReading();
	void prepareBootstraps(const size_t NumBootstraps);
		
	void updateDosage(genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> *GenotypeLikelihoods);
	double calcSpearmanRhoAndBootstrap(std::vector<double> & _bootstrapsRho);
	double bootstrappedRho(const size_t BootstrapIndex) const { return _bootstrapsRho[BootstrapIndex]; }
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
