#include "TSpearmanGWAS.h"

#include "coretools/Files/TOutputFile.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Math/mathFunctions.h"
#include "coretools/Main/TRandomGenerator.h"

#include "genometools/DefaultColumnNames.h"
#include "genometools/TSampleLikelihoods.h"


#include <map>

namespace PopulationTools{

using coretools::TOutputFile;
using coretools::TInputFile;
using coretools::instances::logfile;
using coretools::instances::parameters;

//------------------------------------------------
//TSpearmanGWASPopulation
//------------------------------------------------

void TSpearmanGWASPopulation::resize(size_t Size){
	clear();
	_data.reserve(Size);
	_dosage.reserve(Size);
}

void TSpearmanGWASPopulation::clear(){
	_data.clear();
	_dosage.clear();
}

void TSpearmanGWASPopulation::addSample(double Data) {
	_data.emplace_back(Data);
}

void TSpearmanGWASPopulation::finalizeDataReading(){
	_ranksData = coretools::ranks(_data);
	_meanVarRanksData = coretools::meanVar(_ranksData);	
	_dosage.resize(_data.size());
}

void TSpearmanGWASPopulation::prepareBootstraps(const size_t NumBootstraps){
	_bootstraps.resize(NumBootstraps);
	for(size_t i = 0; i < NumBootstraps; ++i){
		_bootstraps[i].resize(_data.size());
		coretools::instances::randomGenerator().fillRandomOrderOfIndexes(_bootstraps[i]);
	}
}

void TSpearmanGWASPopulation::updateDosage(genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> *GenotypeLikelihoods){
	// calculate dosage for each sample 
	for(size_t i = 0; i < _data.size(); ++i){
		_dosage[i] = GenotypeLikelihoods[i].meanPosteriorGenotype();
	}
}

double TSpearmanGWASPopulation::_calcAbsRho(const double sumPairwiseProductDataDosage, const double productOfMeans, const double sqrtProductOfVariances) const {
	double E_XY = sumPairwiseProductDataDosage / _data.size();
	return std::fabs( (E_XY - productOfMeans) / sqrtProductOfVariances );
}

double TSpearmanGWASPopulation::_sumPairwiseProductBootstrap(const size_t b, const std::vector<double>& ranksDosage){
	double sum = 0.0;
	for(size_t i = 0; i < _data.size(); ++i){
		sum += _ranksData[i] * ranksDosage[ _bootstraps[b][i] ];
	}
	return(sum);
}

double TSpearmanGWASPopulation::calcSpearmanRhoAndBootstrap(std::vector<double> & bootstrapsRho){
	auto ranksDosage = coretools::ranks(_dosage);
	auto meanVarRanksDosage = coretools::meanVar(ranksDosage);

	double productOfMeans = _meanVarRanksData.first * meanVarRanksDosage.first;
	double sqrtProductOfVariances = sqrt(_meanVarRanksData.second * meanVarRanksDosage.second);

	// fill bootstraps
	for(size_t b = 0; b < _bootstraps.size(); ++b){
		bootstrapsRho[b] += _calcAbsRho(_sumPairwiseProductBootstrap(b, ranksDosage), productOfMeans, sqrtProductOfVariances);	
	}
	return _calcAbsRho(coretools::sumPairwiseProduct(_ranksData, ranksDosage), productOfMeans, sqrtProductOfVariances);
}

//------------------------------------------------
//THardyWeinbergTest
//------------------------------------------------

// read data
std::map<std::string, double> TSpearmanGWAS::_readDataIntoMap(std::string_view Filename){	
	TInputFile in(Filename, coretools::FileType::Header);

	// get sample and data columns
	size_t sampleIdx = in.indexOfFirstMatch(genometools::defaultColumnNames_sample);
	int dataIdx = -1;
	if(parameters().exists("dataCol")){
		// user provided column name
		auto dataColName = parameters().get<std::string>("dataCol");
		logfile().list("Reading data from column '", dataColName, "' of file '", Filename, "'. (parameter 'dataCol')");		
		dataIdx = in.index(dataColName);
	} else {
		// use first column in file that is not sample nor population
		for(size_t i = 0; i < in.numCols(); ++i){
			if(std::find(genometools::defaultColumnNames_sample.begin(), genometools::defaultColumnNames_sample.end(), in.header()[i]) == genometools::defaultColumnNames_sample.end()
			&& std::find(genometools::defaultColumnNames_population.begin(), genometools::defaultColumnNames_population.end(), in.header()[i]) == genometools::defaultColumnNames_population.end()){
				//column is not a sample or population
				dataIdx = i;								
				logfile().list("Reading data from column '", in.header()[i], "' of file '", Filename, "'. (specify column with 'dataCol') ");
				break;				
			}			
		}

		// if no column was found: error!
		if(dataIdx < 0){
			UERROR("Sample file '", Filename, "' lacks a data column!");
		}
	}
	
	// read into map
	std::map<std::string, double> data;
	for(; !in.empty(); in.popFront()){
		data.emplace(in.get(sampleIdx), in.get<double>(dataIdx));
	}

	return data;
}


TSpearmanGWAS::TSpearmanGWAS(){	
	//read samples
	std::string sampleFile = parameters().get<std::string>("samples");
	_samples.readSamples(sampleFile);	

	//open VCF
	_vcfFilename = parameters().get<std::string>("vcf");
	logfile().list("Reading vcf from file '" + _vcfFilename + "'.");
	_vcfFile.openStream(_vcfFilename);

	//enable parsers
	_vcfFile.enablePositionParsing();
	_vcfFile.enableVariantParsing();
	_vcfFile.enableVariantQualityParsing();
	_vcfFile.enableFormatParsing();
	_vcfFile.enableSampleParsing();

	//Match samples
	_samples.fillVCFOrder(_vcfFile.parser.samples);		

	// read data into map
	auto data = _readDataIntoMap(sampleFile);
		
	// create populations and add samples in correct order
	_populations.resize(_samples.numPopulations());		
	for(size_t s = 0; s < _samples.numSamples(); ++s){
		_populations[_samples.populationIndex(s)].addSample(data[_samples.sampleName(s)]);
	}

	const auto tmp     = coretools::str::readBeforeLast(_vcfFilename, ".vcf");
	auto tag = parameters().get("out", tmp);
	_outname = tag + "_SpearmanGWAS.txt.gz";
	logfile().list("Writing Spearman GWAS results to file '" + _outname + "'. (parameter 'out')");

	//limit lines?
	_limitLines = parameters().exists("limitLines");
	if(_limitLines){
		_maxNumLines = parameters().get<long>("limitLines");
	} else {
		_maxNumLines = 0;
	}

	// num bootstraps
	_numBootstraps = parameters().get<int>("bootstraps", 1000);
	logfile().listFlush("Preparing ", _numBootstraps, " bootstraps ... ");
	for(auto& p: _populations){
		p.finalizeDataReading();
		p.prepareBootstraps(_numBootstraps);
	}
	logfile().write("done! (parameter 'bootstraps')");
};

void TSpearmanGWAS::run(){
	// open VCF reader	
	genometools::TPopulationLikelihoodReaderLocus reader(false);
	reader.openVCF(_vcfFilename);
	logfile().endIndent();

	//open output file	
	std::vector<std::string> header = {genometools::defaultColumnNames_chromosome[0], genometools::defaultColumnNames_position[0], genometools::defaultColumnNames_reference[0], genometools::defaultColumnNames_alternative[0], "rho", "p"};	
	TOutputFile out(_outname, header);	

	//progress
	coretools::TTimer timer;	

	// likelihood and tmp storage
	genometools::TPopulationLikehoodLocus< genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability > > data(_samples.numSamples());
	std::vector<double> bootstrappedRho(_numBootstraps);

	// run through VCF file
	logfile().startIndent("Parsing VCF file and perform Spearman GWAS:");
	while (reader.readDataFromVCF(data, _samples)) {
		// update dosage for each sample in populations, calculate rho and bootstrap
		double sumRho = 0.0;
		std::fill(bootstrappedRho.begin(), bootstrappedRho.end(), 0.0);
	
		for(size_t p = 0; p < _populations.size(); ++p){
			_populations[p].updateDosage(&data[_samples.startIndex(p)]);
			sumRho += _populations[p].calcSpearmanRhoAndBootstrap(bootstrappedRho);
		}
		
		// determine p value from bootstraps
		size_t nBelow = 0;
		for(size_t b = 0; b < _numBootstraps; ++b){
			if(bootstrappedRho[b] < sumRho){
				++nBelow;
			}
		}
		double p = (double) nBelow / (double) _numBootstraps;

		// output
		reader.writePosition(out);
		out.write(sumRho, p);
		out.endln();		
		if(_limitLines && out.curLine() >= _maxNumLines){
			break;
		}		
	}

	// report final status
	logfile().endIndent();
	reader.concludeFilters();
	logfile().endIndent();
};



}; // namespace PopulationTools
