/*
 * THardyWeinbergTest.cpp
 *
 *  Created on: Jul 31, 2020
 *      Author: wegmannd
 */

#include "THardyWeinbergTest.h"

namespace PopulationTools{

//------------------------------------------------
//THardyWeinbergTest
//------------------------------------------------
THardyWeinbergTest::THardyWeinbergTest(TParameters & Parameters, TLog* logfile, TRandomGenerator* RandonGenerator){
	_logfile = logfile;
	_randonGenerator = RandonGenerator;

	//read samples
	if(Parameters.parameterExists("samples")){
		_samples.readSamples(Parameters.getParameterString("samples"), logfile);
	}

	//open VCF
	_vcfFilename = Parameters.getParameterString("vcf");
	_logfile->startIndent("Reading vcf from file '" + _vcfFilename + "'.");
	_vcfFile.openStream(_vcfFilename);

	//enable parsers
	_vcfFile.enablePositionParsing();
	_vcfFile.enableVariantParsing();
	_vcfFile.enableVariantQualityParsing();
	_vcfFile.enableFormatParsing();
	_vcfFile.enableSampleParsing();

	//Match samples
	if(_samples.hasSamples()){
		_samples.fillVCFOrder(_vcfFile.parser.samples);
	} else {
		_samples.readSamplesFromVCFNames(_vcfFile.parser.samples);
	}
};

void THardyWeinbergTest::testForHardyWeinberg(){
	//traverse VCF
	while(_vcfFile.next()){
		//reset counts
		_populations.clear();

		//add data at current line
		for(uint32_t s = 0; s<_samples.numSamples(); ++s){

		}

			genoFrequencies.estimate(&storage[samples.startIndex(p)], samples.numSamplesInPop(p), glfConverter, epsF);
 			_writeEstimatesOnePop(out, genoFrequencies, genoFrequencies.alleleFrequency, &storage[samples.startIndex(p)], samples.numSamplesInPop(p), MLHWEstimator, BHWEstimator, epsF, writeGenoFreq, doBayesian);

	}
};


}; //end namesapce
