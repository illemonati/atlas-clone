/*
 * TContextQuantifer.cpp
 *
 *  Created on: Jun 7, 2020
 *      Author: phaentu
 */

#include "TContextQuantifer.h"
#include "GenotypeTypes.h"

namespace GenomeTasks{

//----------------------------------------------
// TContextQuantifier
//----------------------------------------------
TContextQuantifier::TContextQuantifier(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_parsed(Parameters, Logfile, RandomGenerator){

};

void TContextQuantifier::_handleAlignment(){
	for(auto& b : _alignment){
		if(b.context != genometools::BaseContext::NN){
			_contextCounts.add(b.recalibratedQualityAsPhredInt.get(), genometools::index(b.context));
		}
	}
};

void TContextQuantifier::quantifyContexts(){
	_contextCounts.clear();

	//traverse BAM and add to counts
	_traverseBAMPassedQC();

	//write counts
	std::string outputFileName = _outputName + "_contextInformation.txt.gz";
	_logfile->list("Writing context information to file '" + outputFileName + "'.");

	std::vector<std::string> contextLabels;

	for(genometools::BaseContext c = genometools::BaseContext::min; c < genometools::BaseContext::max; ++c){
		contextLabels.push_back(genometools::toString(c));
	}

	_contextCounts.writeAsMatrix(outputFileName, "quality", contextLabels);
};

}; //end namespace
