/*
 * TContextQuantifer.cpp
 *
 *  Created on: Jun 7, 2020
 *      Author: phaentu
 */

#include "TContextQuantifer.h"

namespace GenomeTasks{

//----------------------------------------------
// TContextQuantifier
//----------------------------------------------
TContextQuantifier::TContextQuantifier(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_parsed(Parameters, Logfile, RandomGenerator){

};

void TContextQuantifier::_handleAlignment(){
	for(auto& b : _alignment){
		if(b.context != BAM::cNN){
			_contextCounts.add(b.recalibratedQualityAsPhredInt.get(), b.context.get());
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

	for(BAM::BaseContext c = BAM::BaseContext::min(); c < BAM::BaseContext::max(); ++c){
		contextLabels.push_back((std::string) c);
	}

	_contextCounts.writeAsMatrix(outputFileName, "quality", contextLabels);
};

}; //end namespace
