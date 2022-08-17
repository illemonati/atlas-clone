/*
 * TContextQuantifer.cpp
 *
 *  Created on: Jun 7, 2020
 *      Author: phaentu
 */


#include "TContextQuantifer.h"

#include <algorithm>
#include <vector>

#include "GenotypeTypes.h"
#include "PhredProbabilityTypes.h"
#include "TAlignment.h"
#include "TLog.h"
#include "TSequencedBase.h"

namespace GenomeTasks{

using coretools::TParameters;
using coretools::TLog;
using coretools::TRandomGenerator;

//----------------------------------------------
// TContextQuantifier
//----------------------------------------------
TContextQuantifier::TContextQuantifier():TGenome_parsed(){

};

void TContextQuantifier::_handleAlignment(){
	for(auto& b : _alignment){
		if(b.context != genometools::BaseContext::NN){
			_contextCounts.add(b.recalibratedQualityAsPhredInt.get(), coretools::index(b.context));
		}
	}
};

void TContextQuantifier::quantifyContexts(){
	using coretools::instances::logfile;
	_contextCounts.clear();

	//traverse BAM and add to counts
	_traverseBAMPassedQC();

	//write counts
	std::string outputFileName = _outputName + "_contextInformation.txt.gz";
	logfile().list("Writing context information to file '" + outputFileName + "'.");

	std::vector<std::string> contextLabels;

	for(genometools::BaseContext c = genometools::BaseContext::min; c < genometools::BaseContext::max; ++c){
		contextLabels.push_back(genometools::toString(c));
	}

	_contextCounts.writeAsMatrix(outputFileName, "quality", contextLabels);
};

}; //end namespace
