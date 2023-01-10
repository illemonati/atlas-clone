/*
 * TContextQuantifer.cpp
 *
 *  Created on: Jun 7, 2020
 *      Author: phaentu
 */


#include "TContextQuantifer.h"

#include <algorithm>
#include <vector>

#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"
#include "TAlignment.h"
#include "coretools/Main/TLog.h"
#include "TSequencedBase.h"

namespace GenomeTasks{

//----------------------------------------------
// TContextQuantifier
//----------------------------------------------
TContextQuantifier::TContextQuantifier():TGenome_parsed(){

};

void TContextQuantifier::_handleAlignment(){
	for(auto& b : _alignment){
		genometools::BaseContext context = b.context();
		if(context != genometools::BaseContext::NN){					
			_contextCounts.add(b.readGroupID, coretools::index(context));
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
