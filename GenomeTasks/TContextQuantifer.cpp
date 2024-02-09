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
TContextQuantifier::TContextQuantifier() : TBamReadTraverser<ReadType::Parsed>(){

};

void TContextQuantifier::_handleAlignment(BAM::TAlignment& alignment){
	for(auto& b : alignment){
		genometools::BaseContext context = b.context();
		if(context != genometools::BaseContext::NN){					
			_contextCounts.add(b.readGroupID, coretools::index(context));
		}
	}
};

void TContextQuantifier::run(){
	using coretools::instances::logfile;
	_contextCounts.clear();

	//traverse BAM and add to counts
	_traverseBAMPassedQC();

	//write counts
	std::string outputFileName = _genome.outputName() + "_contextInformation.txt.gz";
	logfile().list("Writing context information to file '" + outputFileName + "'.");

	std::vector<std::string> contextLabels;

	for(genometools::BaseContext c = genometools::BaseContext::min; c < genometools::BaseContext::max; ++c){
		contextLabels.push_back(genometools::toString(c));
	}

	_contextCounts.writeAsMatrix(outputFileName, "quality", contextLabels);
};

}; //end namespace
