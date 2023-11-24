/*
 * TGLFWriter.cpp
 *
 *  Created on: Jun 6, 2020
 *      Author: phaentu
 */


#include "TWriteGLF.h"

#include <stdint.h>
#include <vector>

#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "TSite.h"
#include "TWindow.h"

namespace GenomeTasks{
using coretools::instances::logfile;
using coretools::instances::parameters;

//-------------------------------------------
// TGLFWriter
//-------------------------------------------
TWriteGLF::TWriteGLF(){
	if(parameters().exists("printAll")){
		_printAll = true;
		logfile().list("Will write all sites, even those without data. (parameter 'printAll')");
	} else {
		_printAll = false;
		logfile().list("Will only write sites with data. (use 'printAll' to write all sites)");
	}
};

void TWriteGLF::_handleWindow(GenotypeLikelihoods::TWindow& window){
	if(_chrChangedWindow){
		_writer.newChromosome(*_curChromosome);
	}

	//TODO: calculate root mean squared mapping qualities for sites (now just passing 0). Would be helpful in VCFs as well
	logfile().listFlushTime("Adding window to GLF file ...");
	uint32_t pos = 0;
	for(auto& s : window){
		if(!s.empty() || _printAll){
			const auto genoLik = _parser.errorModels().calculateGenotypeLikelihoods(s);
			_writer.writeSite(window.positionOnChr(pos), s.depth(), 0, genoLik);
		}
		++pos;
	}
	logfile().doneTime();
};

void TWriteGLF::run(){
	//open GLF file
	const auto outputFileName = _genome.outputName() + ".glf.gz";
	logfile().list("Will write genotype likelihoods to GLF file '" + outputFileName + "'.");
	_writer.open(outputFileName);

	//traverse BAM
	_traverseBAMWindows();

	//close writer
	_writer.close();
};


}; //end namespace
