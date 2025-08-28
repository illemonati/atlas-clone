/*
 * TGLFWriter.cpp
 *
 *  Created on: Jun 6, 2020
 *      Author: phaentu
 */


#include "TWriteGLFOld.h"

#include <stdint.h>

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
TWriteGLFOld::TWriteGLFOld(){
	if(parameters().exists("printAll")){
		_printAll = true;
		logfile().list("Will write all sites, even those without data. (parameter 'printAll')");
	} else {
		_printAll = false;
		logfile().list("Will only write sites with data. (use 'printAll' to write all sites)");
	}
};

void TWriteGLFOld::_startChromosome(const genometools::TChromosome &Chr) {
	_writer.newChromosome(Chr);
	_curIsHapo = Chr.isHaploid();
}

void TWriteGLFOld::_handleWindow(GenotypeLikelihoods::TWindow &Window) {
	//TODO: calculate root mean squared mapping qualities for sites (now just passing 0). Would be helpful in VCFs as well
	logfile().listFlushTime("Adding window to GLF file ...");
	uint32_t pos = 0;
	for(auto& s : Window){
		if(!s.empty() || _printAll){
			if (_curIsHapo) {
				const auto baseLik = _genome.errorModels().calculateBaseLikelihoods(s);
				_writer.writeSite(Window.positionOnChr(pos), s.depth(), s.rmsMappingQual(), baseLik);
			} else {
				const auto genoLik = _genome.errorModels().calculateGenotypeLikelihoods(s);
				_writer.writeSite(Window.positionOnChr(pos), s.depth(), s.rmsMappingQual(), genoLik);
			}
		}
		++pos;
	}
	logfile().doneTime();
};

void TWriteGLFOld::run(){
	//open GLF file
	const auto outputFileName = _genome.outputName() + ".glf.gz";
	logfile().list("Will write genotype likelihoods to GLF file '" + outputFileName + "'.");

	

	_writer.open(outputFileName, _chromosomes(this->_genome));

	//traverse BAM
	_traverseBAMWindows();

	//close writer
	_writer.close();
};


}; //end namespace
