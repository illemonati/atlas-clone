/*
 * TGLFWriter.cpp
 *
 *  Created on: Jun 6, 2020
 *      Author: phaentu
 */


#include "TWriteGLF.h"

namespace GenomeTasks{

//-------------------------------------------
// TGLFWriter
//-------------------------------------------
TWriteGLF::TWriteGLF(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_windows(Parameters, Logfile, RandomGenerator){
	if(Parameters.parameterExists("printAll")){
		_printAll = true;
		_logfile->list("Will write all sites, even those without data. (parameter 'printAll')");
	} else {
		_printAll = false;
		_logfile->list("Will only write sites with data. (use 'printAll' to write all sites)");
	}
};

void TWriteGLF::_handleWindow(){
	if(_chrChangedWindow){
		_writer.newChromosome(*_curChromosome);
	}

	//TODO: calculate root mean squared mapping qualities for sites (now just passing 0). Would be helpful in VCFs as well
	_logfile->listFlushTime("Adding window to GLF file ...");
	uint32_t pos = 0;
	for(auto& s : _window){
		if(!s.empty() || _printAll){
			_genotypeLikelihoodCalculator.calculateGenotypeLikelihoods(s, _genoLik);
			_writer.writeSite(_window.positionOnChr(pos), s.depth(), 0, _genoLik);
		}
		++pos;
	}
	_logfile->doneTime();
};

void TWriteGLF::writeGLF(){
	//open GLF file
	std::string outputFileName = _outputName + ".glf.gz";
	_logfile->list("Will write genotype likelihoods to GLF file '" + outputFileName + "'");
	_writer.open(outputFileName);

	//traverse BAM
	_traverseBAMWindows();

	//close writer
	_writer.close();
};


}; //end namespace
