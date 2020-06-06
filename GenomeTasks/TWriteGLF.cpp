/*
 * TGLFWriter.cpp
 *
 *  Created on: Jun 6, 2020
 *      Author: phaentu
 */


#include <TWriteGLF.h>

namespace GenomeTasks{

//-------------------------------------------
// TGLFWriter
//-------------------------------------------
TWriteGLF::TWriteGLF(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_windows(Params, Logfile, RandomGenerator){
	if(Params.parameterExists("printAll")){
		_printAll = true;
		_logfile->list("Will write all sites, even those without data. (parameter 'printAll')");
	} else {
		_printAll = false;
		_logfile->list("Will only write sites with data. (use 'printAll' to write all sites)");
	}
};


void TWriteGLF::_handleWindow(){
	if(_chrChangedWindow){
		_writer.newChromosome(_chromosomes.curChromosome());
	}

	//TODO: calculate root mean squared mapping qualities for sites (now just passing 0). Would be helpful in VCFs as well
	_logfile->listFlushTime("Adding window to GLF file ...");
	uint32_t pos = 0;
	for(auto& s : _window){
		if(s.hasData || _printAll){
			_genotypeLikelihoodCalculator.calculateGenotypeLikelihoods(s.bases, _genoLik);
			_writer.writeSite(_window.posInRef(pos), s.depth(), 0, _genoLik);
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
