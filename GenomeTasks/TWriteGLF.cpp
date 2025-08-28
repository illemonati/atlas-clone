/*
 * TGLFWriter.cpp
 *
 *  Created on: Jun 6, 2020
 *      Author: phaentu
 */


#include "TWriteGLF.h"

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
TWriteGLF::TWriteGLF(){
	if(parameters().exists("printAll")){
		_printAll = true;
		_siteTraverser.setMinDepth(0);
		logfile().list("Will write all sites, even those without data. (parameter 'printAll')");
	} else {
		_printAll = false;
		_siteTraverser.setMinDepth(1);
		logfile().list("Will only write sites with data. (use 'printAll' to write all sites)");
	}
};

void TWriteGLF::_traverseSites() {
	for (; !_siteTraverser.endOfChrs(); _siteTraverser.nextChr()) {
		logfile().listFlushTime("Adding chr '", _siteTraverser.curChr().name(), "' of length ", _siteTraverser.curChr().length(), " to GLF file ...");
		_writer.newChromosome(_siteTraverser.curChr());
		for (; !_siteTraverser.endOfCurChr(); _siteTraverser.nextSite()) {
			const auto &site = _siteTraverser.site();
			if (!site.empty() || _printAll) {
				if (_siteTraverser.curChr().isHaploid()) {
					const auto baseLik = _siteTraverser.errorModels().calculateBaseLikelihoods(site);
					_writer.writeSite(_siteTraverser.position().position(), site.depth(), site.rmsMappingQual(),
									  baseLik);
				} else {
					const auto genoLik = _siteTraverser.errorModels().calculateGenotypeLikelihoods(site);
					_writer.writeSite(_siteTraverser.position().position(), site.depth(), site.rmsMappingQual(),
									  genoLik);
				}
			}
		}
		logfile().doneTime();
	}
};

void TWriteGLF::run(){
	//open GLF file
	const auto outputFileName = _siteTraverser.outputName() + ".glf.gz";
	logfile().list("Will write genotype likelihoods to GLF file '" + outputFileName + "'.");

	_writer.open(outputFileName, _siteTraverser.chromosomes());

	//traverse BAM
	_traverseSites();

	//close writer
	_writer.close();
};


}; //end namespace
