/*
 * TIlluminaIdentifier.cpp
 *
 *  Created on: Feb 9, 2023
 *      Author: raphael
 */


#include "TIlluminaIdentifier.h"
#include "coretools/Files/TFile.h"
#include "TBamFile.h"
#include <string>
#include "coretools/Main/TLog.h"
#include "coretools/Strings/stringFunctions.h"

namespace GenomeTasks{

using coretools::instances::logfile;

TIlluminaIdentifier::TIlluminaIdentifier(): _out(_genome.outputName() + "_IlluminaReadGroupsCorrected.bam", _genome.bamFile()){
    // initialize TGenome_basic stuff
	// fill map with platform unit (key) and read group name
for (size_t i = 0; i < _genome.bamFile().readGroups().size(); i++){
        if(_genome.bamFile().readGroups()[i].sequencingTechnology_PL == "Illumina" || _genome.bamFile().readGroups()[i].sequencingTechnology_PL == "ILLUMINA"){
            rgPU_rgID.insert(std::pair<std::string, std::string>(_genome.bamFile().readGroups()[i].platformUnit_PU, _genome.bamFile().readGroups()[i].name_ID));
        } else {
            logfile().list("Readgroup '" + _genome.bamFile().readGroups()[i].name_ID + "' is not Illumina-sequenced and will therefore not be used.");
        }
    }    
};

void TIlluminaIdentifier::_handleAlignment(){
//TODO: put in checks in case alignment names are not illumina
    //to get read group platform unit from read name, split before the 4th colon
    size_t pos = coretools::str::findNthInstanceInString(4, _genome.bamFile().curName(), ':');
    std::string rgPUinName = _genome.bamFile().curName().substr(0, pos);

    if(rgPUinName != _genome.bamFile().readGroups().getReadGroup(_genome.bamFile().curReadGroupID()).platformUnit_PU){
        size_t newId = _genome.bamFile().readGroups().getId(rgPU_rgID[rgPUinName]);
        if(newId == BAM::TReadGroups::noReadGroupId) UERROR("Illumina read group name of read '" + _genome.bamFile().curName() + "' not found in header!");
        _genome.bamFile().curSetNewReadGroup(newId);
        _counter++;
    }
    _genome.bamFile().writeCurAlignment(_out);
}

void TIlluminaIdentifier::run(){
    while(_genome.bamFile().readNextAlignment()){
        _handleAlignment();
    }
    
    //print logFile here
    if(_counter == 1)
        logfile().startIndent("The read group of " + coretools::str::toString(_counter) + " read was corrected.");
    else
        logfile().startIndent("The read groups of " + coretools::str::toString(_counter) + " reads were corrected.");
    logfile().endIndent();
}

} // end namespace genometasks
