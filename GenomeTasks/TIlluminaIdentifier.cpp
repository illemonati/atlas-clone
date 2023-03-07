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

TIlluminaIdentifier::TIlluminaIdentifier():TGenome_basic(){
    // initialize TGenome_basic stuff
	// fill map with platform unit (key) and read group name
for (size_t i = 0; i < _bamFile.readGroups().size(); i++){
        if(_bamFile.readGroups()[i].sequencingTechnology_PL == "Illumina" || _bamFile.readGroups()[i].sequencingTechnology_PL == "ILLUMINA"){
            rgPU_rgID.insert(std::pair<std::string, std::string>(_bamFile.readGroups()[i].platformUnit_PU, _bamFile.readGroups()[i].name_ID));
        } else {
            logfile().list("Readgroup '" + _bamFile.readGroups()[i].name_ID + "' is not Illumina-sequenced and will therefore not be used.");
        }
    }    
    _out.open(_outputName + "_IlluminaReadGroupsCorrected.bam", _bamFile);
};

void TIlluminaIdentifier::_handleAlignment(){
//TODO: put in checks in case alignment names are not illumina
    //to get read group platform unit from read name, split before the 4th colon
    size_t pos = coretools::str::findNthInstanceInString(4, _bamFile.curName(), ':');
    std::string rgPUinName = _bamFile.curName().substr(0, pos);

    if(rgPUinName != _bamFile.readGroups().getReadGroup(_bamFile.curReadGroupID()).platformUnit_PU){
        size_t newId = _bamFile.readGroups().getId(rgPU_rgID[rgPUinName]);
        if(newId == BAM::TReadGroups::noReadGroupId) UERROR("Illumina read group name of read '" + _bamFile.curName() + "' not found in header!");
        _bamFile.curSetNewReadGroup(newId);
        _counter++;
    }
    _bamFile.writeCurAlignment(_out);
}

void TIlluminaIdentifier::run(){
    while(_bamFile.readNextAlignment()){
        _handleAlignment();
    }
    
    //print logFile here
    if(_counter == 1)
        logfile().startIndent("The read group of " + coretools::str::toString(_counter) + " read was corrected.");
    else
        logfile().startIndent("The read groups of " + coretools::str::toString(_counter) + " reads were corrected.");
    logfile().endIndent();

    _out.close();
}

} // end namespace genometasks
