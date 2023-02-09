/*
 * TIlluminaIdentifier.cpp
 *
 *  Created on: Feb 9, 2023
 *      Author: raphael
 */


#include "TIlluminaIdentifier.h"
#include "coretools/Files/TFile.h"
#include "TBamFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Strings/stringFunctions.h"

namespace GenomeTasks{

using coretools::instances::parameters;
using coretools::instances::logfile;

TIlluminaIdentifier::TIlluminaIdentifier():TGenome_filtered(){
    // initialize TGenome_basic stuff
	// fill map with platform unit (key) and read group name
    for (const auto &s: _bamFile.readGroups()){
        rgPU_rgID.insert(std::pair<s.platformUnit_PU, s.nameID>);
    }
    _bamFile.open("test");
};

void TIlluminaIdentifier::_handleAlignment(){
//TODO: put in checks in case alignment names are not illumina
    //to get read group platform unit from read name, split before the 4th colon
    size_t pos = findNthInstanceInString(4, _bamFile.curName(), ':');
    std::string rgPUinName = _bamFile.curName().substr(0, pos-1);

    if(rgPUinName != _bamFile.readGroups().getName(_bamFile.curReadGroupID())){
        size_t newId = _bamFile.readGroups().getId(rgPUinName);
        if(newId == noReadGroupId) UERROR("Illumina read group name of read " + _bamFile.curName() + " not found in header!");
        _bamFile.curSetNewReadGroup(newId);
    }
    _bamFile.writeCurAlignment();
}

void TIlluminaIdentifer::identify(){
    _traverseBAMPassedQC();
    //print logFile here

}

} // end namespace genometasks
