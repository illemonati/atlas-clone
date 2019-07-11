#include "TReadList.h"

void TReadList::_write(const std::string &Name, const std::string &Error, const bool &isSecondMate){
    if(_writeToFile){
        if(isSecondMate){
            ignoredReads << "Read " << Name << " (second mate): " << Error << "\n";
        } else {
            ignoredReads << "Read " << Name << " (first mate): " << Error << "\n";
        }
    }
}

TReadList::TReadList(){
    _writeToFile = false;
}

TReadList::TReadList(const bool &writeToFile, std::string &filename){
    _writeToFile = writeToFile;
}

void TReadList::setWriteReadListToFileToTrue(TLog *logfile, std::string &Filename){
    _writeToFile = true;
    filename = Filename + "_ignoredReads.txt.gz";
    logfile->list("Writing sequencing depth estimates to '" + filename + "'");
    ignoredReads.open(filename.c_str());
    if(!ignoredReads) throw "Failed to open output file '" + filename + "'!";
}

void TReadList::addToReadList(TAlignment &alignment, const std::string &errorMessage){
    //TODO: should check if read already exists in blackfile (could be case in paired-end data) -> remove
    readList.emplace(alignment.name(), 1);
    _write(alignment.name(), errorMessage, alignment.isSecondMate);
}

void TReadList::addToReadList(std::string &alignmentName, const std::string &errorMessage){
    //TODO: should check if read already exists in blackfile (could be case in paired-end data) -> remove
    readList.emplace(alignmentName, 1);
    if(_writeToFile){
        ignoredReads << "Read " << alignmentName << ": " << errorMessage << "\n";
    }
}

void TReadList::removeFromReadList(TAlignment &alignment, const std::string &errorMessage){
    readList.erase(alignment.name());
    _write(alignment.name(), errorMessage, alignment.isSecondMate);
}

bool TReadList::isInReadList(const std::string &alignmentName){
    if(readList.count(alignmentName) > 0)
        return true;
    return false;
}
