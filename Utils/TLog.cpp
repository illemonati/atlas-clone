#include "TLog.h"

TLog::TLog(){
    isFile=false;
    isVerbose=true;
    printWarnings=true;
    numIndent=0;
    lastLineStartInFile=0;
    fillIndentString();
    numberingLevel = -1;
    gettimeofday(&start, NULL);
}

void TLog::close(){
    newLine();
    if(isFile) file.close();
    isFile=false;
}

TLog::~TLog(){ close();}

void TLog::printTimeSinceStartOfProgam(){
    gettimeofday(&end, NULL);
    float runtime=(end.tv_sec  - start.tv_sec);
    if(isFile) file << std::endl << ">>>>>>> Time since start of program: " << runtime << "s" << std::endl;
    if(isVerbose) std::cout << std::endl << ">>>>>>> Time since start of program: " << runtime << "s" << std::endl;
}

void TLog::openFile(std::string Filename){
    list("Writing log to '" + Filename + "'");
    filename=Filename;
    file.open(filename.c_str());
    if(!file) throw "Unable to open logfile '"+ filename +"'!";
    isFile=true;
    lastLineStartInFile=file.tellp();
}

void TLog::setVerbose(bool Verbose){ isVerbose=Verbose;}

bool TLog::verbose(){return isVerbose;}

void TLog::suppressWarings(){printWarnings=false;}

void TLog::showWarings(){printWarnings=true;}

void TLog::newLine(){
    if(isFile){
        file << std::endl;
        lastLineStartInFile=file.tellp();
    }
    std::cout << std::endl;
}

std::string TLog::getFilename(){
    if(isFile) return filename;
    else return "";
}

void TLog::fillIndentString(){
    indetOnlyTabs="";
    for(int i=0; i<numIndent; ++i) indetOnlyTabs += "   ";
    indent = indetOnlyTabs + "- ";
}

void TLog::addIndent(int n){
    numIndent+=n;
    fillIndentString();
}

void TLog::removeIndent(int n){
    numIndent-=n;
    if(numIndent<0) numIndent=0;
    fillIndentString();
}

void TLog::clearIndent(){
    numIndent = 0;
    fillIndentString();
}

void TLog::endIndent(){
    removeIndent();
}

void TLog::addNumberingLevel(){
    ++numberingLevel;
    numberingIndex.push_back(1);
    addIndent();
}

void TLog::removeNumberingLevel(){
    if(numberingLevel>=0){
        numberingIndex.erase(numberingIndex.end()-1);
        --numberingLevel;
        removeIndent();
    }
}

void TLog::endNumbering(){
    removeNumberingLevel();
}

void TLog::done(){
    if(isFile){
        file << " done!" << std::endl;
        lastLineStartInFile=file.tellp();
    }
    if(isVerbose) std::cout << " done!" << std::endl;
}
