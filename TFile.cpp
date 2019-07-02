#include "TFile.h"

//TFile
TFile::TFile(){
    _isOpen = false;
    _lineNum = 0;
}


//TOutputFile
void TOutputFile::_setName(const std::string &Name){
    _name = Name;
    if(_name.empty())
        throw "Please provide a valid file name!";
}

void TOutputFile::_openFile(){
    throw "void _openFile() not implemented for base class TOutputFile!";
}

void TOutputFile::_closeFile(){
    throw "void _closeFile() not implemented for base class TOutputFile!";
}

TOutputFile::TOutputFile():TFile(){
    _headerWritten = false;
    _filePointer = NULL;
    _numCols = 0;
    _colsWritten = 0;
}

void TOutputFile::open(const std::string Name){
    _setName(Name);
    _openFile();
    _isOpen = true;
    _headerWritten = false;
    _numCols = 0;
    _colsWritten = 0;
    _lineNum = 1;
}

void TOutputFile::close(){
    _closeFile();
}

TOutputFile::~TOutputFile(){}

void TOutputFile::writeHeader(std::vector<std::string> header){
    if(_headerWritten)
        throw "Header already written to file '" + _name + "'!";

    _numCols = header.size();
    _write(header);
    endLine();
    _headerWritten = true;
}

void TOutputFile::noHeader(int NumCols){
    _headerWritten = true;
    _numCols = NumCols;
}

void TOutputFile::setPrecision(int precision){
    *_filePointer << std::setprecision(precision);
}

void TOutputFile::endLine(){
    if(_colsWritten != _numCols)
        throw "Can not end line in file '" + _name + "': missing values on line " + toString(_lineNum) + "!";
    *_filePointer << "\n";
    _colsWritten = 0;
    ++_lineNum;
}

TOutputFile &TOutputFile::operator<<(TOutputFile::TOutputFileManipulator manip){
    return manip(*this);
}

TOutputFile &TOutputFile::endl(TOutputFile &of){
    of.endLine();
    return of;
}

TOutputFile &TOutputFile::operator<<(TOutputFile::StandardEndLine manip){
    this->endLine();
    return *this;
}



//TOutputFilePlain
void TOutputFilePlain::_openFile(){
    file.open(_name);
    if(!file)
        throw "Failed to open file '" + _name + "' for writing!";
}

void TOutputFilePlain::_closeFile(){
    if(_isOpen){
        file.close();
        _isOpen = false;
    }
}

TOutputFilePlain::TOutputFilePlain():TOutputFile(){
    _filePointer = &file;
}

TOutputFilePlain::TOutputFilePlain(const std::string Name):TOutputFile(){
    _filePointer = &file;
    open(Name);
}


//TOutputFileZipped
void TOutputFileZipped::_openFile(){
    zippedFile.open(_name.c_str());
    if(!zippedFile)
        throw "Failed to open file '" + _name + "' for writing!";
}

void TOutputFileZipped::_closeFile(){
    if(_isOpen){
        zippedFile.close();
        _isOpen = false;
    }
}

TOutputFileZipped::TOutputFileZipped():TOutputFile(){
    _filePointer = &zippedFile;
}

TOutputFileZipped::TOutputFileZipped(const std::string Name):TOutputFile(){
    _filePointer = &zippedFile;
    open(Name);
}
