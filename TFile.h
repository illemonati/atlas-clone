/*
 * TFile.h
 *
 *  Created on: Jan 8, 2019
 *      Author: phaentu
 */

#ifndef TFILE_H_
#define TFILE_H_

#include "stringFunctions.h"
#include <iomanip>
#include "gzstream.h"

class TFile{
protected:
	std::string _name;
	bool _isOpen;
	int _lineNum;

public:
	TFile(){
		_isOpen = false;
		_lineNum = 0;
	};
};

class TInputFile:public TFile{
protected:

public:

};


//-------------------------------------------------------
// TOutputFile
//-------------------------------------------------------
class TOutputFile:public TFile{
protected:
	std::ostream* _filePointer;
	bool _headerWritten;
	unsigned int _numCols;
	unsigned int _colsWritten;

	void _setName(const std::string & Name){
		_name = Name;
		if(_name.empty())
			throw "Please provide a valid file name!";
	};

	template<class T>
	void _write(const std::vector<T> & data){
		if(!_isOpen)
			throw "Can not write to file: file not open!";

		//write first of line
		typename std::vector<T>::const_iterator it=data.begin();
		if(_colsWritten == 0){
			*_filePointer << *it;
			it++;
		}

		//write others
		for(; it!=data.end(); it++){
			*_filePointer << "\t" << *it;
		}

		_colsWritten += data.size();
	};

	virtual void _openFile(){
		throw "void _openFile() not implemented for base class TOutputFile!";
	};

	virtual void _closeFile(){
		throw "void _closeFile() not implemented for base class TOutputFile!";
	};

public:
	TOutputFile():TFile(){
		_headerWritten = false;
		_filePointer = NULL;
		_numCols = 0;
		_colsWritten = 0;
	};

	void open(const std::string Name){
		_setName(Name);
		_openFile();
		_isOpen = true;
		_headerWritten = false;
		_numCols = 0;
		_colsWritten = 0;
		_lineNum = 1;
	};

	void close(){
		_closeFile();
	};

	virtual ~TOutputFile(){};

	//writing header
	void writeHeader(std::vector<std::string> header){
		if(_headerWritten)
			throw "Header already written to file '" + _name + "'!";

		_numCols = header.size();
		_write(header);
		endLine();
		_headerWritten = true;
	};

	void noHeader(int NumCols){
		_headerWritten = true;
		_numCols = NumCols;
	};

	//writing data
	void setPrecision(int precision){
		*_filePointer << std::setprecision(precision);
	};

	template<typename T>
	void writeLine(std::vector<T> & data){
		if(_colsWritten != 0)
			throw "Can not write line: Last line was not finished in file '" + _name + "'!";
		if(data.size() != _numCols)
			throw "Can not write line: wrong number of columns (" + toString(data.size()) + " instead of " + toString(_numCols) + ")!";
		_write(data);
		endLine();
	};

	template<typename T>
	void write(T data){
		std::cout << "numCols " << _numCols << std::endl;
		std::cout << "_colsWritten " << _colsWritten << std::endl;
		if(_colsWritten == _numCols)
			throw "Extra value on line " + toString(_lineNum) + " in file '" + _name + "'!";

		if(_colsWritten > 0) *_filePointer << "\t";
		*_filePointer << data;
		++_colsWritten;
	};

	template<typename T>
	void write(const std::vector<T> & data){
		if(_colsWritten + data.size() > _numCols)
			throw "Extra values on line " + toString(_lineNum) + " in file '" + _name + "'!";

		_write(data);
	};

	template<typename T>
	TOutputFile& operator<<(const T data){
		write(data);
		return *this;
	};

	//new line
	void endLine(){
		if(_colsWritten != _numCols)
			throw "Can not end line in file '" + _name + "': missing values on line " + toString(_lineNum) + "!";
		*_filePointer << "\n";
		_colsWritten = 0;
		++_lineNum;
	};

	//overloading std::endl
	typedef TOutputFile& (*TOutputFileManipulator)(TOutputFile&);

	// take in a function with the custom signature
	TOutputFile& operator<<(TOutputFile::TOutputFileManipulator manip){
		return manip(*this);
	};

	static TOutputFile& endl(TOutputFile& of){
		of.endLine();
		return of;
	};

	// this is the type of std::cout
	typedef std::basic_ostream<char, std::char_traits<char> > CoutType;

	// this is the function signature of std::endl
	typedef CoutType& (*StandardEndLine)(CoutType&);

	// define an operator<< to take in std::endl
	TOutputFile& operator<<(StandardEndLine manip){
	    this->endLine();
	    return *this;
	};
};


//-------------------------------------------------------
// TOutputFilePlain
//-------------------------------------------------------
class TOutputFilePlain:public TOutputFile{
private:
	std::ofstream file;

	void _openFile(){
		file.open(_name);
		if(!file)
			throw "Failed to open file '" + _name + "' for writing!";
	};

	void _closeFile(){
		if(_isOpen){
			file.close();
			_isOpen = false;
		}
	};

public:
	TOutputFilePlain():TOutputFile(){
		_filePointer = &file;
	};

	TOutputFilePlain(const std::string Name):TOutputFile(){
		_filePointer = &file;
		open(Name);
	};
};

//-------------------------------------------------------
// TOutputFileZipped
//-------------------------------------------------------
class TOutputFileZipped:public TOutputFile{
private:
	gz::ogzstream zippedFile;

	void _openFile(){
		zippedFile.open(_name.c_str());
		if(!zippedFile)
			throw "Failed to open file '" + _name + "' for writing!";
	};

	void _closeFile(){
		if(_isOpen){
			zippedFile.close();
			_isOpen = false;
		}
	};

public:
	TOutputFileZipped():TOutputFile(){
		_filePointer = &zippedFile;
	};

	TOutputFileZipped(const std::string Name):TOutputFile(){
		_filePointer = &zippedFile;
		open(Name);
	};
};


#endif /* TFILE_H_ */
