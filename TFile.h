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
    TFile();
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

    void _setName(const std::string & Name);

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
    }

    virtual void _openFile();
    virtual void _closeFile();

public:
    TOutputFile();
    virtual ~TOutputFile();

    void open(const std::string Name);
    void close();

	//writing header
    void writeHeader(std::vector<std::string> header);

    void noHeader(int NumCols);

	//writing data
    void setPrecision(int precision);

	template<typename T>
	void writeLine(std::vector<T> & data){
		if(_colsWritten != 0)
			throw "Can not write line: Last line was not finished in file '" + _name + "'!";
		if(data.size() != _numCols)
			throw "Can not write line: wrong number of columns (" + toString(data.size()) + " instead of " + toString(_numCols) + ")!";
		_write(data);
		endLine();
    }

	template<typename T>
	void write(T data){
		if(_colsWritten == _numCols)
			throw "Extra value on line " + toString(_lineNum) + " in file '" + _name + "'!";

		if(_colsWritten > 0) *_filePointer << "\t";
		*_filePointer << data;
		++_colsWritten;
    }

	template<typename T>
	void write(const std::vector<T> & data){
		if(_colsWritten + data.size() > _numCols)
			throw "Extra values on line " + toString(_lineNum) + " in file '" + _name + "'!";

		_write(data);
    }

	template<typename T>
	TOutputFile& operator<<(const T data){
		write(data);
		return *this;
    }

	//new line
    void endLine();

	//overloading std::endl
	typedef TOutputFile& (*TOutputFileManipulator)(TOutputFile&);

	// take in a function with the custom signature
    TOutputFile& operator<<(TOutputFile::TOutputFileManipulator manip);

    static TOutputFile& endl(TOutputFile& of);

	// this is the type of std::cout
	typedef std::basic_ostream<char, std::char_traits<char> > CoutType;

	// this is the function signature of std::endl
	typedef CoutType& (*StandardEndLine)(CoutType&);

	// define an operator<< to take in std::endl
    TOutputFile& operator<<(StandardEndLine manip);
};


//-------------------------------------------------------
// TOutputFilePlain
//-------------------------------------------------------
class TOutputFilePlain:public TOutputFile{
private:
	std::ofstream file;

    void _openFile();
    void _closeFile();

public:
    TOutputFilePlain();
    TOutputFilePlain(const std::string Name);

};

//-------------------------------------------------------
// TOutputFileZipped
//-------------------------------------------------------
class TOutputFileZipped:public TOutputFile{
private:
	gz::ogzstream zippedFile;

    void _openFile();
    void _closeFile();

public:
    TOutputFileZipped();
    TOutputFileZipped(const std::string Name);

};


#endif /* TFILE_H_ */
