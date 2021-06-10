/*
 * TGLF.h
 *
 *  Created on: Jul 23, 2017
 *      Author: phaentu
 */

#ifndef TGLF_H_
#define TGLF_H_

#include <cstring>
#include "gzstream.h"
#include <algorithm>
#include <vector>
#include <map>
#include <TPopulationLikelihoodLocus.h>
#include "TParameters.h"
#include "stringFunctions.h"
#include "Types.h"
#include "TFastaBuffer.h"
#include "TGenotypeData.h"
#include "TChromosomes.h"
#include "TTask.h"

namespace GLF{

using namespace GenotypeLikelihoods;

/*
//----------------------------------------------------
// TGlfConverter
// class to converted likelihoods to uint16 and back
//----------------------------------------------------
class TGlfConverter{
private:
	uint16_t _maxVal;
	double _minLikelihood;
	std::vector<double> _likelihoodMap;
	double _logOf10DividedByMinus1000 = log(10.) / -1000.;

public:
	TGlfConverter();

	uint16_t maxValue() const{ return _maxVal; };
	uint16_t toGlfFormat(double scaledLikelihood) const;
	uint16_t log10ToGlfFormat(double log10ScaledLikelihood) const;
	uint16_t phredToGlfFormat(BAM::PhredIntErrorRate phred) const;
	double toScaledLikelihood(uint16_t glfValue) const;
	double operator[](uint16_t glfValue) const { return toScaledLikelihood(glfValue); }
	BAM::PhredIntErrorRate toPhred(uint16_t glfValue) const;
	double toLog10(uint16_t glfValue) const;
	BAM::LogErrorRate toLog(uint16_t glfValue) const;
};

*/

//----------------------------------------------------
//TGlfChromosome
//struct to store info on chromosome
//TODO: derive from TChromosome??
//----------------------------------------------------
class TGlfChromosome{
private:
	void _setPloidy(uint8_t Ploidy){
		if(Ploidy < 1 || Ploidy > 2)
			throw "Currently GLFs only support ploidies 1 and 2 (not " + coretools::str::toString(Ploidy) + ")!";
		ploidy = Ploidy;
		if(ploidy == 1){
			numLikelihoodValues = 4;
		} else {
			numLikelihoodValues = 10;
		}
	};

public:
	//current chromosome
	std::string name;
	uint32_t refId;
	uint32_t length;
	uint8_t ploidy;
	uint8_t numLikelihoodValues; //depends on ploidy
	uint8_t maxNumLikelihoodValues; //maximum possible

	TGlfChromosome(){
		refId = 0;
		length = 0;
		ploidy = 2;
		numLikelihoodValues = 10;
		maxNumLikelihoodValues = 10;
		name = "";
	};

	TGlfChromosome(const std::string& Name, uint32_t Length, uint8_t Ploidy){
		name = Name;
		length = Length;
		refId = 0;
		maxNumLikelihoodValues = 10;
        _setPloidy(Ploidy);
	};

	TGlfChromosome(const TGlfChromosome & other){
		update(other);
	};

	void update(const std::string& Name, const uint16_t RefId, const uint32_t Length, const uint8_t Ploidy){
		name = Name;
		refId = RefId;
		length = Length;
        _setPloidy(Ploidy);
	};

	void update(const TGlfChromosome & other){
		name = other.name;
		refId = other.refId;
		length = other.length;
		ploidy = other.ploidy;
		numLikelihoodValues = other.numLikelihoodValues;
		maxNumLikelihoodValues = other.maxNumLikelihoodValues;
	};

	void clear(){
		name = "";
		refId = 0;
		length = 0;
		ploidy = 2;
	};
};

//----------------------------------------------------
//TGlfHandle
//----------------------------------------------------
class TGlfHandle{
protected:
	std::string _filename;
	gzFile _gzfp;
	bool _isOpen;
	uint32_t _offset;
	std::string _version;
	uint8_t _zero8, _one8;
	uint32_t _zero32;
	std::string _header;
	uint64_t _positionInFile;
	TGlfChromosome _curChr;

public:
	TGlfHandle(){
        _isOpen = false;
        _gzfp = nullptr;
        _offset = 0;
        _version = "GLF2"; //change to next version if older files will not work!
		_zero8 = 0;
        _one8 = 1;
        _zero32 = 0;
        _positionInFile = 0;
	};

	void close(){
		if(_isOpen){
			gzclose(_gzfp);
            _isOpen = false;
		}
	};

	std::string name(){
		return _filename;
	};

	std::string chr() const{
		return _curChr.name;
	};

	uint32_t refId() const{
		return _curChr.refId;
	};

	uint32_t chrLength() const{
		return _curChr.length;
	};

	bool chrIsHaploid() const{
		return _curChr.ploidy == 1;
	};

	uint8_t chrPloidy() const{
		return _curChr.ploidy;
	};

	uint8_t chrNumLikelihoodValues() const{
		return _curChr.numLikelihoodValues;
	};
};

//----------------------------------------------------
//TGlfWriter
//----------------------------------------------------
class TGlfWriter:public TGlfHandle{
private:
	long _oldPos;
	uint8_t _recordType1;
	BAM::HighPrecisionPhredIntErrorRate* _glfValues; //tmp used for writing

	void _init();
	void _writeHeader();

	template <typename T>
	inline void _write(const T* buf, size_t len){
        _positionInFile += gzwrite(_gzfp, buf, len);
	};

public:
	TGlfWriter(){
        _init();
	};
	TGlfWriter(const std::string& Filename){
        _init();
		open(Filename, "");
	};

	~TGlfWriter(){
		close();
		delete[] _glfValues;
	};

	//open & close streams
	void open(const std::string& Filename);
	void open(const std::string& Filename, const std::string& Header);
	void newChromosome(const BAM::TChromosome & chromosome);
	void writeSite(long pos, uint32_t depth, uint8_t RMS_mappingQual, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods);
};


//----------------------------------------------------
//TGlfReader
//----------------------------------------------------
class TGlfReader:public TGlfHandle{
private:
    // file parsing
	bool _reachedEndOfChr;
	uint32_t _HeaderLen;
	uint8_t _tmpInt8;
	int _SNPRecordSize;
	uint8_t _tmpRecordStorage[19];
	int _lenRead;
	bool _eof;

	// about site
	int _recordType;
	uint32_t _position;
	uint16_t _depth;
	int _RMS_mappingQual;
	BAM::HighPrecisionPhredIntErrorRate _genotypeLikelihoodsGLF[10];
	BAM::HighPrecisionPhredIntErrorRate* _genotypeLikelihoodsGLF_missingData;

	// about chromosomes
	std::map< uint32_t, TGlfChromosome > _chromosomesAlreadyParsed;

	// initializing
	void _init();
	void _open();

	// read
    template <typename T> inline bool _read(T* buf, size_t len){
        _lenRead = gzread(_gzfp, buf, len);
        if(!_lenRead) return false;
        _positionInFile += _lenRead;
        return true;
    };
	bool _readChr();
	bool _readRecordType();
	void _readSNPRecord();
	inline void _skipRecord(){
        _read(&_tmpRecordStorage, _SNPRecordSize); //just parse data of one SNP record into garbage
	};

    bool _jumpToEndOfChr();

public:
	TGlfReader(){
        _init();
	};
	TGlfReader(const std::string& Filename){
        _init();
		open(Filename);
	};
	~TGlfReader(){
		close();
		delete[] _genotypeLikelihoodsGLF_missingData;
	};

	//get details
	bool eof() const{ return _eof;};
	TGlfChromosome* pointerToChr(const uint32_t & refId);
	bool fillPointerToChr(uint32_t refId, TGlfChromosome* & chr);
	uint32_t position() const{ return _position; };
	uint16_t depth() const{ return _depth; };
	BAM::HighPrecisionPhredIntErrorRate* pointerToGenotypeLikelihoodsGLF(){ return _genotypeLikelihoodsGLF; };
    void fillGenotypeLikelihoodsGLF(BAM::HighPrecisionPhredIntErrorRate* destination);
    //void fillGenotypeLikelihoods(BAM::ErrorRate* destination);

	//open file and parse header
	void setFilename(const std::string& Filename);
	void open(const std::string& Filename);
	void rewind();

	// parse file
	bool readNext();
	bool jumpToNextChr();
	bool readNextWindow(std::vector<uint16_t*> & genoLikelihoods, const uint32_t refId, const uint32_t start, const uint32_t end);

	//printing
	void printChr();
	void printSite();
	void printToEnd();
};


//------------------------------------------------
// Tasks
//------------------------------------------------
class TTask_printGLF:public coretools::TTask{
public:
	TTask_printGLF(){ _explanation = "Printing a GLF file to screen"; };

	void run(coretools::TParameters & parameters, coretools::TLog* logfile){
		std::string glf = parameters.getParameter<std::string>("glf");
		TGlfReader reader(glf);
		reader.printToEnd();
	};
};

}; //end namespace

#endif /* TGLF_H_ */
