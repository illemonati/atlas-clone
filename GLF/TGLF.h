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
#include <array>
#include <TPopulationLikelihoodLocus.h>
#include "TParameters.h"
#include "stringFunctions.h"
#include "GenotypeTypes.h"
#include "TFastaBuffer.h"
#include "TGenotypeData.h"
#include "TChromosomes.h"
#include "TTask.h"

namespace GLF{


typedef std::array< genometools::HighPrecisionPhredIntProbability, 10 > GLFLikelihoods;

//----------------------------------------------------
//TGlfChromosome
//struct to store info on chromosome
//TODO: derive from TChromosome??
//----------------------------------------------------
class TGlfChromosome{
private:
	std::string _name;
	uint32_t _refId;
	uint32_t _length;
	bool _isHaploid;
	uint8_t _numLikelihoodValues; //depends on ploidy

	void _setPloidy(const uint8_t & Ploidy);

public:

	TGlfChromosome();
	TGlfChromosome(const std::string & Name, uint32_t Length, const uint8_t & Ploidy);

	std::string name() const { return _name; };
	uint32_t refId() const { return _refId; };
	uint32_t length() const { return _length; };
	bool isHaploid() const { return _isHaploid; };
	uint8_t ploidy() const { return 2 - _isHaploid; };
	uint8_t numLikelihoodValues() const { return _numLikelihoodValues; };

	void update(const std::string & Name, uint16_t RefId, uint32_t Length, const uint8_t & Ploidy);

	void clear();
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
		return _curChr.name();
	};

	uint32_t refId() const{
		return _curChr.refId();
	};

	uint32_t chrLength() const{
		return _curChr.length();
	};

	bool chrIsHaploid() const{
		return _curChr.isHaploid();
	};

	uint8_t chrNumLikelihoodValues() const{
		return _curChr.numLikelihoodValues();
	};
};

//----------------------------------------------------
//TGlfWriter
//----------------------------------------------------
class TGlfWriter:public TGlfHandle{
private:
	long _oldPos;
	uint8_t _recordType1;
	GLFLikelihoods _glfValues; //tmp used for writing

	void _init();
	void _writeHeader();

	template <typename T> void _write(T var){
        _positionInFile += gzwrite(_gzfp, &var, sizeof(T));
	};
	void _write(const void * buf, size_t len){
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
	GLFLikelihoods _genotypeLikelihoodsGLF;
	GLFLikelihoods _genotypeLikelihoodsGLF_missingData;

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
	};

	//get details
	bool eof() const{ return _eof;};
	TGlfChromosome* pointerToChr(uint32_t refId);
	bool fillPointerToChr(uint32_t refId, TGlfChromosome* & chr);
	uint32_t position() const{ return _position; };
	uint16_t depth() const{ return _depth; };
	const GLFLikelihoods& genotypeLikelihoodsGLF() const { return _genotypeLikelihoodsGLF; };

	//open file and parse header
	void setFilename(const std::string& Filename);
	void open(const std::string& Filename);
	void rewind();

	// parse file
	bool readNext();
	bool jumpToNextChr();
	bool readNextWindow(std::vector<GLFLikelihoods> & genoLikelihoods, uint32_t refId, uint32_t start, uint32_t end);

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

	void run(){
		using coretools::instances::parameters;
		std::string glf = parameters().getParameter<std::string>("glf");
		TGlfReader reader(glf);
		reader.printToEnd();
	};
};

}; //end namespace

#endif /* TGLF_H_ */
