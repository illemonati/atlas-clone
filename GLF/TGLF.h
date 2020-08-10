/*
 * TGLF.h
 *
 *  Created on: Jul 23, 2017
 *      Author: phaentu
 */

#ifndef TGLF_H_
#define TGLF_H_

#include <stdio.h>
#include <zlib.h>
#include "gzstream.h"
#include <algorithm>
#include <string.h>
#include <TPopulationLikelihoodLocus.h>
#include <vector>
#include "TParameters.h"
#include "stringFunctions.h"
#include "TGenotypeMap.h"
#include "TRandomGenerator.h"
#include "TFastaBuffer.h"
#include "TGenotypeData.h"
#include "TChromosomes.h"
#include "TTask.h"

using namespace GenotypeLikelihoods;

//----------------------------------------------------
// TGlfConverter
// class to converted likelihoods to uint16 and back
//----------------------------------------------------
class TGlfConverter{
private:
	uint16_t maxVal;
	double minLikelihood;
	std::vector<double> likelihoodMap;
	double logOf10DividedByMinus1000 = log(10) / -1000.0;

public:
	TGlfConverter();

	uint16_t maxValue(){ return maxVal; };
	uint16_t toGlfFormat(double scaledLikelihood);
	uint16_t log10ToGlfFormat(double log10ScaledLikelihood);
	uint16_t phredToGlfFormat(uint8_t phred);
	double toScaledLikelihood(uint16_t glfValue);
	double operator[](uint16_t glfValue){ return toScaledLikelihood(glfValue); }
	uint8_t toPhred(uint16_t glfValue);
	double toLog10(uint16_t glfValue);
	double toLog(uint16_t glfValue);
};

//----------------------------------------------------
//TGlfChromosome
//struct to store info on chromosome
//TODO: derive from TChromosome??
//----------------------------------------------------
class TGlfChromosome{
private:
	void setPloidy(uint8_t Ploidy){
		if(Ploidy < 1 || Ploidy > 2)
			throw "Currently GLFs only support ploidies 1 and 2 (not " + toString(Ploidy) + ")!";
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

	TGlfChromosome(std::string Name, uint32_t Length, uint8_t Ploidy){
		name = Name;
		length = Length;
		setPloidy(Ploidy);
	};

	TGlfChromosome(const TGlfChromosome & other){
		update(other);
	};

	void update(const std::string Name, const uint16_t RefId, const uint32_t Length, const uint8_t Ploidy){
		name = Name;
		refId = RefId;
		length = Length;
		setPloidy(Ploidy);
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
	std::string filename;
	gzFile gzfp;
	bool isOpen;
	uint32_t offset;
	std::string version;
	uint8_t zero8, one8;
	uint32_t zero32;
	std::string header;
	uint64_t positionInFile;
	TGlfChromosome curChr;

public:
	TGlfHandle(){
		isOpen = false;
		gzfp = NULL;
		offset = 0;
		version = "GLF2"; //change to next version if older files will not work!
		zero8 = 0;
		one8 = 1;
		zero32 = 0;
		positionInFile = 0;
	};

	void close(){
		if(isOpen){
			gzclose(gzfp);
			isOpen = false;
		}
	};

	std::string name(){
		return filename;
	};

	std::string chr(){
		return curChr.name;
	};

	uint32_t refId(){
		return curChr.refId;
	};

	uint32_t chrLength(){
		return curChr.length;
	};

	bool chrIsHaploid(){
		return curChr.ploidy == 1;
	};

	uint8_t chrPloidy(){
		return curChr.ploidy;
	};

	uint8_t chrNumLikelihoodvalues(){
		return curChr.numLikelihoodValues;
	};
};

//----------------------------------------------------
//TGlfWriter
//----------------------------------------------------
class TGlfWriter:public TGlfHandle{
private:
	long oldPos;
	uint8_t recordType1;
	uint16_t* glfValues; //tmp used for writing
	TGlfConverter converter;

	void init();
	void writeHeader();

	template <typename T>
	inline void write(const T* buf, size_t len){
		positionInFile += gzwrite(gzfp, buf, len);
	};

public:
	TGlfWriter(){
		init();
	};
	TGlfWriter(std::string Filename){
		init();
		open(Filename, "");
	};

	~TGlfWriter(){
		close();
		delete[] glfValues;
	};

	//open & close streams
	void open(std::string Filename);
	void open(std::string Filename, std::string Header);
	void newChromosome(const BAM::TChromosome & chromosome);
	void writeSite(long pos, uint32_t depth, uint8_t RMS_mappingQual, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods);
};


//----------------------------------------------------
//TGlfReader
//----------------------------------------------------
class TGlfReader:public TGlfHandle{
private:
	bool reachedEndOfChr;
	uint32_t HeaderLen;
	uint32_t offset;
	uint8_t tmpInt8;
	int SNPRecordSize;
	uint8_t tmpRecordStorage[19];
	int _lenRead;
	bool _eof;

	int recordType;
	uint32_t _position;
	uint16_t _depth;
	int _RMS_mappingQual;
	uint16_t _genotypeLikelihoodsGLF[10];

	uint16_t* _genotypeLikelihoodsGLF_missingData;
	std::map< uint32_t, TGlfChromosome > _chromosomesAlreadyParsed;

	void init();
	template <typename T>
	inline bool read(T* buf, size_t len){
		_lenRead = gzread(gzfp, buf, len);
		if(!_lenRead) return false;
		positionInFile += _lenRead;
		return true;
	};
	void open();
	bool readChr();
	bool readRecordType();
	void readSNPRecord();
	inline void skipRecord(){
		read(&tmpRecordStorage, SNPRecordSize); //just parse data of one SNP record into garbage
	};

public:
	TGlfReader(){
		init();
	};
	TGlfReader(std::string Filename){
		init();
		open(Filename);
	};
	~TGlfReader(){
		close();
		delete[] _genotypeLikelihoodsGLF_missingData;
	};

	//get details
	bool eof(){ return _eof;};
	TGlfChromosome* pointerToChr(uint32_t refId);
	bool fillPointerToChr(uint32_t refId, TGlfChromosome* & chr);
	uint32_t position(){ return _position; };
	uint16_t depth(){ return _depth; };
	uint16_t* pointerToGenotypeLikelihoodsGLF(){ return _genotypeLikelihoodsGLF; };

	//open file and parse header
	void setFilename(std::string Filename);
	void open(std::string Filename);
	void rewind();
	bool readNext();
	bool jumpToEndOfChr();
	bool jumpToNextChr();
	bool readNextWindow(std::vector<uint16_t*> & genoLikelihoods, const uint32_t refId, const uint32_t start, const uint32_t end);
	void fillGenotypeLikelihoodsGLF(uint16_t* destination);
	void fillGenotypeLikelihoods(double* destination, TGlfConverter* converter);

	//printing
	void printChr();
	void printSite();
	void printToEnd();
};


//------------------------------------------------
// Tasks
//------------------------------------------------
class TTask_printGLF:public TTask{
public:
	TTask_printGLF(){ _explanation = "Printing a GLF file to screen"; };

	void run(TParameters & parameters, TLog* logfile){
		std::string glf = parameters.getParameterString("glf");
		TGlfReader reader(glf);
		reader.printToEnd();
	};
};


#endif /* TGLF_H_ */
