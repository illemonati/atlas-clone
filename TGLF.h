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
#include <vector>
//#include <bitset>

//----------------------------------------------------
//TGlfHandle
//----------------------------------------------------
class TGlfHandle{
protected:
	std::string filename;
	gzFile gzfp;
	bool isOpen;
	long positionInFile;
	uint32_t offset;
	std::string version;
	uint8_t zero8, one8;
	uint32_t zero32;
	std::string curChr;
	std::string header;

public:
	TGlfHandle(){
		isOpen = false;
		gzfp = NULL;
		positionInFile = 0;
		offset = 0;
		version = "GLF3";
		zero8 = 0;
		one8 = 1;
		zero32 = 0;
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
		return curChr;
	};
};

//----------------------------------------------------
//TGlfWriter
//----------------------------------------------------
class TGlfWriter:public TGlfHandle{
private:
	long oldPos;
	uint8_t recordType1;

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
	void init(){
		oldPos = 0;
		recordType1 = one8 << 4;
	};
	~TGlfWriter(){
		close();
	};

	//open & close streams
	void open(std::string Filename, std::string Header);
	void newChromosome(std::string name, uint32_t length);
	void writeSite(long pos, uint32_t depth, uint8_t RMS_mappingQual, uint8_t* genoQualities, uint32_t & maxLL);
};


//----------------------------------------------------
//TGlfReader
//----------------------------------------------------
class TGlfReader:public TGlfHandle{
private:
	bool reachedEndOfChr;
	uint32_t HeaderLen;
	uint32_t offset;
	uint32_t depth_mask;
	uint32_t tmpInt32;
	uint8_t tmpInt8;
	uint8_t tmpInt8_10[10];
	uint8_t tmpRecordStorage[19];
	long _chrLength;
	int _lenRead;
	int _i;
	bool _eof;
	int* genotypeQualitiesMissingData;
	std::vector<std::string> chromosomesAlreadyParsed;

	void init();
	template <typename T>
	inline bool read(T* buf, size_t len){
		_lenRead = gzread(gzfp, buf, len);
		if(!_lenRead) return false;
		positionInFile += _lenRead;
		return true;
	};
	bool readChr();
	bool chromosomeParsed(std::string & chr);
	bool readRecordType();
	void readSNPRecord();
	inline void skipRecord(){
		read(&tmpRecordStorage, 152); //just parse data of one SNP record into garbage
	};

public:
	int recordType;
	long position;
	int maxLL;
	int depth;
	int RMS_mappingQual;
	int genotypeQualities[10];

	TGlfReader(){
		init();
	};
	TGlfReader(std::string Filename){
		init();
		open(Filename);
	};
	~TGlfReader(){
		close();
		delete[] genotypeQualitiesMissingData;
	};

	//get detailes
	long chrLength(){ return _chrLength; };
	bool eof(){ return _eof;};

	//open file and parse header
	void open(std::string Filename);
	void rewind();
	bool readNext();
	bool jumpToEndOfChr();
	bool readNextWindow(std::vector<int*> & genoLikelihoods, std::string chr, long start, long end);
	void fillGenotypeQualities(int* destination);

	//printing
	void printChr();
	void printSite();
	void printToEnd();

};

#endif /* TGLF_H_ */
