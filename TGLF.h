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
#include "TParameters.h"
#include "stringFunctions.h"
#include "TGenotypeMap.h"
#include "TRandomGenerator.h"


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
	std::string curChr;
	int curChrNumber;
	std::string header;
	long positionInFile;

public:
	TGlfHandle(){
		isOpen = false;
		gzfp = NULL;
		offset = 0;
		version = "GLF3";
		zero8 = 0;
		one8 = 1;
		zero32 = 0;
		curChrNumber = 0;
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
		return curChr;
	};

	int chrNumber(){
		return curChrNumber;
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
	int SNPRecordSize;
	uint8_t tmpRecordStorage[19];
	long _chrLength;
	int _lenRead;
	bool _eof;
	uint8_t genotypeQualitiesMissingData[10];
	std::vector< std::pair<std::string, long> > chromosomesAlreadyParsed;

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
	bool chromosomeParsed(std::string & chr);
	bool readRecordType();
	void readSNPRecord();
	inline void skipRecord(){
		read(&tmpRecordStorage, SNPRecordSize); //just parse data of one SNP record into garbage
	};

public:
	int recordType;
	long position;
	int maxLL;
	int depth;
	int RMS_mappingQual;
	uint8_t genotypeQualities[10];

	TGlfReader(){
		init();
	};
	TGlfReader(std::string Filename){
		init();
		open(Filename);
	};
	~TGlfReader(){
		close();
		//delete[] genotypeQualitiesMissingData;
	};

	//get details
	long chrLength(){ return _chrLength; };
	bool eof(){ return _eof;};
	std::string getNameOfParsedChr(int chrNumber);
	long getLengthOfParsedChr(int chrNumber);

	//open file and parse header
	void setFilename(std::string Filename);
	void open(std::string Filename);
	void rewind();
	bool readNext();
	bool jumpToEndOfChr();
	bool jumpToNextChr();
	bool readNextWindow(std::vector<uint8_t*> & genoLikelihoods, std::string chr, long start, long end);
	void fillGenotypeQualities(uint8_t* destination);

	//printing
	void printChr();
	void printSite();
	void printToEnd();

};

//----------------------------------------------------
//TGlfMultiReader
//----------------------------------------------------
class TGlfMultiReader{
private:
	int numGLFs;
	std::vector<std::string> GLFNames;
	TGlfReader* GLFs;
	bool readersOpened;
	TGenotypeMap genoMap;

	void _openGLFs(TLog* logfile);

	//active files
	//Object will loop only over active files
	int numActiveFiles;
	bool* GLFIsActive;
	std::vector<int> activeGLFs;
	std::vector<TGlfReader*> pointerToActiveGLFs;
	int _getGLFIndexFromName(const std::string & name);
	void _setActive(const int index);
	void _setAllInactive();
	int _minChrNumberActiveFiles();
	void _setCurChrName();
	void _prepareParsing();

	//Moving along active files
	long _position;
	int _curChrNumber;
	long _curChrLength;
	std::string _curChrName;
	int _numActiveFilesWithData;
	uint8_t genotypeQualitiesMissingData[10];

	bool moveToNextChromosome();

public:
	uint8_t** data;
	bool* hasData;
	bool dataInitialized;

	TGlfMultiReader();
	TGlfMultiReader(std::vector<std::string> FileNames, TLog* logfile);
	TGlfMultiReader(TParameters & params, TLog* logfile);
	void init();

	~TGlfMultiReader();

	void openGLFs(const std::vector<std::string> & Filenames, TLog* logfile);
	void openGLFs(TParameters & params, TLog* logfile);
	void closeGLF();

	//set active / inactive
	void setActive(const int index);
	void setActive(const std::string & name);
	void setActive(const int index1, const int index2);
	void setActive(const std::string & name1, const std::string & name2);
	void setActive(std::vector<int> & indexes);
	void setActive(std::vector<std::string> & names);
	void setAllActive();

	//access data
	int numSamples(){ return numGLFs; };
	int numActiveSamples(){ return numActiveFiles; };
	int chrNumber(){return _curChrNumber;};
	std::string chr(){return _curChrName;};
	long position(){return _position;};

	//parse
	bool readNext();
	void print();
	void writeSampleNamesOfActiveFiles(gz::ogzstream & out, std::string sep);
	void writeVCFHeader(gz::ogzstream & vcf);
	void writeSiteToVCF(gz::ogzstream & vcf, int variantQuality, int & refHomIndex, int & hetIndex, int & altHomIndex, TRandomGenerator & randomGenerator);

};

#endif /* TGLF_H_ */
