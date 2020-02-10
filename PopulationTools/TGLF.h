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
#include "../gzstream.h"
#include <algorithm>
#include <string.h>
#include <TPopulationLikelihoodLocus.h>
#include <vector>
#include "../TParameters.h"
#include "../stringFunctions.h"
#include "../TGenotypeMap.h"
#include "../TRandomGenerator.h"
#include "../TFastaBuffer.h"

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
	uint32_t number;
	std::string name;
	uint16_t refId;
	uint32_t length;
	uint8_t ploidy;
	uint8_t numLikelihoodValues; //depends on ploidy
	uint8_t maxNumLikelihoodValues; //maximum possible

	TGlfChromosome(){
		number = 0;
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
		number = 1;
	};

	TGlfChromosome(const TGlfChromosome & other){
		name = other.name;
		refId = other.refId;
		length = other.length;
		ploidy = other.ploidy;
		number = other.number;
		numLikelihoodValues = other.numLikelihoodValues;
		maxNumLikelihoodValues = other.maxNumLikelihoodValues;
	};

	void update(const std::string Name, const uint16_t RefId, const uint32_t Length, const uint8_t Ploidy){
		name = Name;
		refId = RefId;
		length = Length;
		setPloidy(Ploidy);
		++number;
	};

	void clear(){
		name = "";
		refId = 0;
		number = 0;
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
	long positionInFile;
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

	int chrNumber(){
		return curChr.number;
	};

	uint16_t refId(){
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
	void open(std::string Filename, std::string Header);
	void newChromosome(const std::string name, const uint16_t refId, const uint32_t length, const uint8_t ploidy);
	void writeSite(long pos, uint32_t depth, uint8_t RMS_mappingQual, double* genotypeLikelihoods);
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
	uint16_t* genotypeQualitiesMissingData;
	std::vector< TGlfChromosome > chromosomesAlreadyParsed;

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
	//uint16_t maxLL;
	uint16_t depth;
	int RMS_mappingQual;
	uint16_t genotypeQualities[10];

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

	//get details
	bool eof(){ return _eof;};
	std::string getNameOfParsedChr(uint32_t chrNumber);
	uint16_t getRefIdOfParsedChr(uint32_t chrNumber);
	uint32_t getLengthOfParsedChr(uint32_t chrNumber);

	//open file and parse header
	void setFilename(std::string Filename);
	void open(std::string Filename);
	void rewind();
	bool readNext();
	bool jumpToEndOfChr();
	bool jumpToNextChr();
	bool readNextWindow(std::vector<uint16_t*> & genoLikelihoods, std::string chr, long start, long end);
	void fillGenotypeQualities(uint16_t* destination);
	void fillGenotypeLikelihoods(double* destination, TGlfConverter* converter);

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
	TGlfConverter converter;

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
	int32_t _position; //needs to be signed as it is initialized to -1 such that the first position will be 0
	uint16_t _curChrNumber;
	uint16_t _curRefId;
	uint32_t _curChrLength;
	std::string _curChrName;
	int _numActiveFilesWithData;
	uint16_t genotypeQualitiesMissingData[10];
	int minDepth;

	//reference
	bool hasReference;
	TFastaBuffer fastaBuffer;

	bool moveToNextChromosome();

	void writeDiploidIndividualToVCF(const int & ind, gz::ogzstream & vcf, const Base & major, const Base & minor, const std::vector<std::string> & genotypeStrings, TRandomGenerator* randomGenerator, const bool & usePhredLikelihoods);
	void writeHaploidIndividualToVCF(const int & ind, gz::ogzstream & vcf, const Base & major, const Base & minor, const std::vector<std::string> & genotypeStrings, TRandomGenerator* randomGenerator, const bool & usePhredLikelihoods);

public:
	uint16_t** data;
	bool* hasData;
	bool* isHaploid;
	bool dataInitialized;

	TGlfMultiReader();
	TGlfMultiReader(std::vector<std::string> FileNames, TLog* logfile);
	TGlfMultiReader(TParameters & params, TLog* logfile);
	void init();

	~TGlfMultiReader();

	void openGLFs(const std::vector<std::string> & Filenames, TLog* logfile);
	void openGLFs(TParameters & params, TLog* logfile);
	void closeGLF();
	void setDepthFilter(int MinDepth, TLog* logfile);
	void addReference(BamTools::Fasta* Reference);

	//set active / inactive
	void setActive(const int index);
	void setActive(const std::string & name);
	void setActive(const int index1, const int index2);
	void setActive(const std::string & name1, const std::string & name2);
	void setActive(std::vector<int> & indexes);
	void setActive(std::vector<std::string> & names);
	void setAllActive();

	//parse
	bool readNext();
	void print();
	void fill(TPopulationLikehoodLocus & data, const int alleleicCombination);
	void writeSampleNamesOfActiveFiles(gz::ogzstream & out, std::string sep);
	void writeVCFHeader(gz::ogzstream & vcf, bool usePhredLikelihoods);
	void writeSiteToVCF(gz::ogzstream & vcf, const int & varianTQuality, const Base major, const Base minor, TRandomGenerator* randomGenerator, const bool & usePhredLikelihoods);

	//access data
	int numSamples(){ return numGLFs; };
	int numActiveSamples(){ return numActiveFiles; };
	int numActiveSamplesWithData(){ return _numActiveFilesWithData; };
	int chrNumber(){return _curChrNumber;};
	std::string chr(){return _curChrName;};
	uint32_t position(){return _position;};
	Base refBase();
};

#endif /* TGLF_H_ */
