/*
 * TGlfMultiReader.h
 *
 *  Created on: Feb 11, 2020
 *      Author: wegmannd
 */

#ifndef GLF_TGLFMULTIREADER_H_
#define GLF_TGLFMULTIREADER_H_

#include "TGLF.h"
#include "TFastaBuffer.h"

//----------------------------------------------------
//TMultiGLFData
//----------------------------------------------------
struct TMultiGLFDataSample{
	uint16_t* genotypeLikelihoodsGLF;
	bool hasData;
	bool isHaploid;
	uint16_t depth;
};

class TMultiGLFData{
private:
	GenotypeLikelihoods::TGenotypeMap genoMap;
	TGlfConverter converter;

	void _free();

public:
	uint32_t size;
	TMultiGLFDataSample* samples;

	TMultiGLFData();
	TMultiGLFData(int Size);
	~TMultiGLFData(){
		_free();
	};

	void resize(uint32_t Size);
	void fill(PopulationTools::TPopulationLikehoodLocus & storage, const int alleleicCombination);
	uint32_t totalDepth();
};

//----------------------------------------------------
//TGlfMultiReaderVcfLocusDefinition
//----------------------------------------------------
class TGlfMultiReaderVcf{
private:
	std::vector<std::string> genotypeStrings = {"0", "1", "0/0", "0/1", "1/1"};
	bool _usePhredScaledLikelihoods;

	TRandomGenerator* randomGenerator;
	GenotypeLikelihoods::TGenotypeMap genoMap;
	TGlfConverter converter;
	gz::ogzstream vcf;
	bool vcfOpened;

	Base ref, alt;
	char ref_char, alt_char;
	std::string genotypeString[5];
	uint8_t refHomIndex, hetIndex, altHomIndex;

	void _openVCF(const std::string & filename, const std::string & source, std::vector<std::string> & sampleNames);
	void _closeVCF();
	void _setMajorMinor(const Base & refAllele, const Base & altAllele);
	void writeLikelihood(const uint16_t & likGlf);
	void writeDiploidIndividualToVCF(TMultiGLFDataSample & sample);
	void writeHaploidIndividualToVCF(TMultiGLFDataSample & sample);

public:
	TGlfMultiReaderVcf(const std::string filename, const std::string source, std::vector<std::string> & sampleNames, TRandomGenerator* RandomGenerator);
	~TGlfMultiReaderVcf(){
		_closeVCF();
	}

	void usePhredScaledLikelihoods(){ _usePhredScaledLikelihoods = true; };
	void writeSite(const std::string & chrName, const uint32_t & position, const int & varianTQuality, TMultiGLFData & data, const Base Ref, const Base Alt);
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
	GenotypeLikelihoods::TGenotypeMap genoMap;

	void _openGLFs(TLog* logfile);

	//active files
	//Object will loop only over active files
	bool _onlyJumpToPositionsWithData;
	int numActiveFiles;
	bool* GLFIsActive;
	std::vector<int> activeGLFs;
	std::vector<TGlfReader*> pointerToActiveGLFs;
	int _getGLFIndexFromName(const std::string & name);
	void _setActive(const int index);
	void _setAllInactive();
	void _prepareParsing();
	bool _jumpToNextPositionWithData();
	void _updateChromosomeInfo();

	//Moving along active files
	uint32_t _position, _nextPosition; //next is anticipated position, used to advance
	uint32_t _curRefId;
	TGlfChromosome _curChr;
	int _numActiveFilesWithData;
	uint16_t genotypeQualitiesMissingData[10];
	int minDepth;

	//reference
	bool hasReference;
	BAM::TFastaBuffer fastaBuffer;

	bool moveToNextChromosome();

	void writeDiploidIndividualToVCF(const int & ind, gz::ogzstream & vcf, const Base & major, const Base & minor, const std::vector<std::string> & genotypeStrings, TRandomGenerator* randomGenerator, const bool & usePhredLikelihoods);
	void writeHaploidIndividualToVCF(const int & ind, gz::ogzstream & vcf, const Base & major, const Base & minor, const std::vector<std::string> & genotypeStrings, TRandomGenerator* randomGenerator, const bool & usePhredLikelihoods);

public:
	TGlfConverter converter;
	TMultiGLFData data;

	TGlfMultiReader();
	TGlfMultiReader(std::vector<std::string> FileNames, TLog* logfile);
	TGlfMultiReader(TParameters & params, TLog* logfile);
	void init();

	~TGlfMultiReader();

	void openGLFs(const std::vector<std::string> & Filenames, TLog* logfile);
	void openGLFs(TParameters & params, TLog* logfile);
	void closeGLF();
	void setDepthFilter(int MinDepth, TLog* logfile);
	void addReference(const std::string FastaFile);
	void onlyJumpToPositionsWithData(){ _onlyJumpToPositionsWithData = true; };

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
	void writeSampleNamesOfActiveFiles(gz::ogzstream & out, std::string sep);
	void fillSampleNamesOfActiveFiles(std::vector<std::string> & vec);

	//access data
	int numSamples(){ return numGLFs; };
	int numActiveSamples(){ return numActiveFiles; };
	int numActiveSamplesWithData(){ return _numActiveFilesWithData; };
	std::string chr(){return _curChr.name; };
	uint32_t position(){return _position; };
	Base refBase();
};



#endif /* GLF_TGLFMULTIREADER_H_ */
