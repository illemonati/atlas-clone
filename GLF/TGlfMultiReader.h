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
#include "TParameters.h"

namespace GLF{

using coretools::TRandomGenerator;
using coretools::TParameters;
using coretools::TLog;

//----------------------------------------------------
//TMultiGLFData
//----------------------------------------------------
struct TMultiGLFDataSample{
	genometools::HighPrecisionPhredIntProbability* genotypeLikelihoodsGLF; //points to data TGlfReader
	bool hasData;
	bool isHaploid;
	uint16_t depth;
};

class TMultiGLFData{
public:
	std::vector<TMultiGLFDataSample> samples;

	TMultiGLFData() = default;
	TMultiGLFData(const uint32_t & Size);
	~TMultiGLFData(){};

	void resize(const uint32_t & Size);
	size_t size() const { return samples.size(); };
	void fill(PopulationTools::TPopulationLikehoodLocus & storage, const genometools::AllelicCombination & alleleicCombination) const;
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
	gz::ogzstream vcf;
	bool vcfOpened;

	genometools::Base _ref, _alt;
	//char ref_char, alt_char;
	std::string _genotypeString[5];
	genometools::Genotype _refHom, _het, _altHom;

	void _openVCF(const std::string & filename, const std::string & source, std::vector<std::string> & sampleNames);
	void _closeVCF();
	void _setMajorMinor(const genometools::Base & refAllele, const genometools::Base & altAllele);
	void writeLikelihood(const genometools::HighPrecisionPhredIntProbability & likGlf);
	void writeDiploidIndividualToVCF(TMultiGLFDataSample & sample);
	void writeHaploidIndividualToVCF(TMultiGLFDataSample & sample);

public:
	TGlfMultiReaderVcf(const std::string filename, const std::string source, std::vector<std::string> & sampleNames, TRandomGenerator* RandomGenerator);
	~TGlfMultiReaderVcf(){
		_closeVCF();
	}

	void usePhredScaledLikelihoods(){ _usePhredScaledLikelihoods = true; };
	void writeSite(const std::string & chrName, const uint32_t & position, const genometools::PhredIntProbability & varianTQuality, TMultiGLFData & data, const genometools::Base Ref, const genometools::Base Alt);
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
	genometools::HighPrecisionPhredIntProbability genotypeQualitiesMissingData[10];
	int minDepth;

	//reference
	bool hasReference;
	BAM::TFastaBuffer fastaBuffer;

	bool moveToNextChromosome();

	void writeDiploidIndividualToVCF(const int & ind, gz::ogzstream & vcf, const genometools::Base & major, const genometools::Base & minor, const std::vector<std::string> & genotypeStrings, TRandomGenerator* randomGenerator, const bool & usePhredLikelihoods);
	void writeHaploidIndividualToVCF(const int & ind, gz::ogzstream & vcf, const genometools::Base & major, const genometools::Base & minor, const std::vector<std::string> & genotypeStrings, TRandomGenerator* randomGenerator, const bool & usePhredLikelihoods);

public:
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
	void onlyJumpToPositionsWithData(const bool & set = true){ _onlyJumpToPositionsWithData = set; };


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
	genometools::Base refBase();
};

}; //end namespace GLF

#endif /* GLF_TGLFMULTIREADER_H_ */
