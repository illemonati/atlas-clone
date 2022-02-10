/*
 * TGlfMultiReader.h
 *
 *  Created on: Feb 11, 2020
 *      Author: wegmannd
 */

#ifndef GLF_TGLFMULTIREADER_H_
#define GLF_TGLFMULTIREADER_H_


#include "GenotypeTypes.h"
#include "TGLF.h"
#include "TFastaBuffer.h"
#include "TParameters.h"
#include "stringFunctions.h"

namespace GLF{

using coretools::TRandomGenerator;
using coretools::TParameters;
using coretools::TLog;
using genometools::HighPrecisionPhredIntProbability;

//----------------------------------------------------
// TMultiGLFDataSample
//----------------------------------------------------
class TMultiGLFDataSample{
private:
	HighPrecisionPhredIntProbability* _genotypeLikelihoodsGLF; //points to data TGlfReader
	bool _hasData;
	bool _isHaploid;
	uint16_t _depth;

public:
	TMultiGLFDataSample(){
		_hasData = false;
		_isHaploid = false;
		_depth = 0;
		_genotypeLikelihoodsGLF = nullptr;
	};

	void setDiploid(HighPrecisionPhredIntProbability* GLs, uint16_t Depth){
		_genotypeLikelihoodsGLF = GLs;
		_hasData = true;
		_isHaploid = false;
		_depth = Depth;
	};

	void setHaploid(HighPrecisionPhredIntProbability* GLs, uint16_t Depth){
		_genotypeLikelihoodsGLF = GLs;
		_hasData = true;
		_isHaploid = true;
		_depth = Depth;
	};

	void setMissingDiploid(){
		_hasData = false;
		_depth = 0;
		_isHaploid = false;
	};

	void setMissingHaploid(){
		_hasData = false;
		_depth = 0;
		_isHaploid = true;
	};

	uint16_t depth() const { return _depth; };
	bool hasData() const { return _hasData; };

	bool isHaploid() const {
		return _isHaploid;
	};

	const HighPrecisionPhredIntProbability operator[](const genometools::Genotype & G) const {
		if(!_hasData){
			return HighPrecisionPhredIntProbability::highest();
		}
		if(_isHaploid){
			throw std::runtime_error("HighPrecisionPhredIntProbability TMultiGLFDataSample::operator[](const genometools::Genotype & G): sample is haploid!");
		}
		return _genotypeLikelihoodsGLF[genometools::index(G)];
	};

	const HighPrecisionPhredIntProbability operator[](const genometools::Base & B) const{
		if(!_hasData){
			return HighPrecisionPhredIntProbability::highest();
		}
		if(!_isHaploid){
			throw std::runtime_error("HighPrecisionPhredIntProbability TMultiGLFDataSample::operator[](const genometools::Base & B): sample is diploid!");
		}
		return _genotypeLikelihoodsGLF[genometools::index(B)];
	};
};

//-------------------------------------
// TGenotypeLikelihoodsOneAllelicCombination
//-------------------------------------
class TMultiGLFDataSampleOneAllelicCombination{
private:
	bool _isMissing;
	bool _isHaploid;
	std::array<HighPrecisionPhredIntProbability, 3> _GLs;
	constexpr static const HighPrecisionPhredIntProbability _maxGTL = HighPrecisionPhredIntProbability::highest();

public:
	TMultiGLFDataSampleOneAllelicCombination() = default;
	~TMultiGLFDataSampleOneAllelicCombination() = default;

	void setDiploid(const HighPrecisionPhredIntProbability & GL_homoFirst, const HighPrecisionPhredIntProbability & GL_het, const HighPrecisionPhredIntProbability & GL_homoSecond){
		_GLs[0] = GL_homoFirst;
		_GLs[1] = GL_het;
		_GLs[2] = GL_homoSecond;
		_isHaploid = false;
		_isMissing = false;
	};

	void setHaploid(const HighPrecisionPhredIntProbability & GL_first, const HighPrecisionPhredIntProbability & GL_second){
		_GLs[0] = GL_first;
		_GLs[1] = GL_second;
		_isHaploid = true;
		_isMissing = false;
	};

	void setMissingDiploid(){
		_isHaploid = false;
		_isMissing = true;
	};

	void setMissingHaploid(){
		_isHaploid = true;
		_isMissing = true;
	};

	bool isMissing() const {
		return _isMissing;
	};

	bool isHaploid() const {
		return _isHaploid;
	};

	constexpr const HighPrecisionPhredIntProbability& operator[](const genometools::BiallelicGenotype & Genotype) const {
		 //check
		 if(_isHaploid){
			 if(!genometools::isHaploid(Genotype)){
				throw std::runtime_error("constexpr const HighPrecisionPhredIntProbability& TGenotypeLikelihoodsOneAllelicCombination::operator[](const genometools::BiallelicGenotype & Genotype) const: sample is haploid!");
			 }
		 } else {
			 if(!genometools::isDiploid(Genotype)){
				throw std::runtime_error("constexpr const HighPrecisionPhredIntProbability& TGenotypeLikelihoodsOneAllelicCombination::operator[](const genometools::BiallelicGenotype & Genotype) const: sample is Diploid!");
			 }
		 }

		 //return
		 if(_isMissing){
			 return _maxGTL;
		 } else {
			 return _GLs[genometools::altAlleleCounts(Genotype)];
		 }
	};
};

using TMultiGLFDataOneAllelicCombination = std::vector<TMultiGLFDataSampleOneAllelicCombination>;

//------------------------------------------------
// TMultiGLFData
//------------------------------------------------
class TMultiGLFData{
private:
	std::vector<TMultiGLFDataSample> samples;

public:
	TMultiGLFData() = default;
	TMultiGLFData(uint32_t Size);
	~TMultiGLFData(){};

	void resize(uint32_t Size);
	size_t size() const { return samples.size(); };
	TMultiGLFDataSample& operator[](uint32_t Sample){ return samples[Sample]; };
	const TMultiGLFDataSample& operator[](uint32_t Sample) const { return samples[Sample]; };
	void fill(TMultiGLFDataOneAllelicCombination & storage, const genometools::AllelicCombination & alleleicCombination) const;
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
	void writeLikelihood(const HighPrecisionPhredIntProbability & likGlf);
	void writeDiploidIndividualToVCF(TMultiGLFDataSample & sample);
	void writeHaploidIndividualToVCF(TMultiGLFDataSample & sample);

public:
	TGlfMultiReaderVcf(const std::string filename, const std::string source, std::vector<std::string> & sampleNames, const bool & UsePhredScaledLikelihoods, TRandomGenerator* RandomGenerator);
	~TGlfMultiReaderVcf(){
		_closeVCF();
	}

	void writeSite(const std::string & chrName, uint32_t position, const genometools::PhredIntProbability & varianTQuality, TMultiGLFData & data, const genometools::Base Ref, const genometools::Base Alt);
};

//----------------------------------------------------
//TGlfMultiReader
//----------------------------------------------------
class TGlfMultiReader{
private:
	uint8_t numGLFs;
	std::vector<std::string> GLFNames;
	TGlfReader* GLFs;
	bool readersOpened;

	void _openGLFs(TLog* logfile);

	//active files
	//Object will loop only over active files
	bool _onlyJumpToPositionsWithData;
	uint32_t numActiveFiles;
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
	uint32_t _numActiveFilesWithData;
	uint32_t minDepth;

	//reference
	bool hasReference;
	BAM::TFastaBuffer fastaBuffer;

	bool moveToNextChromosome();

	void writeDiploidIndividualToVCF(int ind, gz::ogzstream & vcf, const genometools::Base & major, const genometools::Base & minor, const std::vector<std::string> & genotypeStrings, TRandomGenerator* randomGenerator, const bool & usePhredLikelihoods);
	void writeHaploidIndividualToVCF(int ind, gz::ogzstream & vcf, const genometools::Base & major, const genometools::Base & minor, const std::vector<std::string> & genotypeStrings, TRandomGenerator* randomGenerator, const bool & usePhredLikelihoods);

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

	//output
	void print();
	void writeSampleNamesOfActiveFiles(gz::ogzstream & out, std::string sep);
	std::vector<std::string> namesOfActiveFiles();
	std::vector<std::string> sampleNamesOfActiveFiles();

	//access data
	uint32_t numSamples(){ return numGLFs; };
	uint32_t numActiveSamples(){ return numActiveFiles; };
	uint32_t numActiveSamplesWithData(){ return _numActiveFilesWithData; };
	std::string chr(){return _curChr.name; };
	uint32_t position(){return _position; };
	genometools::Base refBase();
};

}; //end namespace GLF

#endif /* GLF_TGLFMULTIREADER_H_ */
