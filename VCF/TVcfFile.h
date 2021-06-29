/*
 * TVcfFile.h
 *
 *  Created on: Aug 8, 2011
 *      Author: wegmannd
 */

#ifndef TVCFFILE_H_
#define TVCFFILE_H_

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "TLog.h"
#include "gzstream.h"
#include "TGenotypeMap.h"
#include "TVcfParser.h"

namespace VCF{

typedef void (TVcfParser::*pt2Function)(TVcfLine &);
//---------------------------------------------------------------------------------------------------------
class TVcfFile_base{
public:
	std::istream* myStream;
	bool inputStreamOpend;
	std::ostream* myOutStream;
	bool outputStreamOpend;
	std::string fileFormat;
	//TVcfColumnNumbers cols;
	TVcfParser parser;
	unsigned int numCols;
	long currentLine;
	std::vector< pt2Function > usedParsers;
	bool automaticallyWriteVcf;
	TVcfLine tempLine;
	bool eof;
	double totalFileSize;
	std::string filename, outputFilename;

	//vector<TVcfFilter> filters;
	//bool applyFilters;

	std::vector<std::string> unknownHeader;

	TVcfFile_base();
	TVcfFile_base(std::string & filename, bool zipped);
	virtual ~TVcfFile_base(){
		if(inputStreamOpend) delete myStream;
		if(outputStreamOpend) delete myOutStream;
	};
	void enableAutomaticWriting(){
		if(!outputStreamOpend)
			throw "Can not automatically write VCF: no output stream has been opened!";
		automaticallyWriteVcf = true;
	};

	void openStream(const std::string & filename);
	void openStream(const std::string & filename, const bool & zipped);
	void openOutputStream(const std::string & filename, const bool & zipped);
	void setOutStream(std::ostream & is);

	//which parsers to use?
	void enablePositionParsing(){usedParsers.push_back(&TVcfParser::parsePosition);};
	void enableVariantQualityParsing(){usedParsers.push_back(&TVcfParser::parseQuality);};
	void enableVariantParsing(){usedParsers.push_back(&TVcfParser::parseVariant);};
	void enableInfoParsing(){usedParsers.push_back(&TVcfParser::parseInfo);};
	void enableFormatParsing(){usedParsers.push_back(&TVcfParser::parseFormat);};
	void enableSampleParsing(){usedParsers.push_back(&TVcfParser::parseSamples);};

	//retrieve info
	int sampleNumber(std::string & Name);
	int numSamples();
	std::string sampleName(unsigned int num);
	bool sampleIsMissing(TVcfLine* line, unsigned int sample);
	bool sampleHasUnknownGenotype(TVcfLine* line, unsigned int sample);
	virtual std::string fieldContentAsString(std::string tag, TVcfLine* line, unsigned int sample);
	virtual int fieldContentAsInt(std::string tag, TVcfLine* line, unsigned int sample);
	virtual int depthAsIntNoCheckForMissingSample(std::string tag, TVcfLine* line, unsigned int sample);

	//modify
	void setSampleMissing(TVcfLine* line, unsigned int sample);
	void setSampleHasUndefinedGenotype(TVcfLine* line, unsigned int sample);
	void updateField(TVcfLine* line, std::string & tag, std::string & Data, unsigned int sample);

	//void addFilter(my_string filter);
	//void filterSamples();
	//void printFilters();

	void parseHeaderVCF_4_0();
	void writeHeaderVCF_4_0();
	void addNewHeaderLine(std::string headerLine);
	bool readLine();
	void addFormat(std::string Line){parser.addFormat(Line);};

	//modify header and columns
	void updateInfo(TVcfLine* line, std::string Id, std::string Data);
	void addToInfo(TVcfLine* line, std::string Id, std::string Data);
	void updatePL(TVcfLine* line, std::string Data, int S);
};

class TVcfFileSingleLine:public TVcfFile_base{
public:
	bool written;

	TVcfFileSingleLine(){written=true;};
	TVcfFileSingleLine(std::string & filename, bool zipped);
	virtual ~TVcfFileSingleLine();
	void writeLine();
	bool next();
	//call specific parsers
	void parseInfo(){ parser.parseSamples(tempLine); };
	void parseFormat(){ parser.parseFormat(tempLine); };
	void parseSamples(){ parser.parseSamples(tempLine); };

	//other stuff
	TVcfLine* pointerToVcfLine(){return &tempLine;};
	void updateInfo(std::string Id, std::string Data);
	void addToInfo(std::string Id, std::string Data);
	void updatePL(std::string Data, unsigned int sample);
	virtual std::string fieldContentAsString(std::string tag, unsigned int sample);
	virtual int fieldContentAsInt(std::string tag, unsigned int sample);
	virtual int depthAsIntNoCheckForMissingSample(std::string tag, unsigned int sample);

	template <typename T>
	std::array<T, 3> genotypeLikelihoods(unsigned int & s){
		return parser.genotypeLikelihoods<T>(tempLine, s);
	};

	//variant info
	uint64_t position();
	uint64_t positionZeroBased();
	std::string chr();
	bool variantQualityIsMissing();
	double variantQuality();
	int getNumAlleles();
	std::string getRefAllele();
	std::string getFirstAltAllele();
	std::string getAllele(int num);
	bool isBialleleicSNP();

	//sample info
	void setSampleMissing(unsigned int sample);
	void setSampleHasUndefinedGenotype(unsigned int sample);
	void updateField(std::string tag, std::string & Data, unsigned int sample);
	bool sampleIsMissing(unsigned int sample);
	bool sampleHasUndefinedGenotype(unsigned int sample);
	bool sampleIsHaploid(unsigned int sample);
	bool sampleIsDiploid(unsigned int sample);
	bool sampleIsHomoRef(unsigned int sample);
	bool sampleIsHeteroRefNonref(unsigned int sample);
	genometools::Base getFirstAlleleOfSample(unsigned int num);
	genometools::Base getSecondAlleleOfSample(unsigned int num);
	genometools::BiallelicGenotype sampleBiallelicGenotype(const unsigned int & num);
	genometools::Genotype sampleGenotype(const unsigned int & num);
	float sampleGenotypeQuality(unsigned int sample);
	double sampleDepth(unsigned int sample);
	// int sampleDepth(unsigned int sample);
	bool formatColExists(std::string tag){ return parser.formatColExists(tag, tempLine); };
	std::string getSampleContentAt(std::string tag, unsigned int sample){
		return parser.sampleContentAt(tempLine, tag, sample);
	}
};

}; //end namespace

#endif /* TVCFFILE_H_ */
