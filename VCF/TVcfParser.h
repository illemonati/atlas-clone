/*
 * TVcfParser.h
 *
 *  Created on: Jun 15, 2011
 *      Author: wegmannd
 */

#ifndef TVCFPARSER_H_
#define TVCFPARSER_H_

#include "stringFunctions.h"
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <math.h>

namespace VCF{

//TODO: use header info to check entries

	enum VCF_TYPE {UNKNOWN, INTEGER, FLOAT, FLAG, CHAR,  STRING};

struct GTLikelihoods{
	float AA;
	float AB;
	float BB;
};
class TVcfColumnNumbers{
public:
	int Chr, Pos, Id, Ref, Alt, Qual, Filter, Info, Format, FirstInd;
	TVcfColumnNumbers(){
		FirstInd=999999;
		Chr=-1;
		Pos=-1;
		Id=-1;
		Ref=-1;
		Alt=-1;
		Qual=-1;
		Filter=-1;
		Info=-1;
		Format=-1;
	};

	void set(std::string & tag, int & i){
		if(stringContains(tag, "CHROM")) Chr=i;
		else if(stringContains(tag, "POS")) Pos=i;
		else if(stringContains(tag, "ID")) Id=i;
		else if(stringContains(tag, "REF")) Ref=i;
		else if(stringContains(tag, "ALT")) Alt=i;
		else if(stringContains(tag, "QUAL")) Qual=i;
		else if(stringContains(tag, "FILTER")) Filter=i;
		else if(stringContains(tag, "INFO")) Info=i;
		else if(stringContains(tag, "FORMAT")){
			Format=i;
			//next is first individual!
			FirstInd=i+1;
		}
	};
	void check(){
		if(Chr<0) throw "Error when reading vcf header: column 'CHROM' is missing!";
		if(Pos<0) throw "Error when reading vcf header: column 'POS' is missing!";
		if(Id<0) throw "Error when reading vcf header: column 'ID' is missing!";
		if(Ref<0) throw "Error when reading vcf header: column 'REF' is missing!";
		if(Alt<0) throw "Error when reading vcf header: column 'ALT' is missing!";
		if(Qual<0) throw "Error when reading vcf header: column 'QUAL' is missing!";
		if(Filter<0) throw "Error when reading vcf header: column 'FILTER' is missing!";
		if(Info<0) throw "Error when reading vcf header: column 'INFO' is missing!";
		if(Format<0) throw "Error when reading vcf header: column 'FORMAT' is missing!";
	};
};
//---------------------------------------------------------------------------------------------------------
class TVcfHeaderLine{
public:
	std::string _id;
	int number;
	std::string numberString;
	VCF_TYPE type;
	std::string typeString;
	std::string desc;

	TVcfHeaderLine(){init();};
	TVcfHeaderLine(std::string & Line);
	TVcfHeaderLine(std::string & ID, std::string & Number, VCF_TYPE & Type, std::string & Desc);
	void init();
	void update(std::string & Number, VCF_TYPE & Type, std::string & Desc);
	VCF_TYPE getTypefromString(std::string & s);
	std::string getStringfromType(VCF_TYPE & type);
	std::string getString();
};
//---------------------------------------------------------------------------------------------------------
//FILTER _> TODOD: integrate with header rows and FILTER column -> see TVcfInfo
/*
class TVcfFilter{
public:
	long* currentLine;
	std::string tag;
	std::string subTag;
	bool larger;
	bool sub;
	float val;

	TVcfFilter(std::string filter, long* CurrentLine);
	bool pass(TVcfFormat* format, vector<std::string>* data);
	void print();
};
*/
//---------------------------------------------------------------------------------------------------------
class TVcfSample{
private:

public:
	std::vector<std::string> data;
	std::pair<int, int> genotype;
	bool missing;
	bool hasGenotype;
	bool unknownGenotype;
	bool isHaploid;

	TVcfSample();
	void readData(std::string s);
	void addData(std::string d);
	void updateData(int pos, std::string d);
	void setGenotype(int firstAllele, int secondAllele);
	void setGenotype(int haploidAllele);
	void setMissingGenotype();
	bool parse(std::string s, const int genotypeCol);
	bool checkGenotype(int max);
	std::string getCol(const int col);
	void write(std::ostream & out, unsigned int numFields);
};
//---------------------------------------------------------------------------------------------------------

class TVcfLine{
public:
	long lineNumber;
	std::vector<std::string> data; //used to store read data
	bool positionParsed, variantParsed, idParsed, filterParsed, qualityParsed, infoParsed, formatParsed, samplesParsed;
	uint64_t pos;
	std::string chr;
	double variantQuality;
	bool variantQualityMissing;

	std::vector<std::string> variants; //entry at 0 is reference
	std::map<std::string, std::vector<std::string> > info;
	std::vector<std::string> formatOrdered;
	std::map<std::string, int> format;
	std::vector<TVcfSample> samples;

	std::string id, qual, filter;

    //ID, FILTER and QUAL: these fields are currently NOT PARSED -> TODO
	std::vector<char> bases; //entry at 0 is ref

	TVcfLine();
	TVcfLine(std::string & line, unsigned int & numCols, long & LineNumber);
	~TVcfLine(){};
	void update(std::string & line, unsigned int & numCols, long & LineNumber);
	bool variantExists(const std::string & var);
	bool addVariant(const std::string & var);
	void writeVariant(std::ostream & out);
};

class TVcfParser{
private:
	void savePhredScore(std::string & phredString, uint8_t & phred);
	double readGL(std::string & GLString);
	void saveGLAsPhredScore(std::string & GLString, uint8_t & phred);

public:
	TVcfColumnNumbers cols;
	std::map<std::string, TVcfHeaderLine> info;
	std::map<std::string, TVcfHeaderLine> format;
	std::vector<std::string> samples;
	int maxIndColPlusOne;
	std::vector<TVcfSample>::iterator lineSampleIt;
	std::string genotypeTag;

	TVcfParser(){
		genotypeTag="GT";
		maxIndColPlusOne=-1;
	};

	//parsers
	void parsePosition(TVcfLine & line);
	void parseVariant(TVcfLine & line);
	void parseQuality(TVcfLine & line);
	void parseFormat(TVcfLine & line);
	void parseInfo(TVcfLine & line);
	void parseSamples(TVcfLine & line);

	//other functions
	int getFormatCol(std::string & tag, TVcfLine & line);
	int getFormatCol(TVcfLine & line, std::string tag){return getFormatCol(tag, line);};
	bool formatColExists(std::string & tag, TVcfLine & line);
	int addFormatCol(std::string & tag, TVcfLine & line);
//	void setColNumbers(TVcfColumnNumbers* Cols){cols=Cols;};
	void addInfo(std::string & Line);
	void updateInfo(std::string ID, std::string Number, VCF_TYPE Type, std::string Desc);
	void addFormat(std::string & Line);
	void addSample(std::string & Name);
	void updateInfo(TVcfLine & line, std::string & Id, std::string & Data);
	void updatePL(TVcfLine & tempLine, std::string & Data, unsigned int & sample);
	void addToInfo(TVcfLine & line, std::string & Id, std::string & Data);
	int getSampleNum(std::string & Name);
	std::string getSampleName(unsigned int & sample);
	int getNumSamples();
	void checkSampleNum(TVcfLine & line, unsigned int & sample);

	//get variant info
	std::string getChr(TVcfLine & line);
	uint64_t getPos(TVcfLine & line);
	int getNumAlleles(TVcfLine & line);
	std::string getRefAllele(TVcfLine & line);
	std::string getFirstAltAllele(TVcfLine & line);
	std::string getAllele(TVcfLine & line, int num);
	bool variantQualityIsMissing(TVcfLine & line);
	double variantQuality(TVcfLine & line);

	//modify samples
	void addInfoToSample(TVcfLine & line, unsigned int & sample, std::string & tag, std::string & Data);
	void setSampleMissing(TVcfLine & line, unsigned int & sample);
	void setSampleHasUndefinedGenotype(TVcfLine & line, unsigned int & sample);
	void updateField(TVcfLine & line, std::string & tag, std::string & Data, unsigned int & sample);

	//retrieve sample info
	bool sampleIsHaploid(TVcfLine & line, unsigned int & sample);
	bool sampleIsDiploid(TVcfLine & line, unsigned int & sample);
	bool sampleIsHomoRef(TVcfLine & line, unsigned int & sample);
	bool sampleIsHeteroRefNonref(TVcfLine & line, unsigned int & sample);
	std::string getFirstAlleleOfSample(TVcfLine & line, const unsigned int & sample);
	std::string getSecondAlleleOfSample(TVcfLine & line, const unsigned int & sample);
	//std::string sampleGenotype(TVcfLine & line, unsigned int & sample);
	short sampleGenotype(TVcfLine & line, const unsigned int & sample);
	bool sampleIsMissing(TVcfLine & line, unsigned int & sample);
	bool sampleHasUndefinedGenotype(TVcfLine & line, unsigned int & s);
	float sampleGenotypeQuality(TVcfLine & line, unsigned int & sample);
	GTLikelihoods genotypeLikelihoods(TVcfLine & line, unsigned int & sample);
	GTLikelihoods genotypeLikelihoodsPhred(TVcfLine & line, unsigned int & sample);
	void fillGenotypeLikelihoods(TVcfLine & line, unsigned int & s, float* gtl);
	void fillPhredScore(TVcfLine & line, unsigned int & s, uint8_t & gtl_0, uint8_t & gtl_1, uint8_t & gtl_2);
	void fillLog10GenotypeLikelihoods(TVcfLine & line, unsigned int & s, double & gtl_0, double & gtl_1, double & gtl_2);
	std::string sampleContentAt(TVcfLine & line, std::string & tag, unsigned int & sample);
	std::string sampleContentAtNoCheckForMissingSample(TVcfLine & line, std::string & tag, unsigned int & sample);
	int phred(double x);
	float dePhred(double x);

	//output
	void writeColumnDescriptionHeader(std::ostream & out);
	void writeInfoHeader(std::ostream &  out);
	void writeFormatHeader(std::ostream &  out);
	void writeLine(TVcfLine & line, std::ostream & out);
};

}; //end namespace

#endif /* TVCFPARSER_H_ */
