/*
 * TAlignmentParser.h
 *
 *  Created on: Nov 14, 2017
 *      Author: wegmannd
 */

#ifndef TALIGNMENTPARSER_H_
#define TALIGNMENTPARSER_H_

#include "bamtools/api/BamReader.h"
#include "bamtools/api/BamWriter.h"
#include "bamtools/api/SamSequenceDictionary.h"
#include "bamtools/utils/bamtools_fasta.h"
#include "TGenotypeMap.h"
#include "TReadGroups.h"
#include "TLog.h"
#include "TRecalibration.h"
#include "TPostMortemDamage.h"


//-----------------------------------------------------
//TFastaBuffer
//-----------------------------------------------------
//a buffer class to speed up adding the reference sequence to each read
//This class makes use of the fact that bam files are sorted, hence the buffer can always start at the current position

class TFastaBuffer{
private:
	BamTools::Fasta* reference;
	int bufferSize;
	std::string referenceSequence;

	int curChr;
	long curStart, curEnd;

	void moveTo(const int & chr, const int32_t & pos);

public:
	TFastaBuffer(BamTools::Fasta* Reference);
	void fill(const int & chr, const int32_t & start, const int32_t end, std::string & ref);
};


//-----------------------------------------------------
//TAlignmentParser
//-----------------------------------------------------
class TAlignmentParser{
private:
	//variables
	unsigned int maxSize;
	TGenotypeMap genoMap;
	TQualityMap qualityMap;
	TReadGroups* readGroupTable;
	TLog* logfile;
	bool _keepDuplicates;

	//quality filter
	bool applyQualityFilter;
	int minQual, maxQual;
	int minQualForPrinting, maxQualForPrinting;

	//read trimming
	bool trimReads;
	int trimmingLength3Prime;
	int trimmingLength5Prime;

	//data
	bool initialized;
	bool parsed;
	bool changed;
	int* qualityOriginal; //Note: quality is char as int: quality = (int) bam.quality
	int* qualityRecalibrated;
	double* errorRates;

	//tmp variables
	std::vector<BamTools::CigarOp>::const_iterator cigarIter;
	std::vector<BamTools::CigarOp>::const_iterator cigarEnd;
	std::string tmpString;
	std::string tmpString2;
	std::string referenceSequence;
	std::string::iterator stringIt;

	//reference
	bool hasReference;
	TFastaBuffer* fastaBuffer;

	//functions
	inline int toQual(const char & q);
	void parseBasesQualities();
	void setDistancesFromEnds();
	void filterForBaseQuality();
	void filterForPrintingBaseQuality(std::string & qual);
	void trimRead();
	void fillContext();
	void fillPmdProbabilities(TPMD* pmdObjects);

public:
	//alignment: goal is to make this private!
	BamTools::BamAlignment bamAlignment;

	//details
	int length;
	int chrNumber;
	std::string readGroup;
	int readGroupId;

	int32_t position;
	bool isReverseStrand;
	bool isProperPair;
	bool passedFilters;
	bool recalibrated;

	//per base data
	//TODO: try to move to private if possible
	Base* base;
	char* baseAsChar; //TODO: to be removed, if possible
	BaseContext* context;
	int* quality; //pointer to qualities to be used
	bool* aligned; //whether or not base is aligned to ref. Insertions are not aligned
	int* alignedPos;
	int* distFrom3Prime;
	int* distFrom5Prime;
	double* pmdCT;
	double* pmdGA;

	//soft clipped data
	uint8_t softClippedEntry; //0 means start, 1 means end of read
	int* softClippedLength;
	char** softClippedBase;
	char** softClippedQuality;

	//construction
	TAlignmentParser();
	TAlignmentParser(TReadGroups* readGroupTable, unsigned int MaxSize, TLog* Logfile);
	~TAlignmentParser(){
		clear();
		if(hasReference)
			delete fastaBuffer;
	}
	void init(TReadGroups* readGroupTable, unsigned int MaxSize, TLog* Logfile);
	void clear();
	void keepDuplicates(){_keepDuplicates = true;};
	void setQualityFilters(int minQual, int maxQual);
	void setQualityRangeForPrinting(int minQual, int maxQual);
	void setReadTrimming(int trim3Prime, int trim5Prime);

	//functions to read and parse
	bool readAlignment(BamTools::BamReader & bamReader);
	void parse();
	void addReference(BamTools::Fasta* reference);
	void fillReferenceSequence();

	//functions to access and modify data
	std::string& name(){return bamAlignment.Filename;};
	void recalibrate(TRecalibration & recalObject);
	void recalibrate(TRecalibration & recalObject, TPMD* pmdObjects);
	void binQualityScores();
	void addToPMDTables(TPMDTables & pmdTables);
	double calculatePMDS(double & pi, TPMD* pmdObjects);
	void assessSoftClipping(int & S_left, int & middle, int & S_right);
	void addToQualityTable(TQualityTable & qualTable);

	//functions to modify alignment
	void updateOptionalSamField(std::string tag, float value);
	void downsampleAlignment(double& fraction, TRandomGenerator& randomGenerator);

	//functions to write / print alignment
	void save(BamTools::BamWriter & bamWriter);
	void print();
};



#endif /* TALIGNMENTPARSER_H_ */
