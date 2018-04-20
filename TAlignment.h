/*
 * TAlignment.h
 *
 *  Created on: Apr 19, 2018
 *      Author: wegmannd
 */

#ifndef TALIGNMENT_H_
#define TALIGNMENT_H_

#include "stringFunctions.h"
#include "TBase.h"
#include "BamAlignment.h"


class TAlignment{
private:
	//details
	bool empty;
	int length;
	int chrNumber;
	std::string readGroup;
	int readGroupId;
	int32_t position;
	bool isReverseStrand;
	bool isProperPair;
	bool passedFilters;
	bool recalibrated;

	//data
	BamTools::BamAlignment bamAlignment;
	unsigned int maxSize;
	bool storageInitialized;
	bool parsed;
	bool changed;

	TBase* bases;

	//per base data
/*	Base* base;
	char* baseAsChar; //TODO: to be removed, if possible
	BaseContext* context;
	int* quality; //pointer to qualities to be used
	int* qualityOriginal; //Note: quality is char as int: quality = (int) bam.quality
	int* qualityRecalibrated;
	double* errorRates;
	bool* aligned; //whether or not base is aligned to ref. Insertions are not aligned
	int* alignedPos;
	int* distFrom3Prime;
	int* distFrom5Prime;
	double* pmdCT;
	double* pmdGA;

	//soft clipped data
	int* softClippedLength;
	char** softClippedBase;
	char** softClippedQuality;
*/
	uint8_t softClippedEntry; //0 means start, 1 means end of read


	//reference
	bool hasReference;
	std::string referenceSequence;

	//functions
	void initStorage();
	void freeStorage();
	void clear();

	//functions to read and parse
	void parse(TGenotypeMap & genoMap, TQualityMap & qualityMap);
	void setDistancesFromEnds();
	void parseBasesQualities(TGenotypeMap & genoMap, TQualityMap & qualityMap);
	void fillContext(TGenotypeMap & genoMap);
	void fillPmdProbabilities(TPMD* pmdObjects);


	//functions to access and modify data
	std::string& name(){return bamAlignment.Filename;};
	void filterForBaseQuality(int & minQual, int & maxQual);
	void filterForPrintingBaseQuality(std::string & qual, int & minQualForPrinting, int & maxQualForPrinting);
	void trimRead(int & trimmingLength3Prime, int & trimmingLength5Prime);
	void recalibrate(TRecalibration & recalObject, TQualityMap & qualityMap);
	void recalibrate(TRecalibration & recalObject, TPMD* pmdObjects, TFastaBuffer* fastaBuffer, TQualityMap & qualityMap);
	void binQualityScores(TQualityMap & qualityMap);
	void addToPMDTables(TPMDTables & pmdTables, TFastaBuffer* fastaBuffer, TGenotypeMap & genoMap);
	double calculatePMDS(double & pi, TPMD* pmdObjects, TFastaBuffer* fastaBuffer);
	void assessSoftClipping(int & S_left, int & middle, int & S_right);
	void addToQualityTable(TQualityTable & qualTable);

	//functions to modify alignment
	void updateOptionalSamField(std::string tag, float value);
	void downsampleAlignment(double& fraction, TRandomGenerator& randomGenerator);


public:
	TAlignment();
	TAlignment(unsigned int MaxSize);
	~TAlignment(){
		freeStorage();
	}

	void fill(BamTools::BamAlignment & bamAlignment, int ReadGroupId);
	void setFiltersPassed(bool passed);
	void setReferenceAdded();

	//functions to write / print alignment
	void save(BamTools::BamWriter & bamWriter, TGenotypeMap & genoMap, int & minQualForPrinting, int & maxQualForPrinting);
	void print(TGenotypeMap & genoMap);
};

#endif /* TALIGNMENT_H_ */
