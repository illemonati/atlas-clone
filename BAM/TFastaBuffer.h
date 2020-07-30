/*
 * TFastaBuffer.h
 *
 *  Created on: Feb 10, 2020
 *      Author: wegmannd
 */

#ifndef TFASTABUFFER_H_
#define TFASTABUFFER_H_

#include "../bamtools/utils/bamtools_fasta.h"
#include "TGenotypeMap.h"
#include "TGenomePosition.h"

namespace BAM{

//-----------------------------------------------------
//TFastaBuffer
//-----------------------------------------------------
//a buffer class to speed up adding the reference sequence to each read
//This class makes use of the fact that bam files are sorted, hence the buffer can always start at the current position

class TFastaBuffer{
private:
	TGenotypeMap* _genoMap;
	bool _hasReference;

	//tmp storage (mutable)
	mutable BamTools::Fasta _reference;
	mutable uint32_t _bufferSize;
	mutable std::string _referenceSequence;
	mutable TGenomeWindow _coordinates;

	void _moveTo(const BAM::TGenomePosition Position) const;

public:
	TFastaBuffer(){
		_bufferSize = 100000;
		_referenceSequence = "";
		_hasReference = false;
		_genoMap = nullptr;
	};

	TFastaBuffer(std::string fastaFile, TGenotypeMap* GenoMap){
		initialize(fastaFile, GenoMap);
	};

	TFastaBuffer(std::string fastaFile, TGenotypeMap* GenoMap, const uint32_t BufferSize){
		initialize(fastaFile, GenoMap, BufferSize);
	};

	~TFastaBuffer(){};

	void initialize(std::string fastaFile, TGenotypeMap* GenoMap, const uint32_t BufferSize=1000000);
	bool hasReference() const{ return _hasReference; };

	void fill(const TGenomeWindow & Window, std::string & ref) const;
	void fill(const TGenomePosition & Position, const uint32_t & Length, std::string & ref) const;
	void fill(const TGenomePosition & Start, const TGenomePosition & End, std::string & ref) const;
	Base refAt(const TGenomePosition Position) const;
	char refCharAt(const TGenomePosition Position) const;

	/*
	void fill(const uint16_t & chr, const uint32_t start, const uint32_t end, std::string & ref);
	char refAt(const uint16_t & chr, const int32_t & position);
	*/

	/*
	int getCurChr(){return _curChr;};
	long getCurStart(){return _curStart;};
	long getCurEnd(){return _curEnd;};
	*/
};

}; //end namespace

#endif /* TFASTABUFFER_H_ */
