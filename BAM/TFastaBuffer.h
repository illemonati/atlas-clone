/*
 * TFastaBuffer.h
 *
 *  Created on: Feb 10, 2020
 *      Author: wegmannd
 */

#ifndef TFASTABUFFER_H_
#define TFASTABUFFER_H_

#include <stddef.h>
#include <stdint.h>
#include <string>
#include <vector>
#include "genometools/GenotypeTypes.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "utils/bamtools_fasta.h"

namespace BAM{

//-----------------------------------------------------
//TFastaBuffer
//-----------------------------------------------------
//a buffer class to speed up adding the reference sequence to each read
//This class makes use of the fact that bam files are sorted, hence the buffer can always start at the current position

class TFastaBuffer{
private:
	bool _hasReference;

	//tmp storage (mutable)
	mutable BamTools::Fasta _reference;
	mutable uint32_t _bufferSize;
	mutable std::string _referenceSequence;
	mutable genometools::TGenomeWindow _coordinates;

	void _moveTo(const genometools::TGenomePosition Position) const;
	void _fill(std::vector<genometools::Base> & VecToFill, size_t Length, size_t Offset) const;

public:
	TFastaBuffer(){
		_bufferSize = 100000;
		_referenceSequence = "";
		_hasReference = false;
	};

	TFastaBuffer(std::string fastaFile){
		initialize(fastaFile);
	};

	TFastaBuffer(std::string fastaFile, const uint32_t BufferSize){
		initialize(fastaFile, BufferSize);
	};

	~TFastaBuffer(){};

	void initialize(std::string fastaFile, const uint32_t BufferSize=1000000);
	bool hasReference() const { return _hasReference; };
	explicit operator bool() const { return _hasReference; }

	void fill(const genometools::TGenomeWindow & Window, std::vector<genometools::Base> & VecToFill) const;
	void fill(const genometools::TGenomePosition & Position, uint32_t Length, std::vector<genometools::Base> & VecToFill) const;
	void fill(const genometools::TGenomePosition & Start, const genometools::TGenomePosition & End, std::vector<genometools::Base> & VecToFill) const;
	genometools::Base refAt(const genometools::TGenomePosition Position) const;
	char refCharAt(const genometools::TGenomePosition Position) const;

	/*
	void fill(uint16_t chr, const uint32_t start, const uint32_t end, std::string & ref);
	char refAt(uint16_t chr, const int32_t & position);
	*/

	/*
	int getCurChr(){return _curChr;};
	long getCurStart(){return _curStart;};
	long getCurEnd(){return _curEnd;};
	*/
};

}; //end namespace

#endif /* TFASTABUFFER_H_ */
