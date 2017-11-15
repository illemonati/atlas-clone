/*
 * TAlignmentParser.cpp
 *
 *  Created on: Nov 14, 2017
 *      Author: wegmannd
 */

#include "TAlignmentParser.h"

TAlignmentParser::TAlignmentParser(){
	maxSize = 0;
	readGroupTable = NULL;

	//data
	initialized = false;
	length = 0;
	readGroupId = -1;
	position = 0;
	isReverseStrand = false;
	base = NULL;
	baseAsChar = NULL;
	context = NULL;
	quality = NULL;
	aligned = NULL;
	alignedPos = NULL;
	distFrom3Prime = NULL;
	distFrom5Prime = NULL;

	//tmp variables
	i = 0;
	d = 0;
	k = 0;
	p = 0;
}

TAlignmentParser::TAlignmentParser(TReadGroups* ReadGroupTable, int MaxSize){
	TAlignmentParser();
	init(ReadGroupTable, MaxSize);
};

void TAlignmentParser::init(TReadGroups* ReadGroupTable, int MaxSize){
	clear();

	maxSize = MaxSize;
	readGroupTable = ReadGroupTable;

	base = new Base[maxSize];
	baseAsChar = new char[maxSize];
	context = new BaseContext[maxSize];
	quality = new int[maxSize];
	aligned = new bool[maxSize];
	alignedPos = new int[maxSize];
	distFrom3Prime = new int[maxSize];
	distFrom5Prime = new int[maxSize];

	initialized = true;
};

void TAlignmentParser::clear(){
	if(initialized){
		delete[] base;
		delete[] quality;
		delete[] aligned;
		delete[] alignedPos;
		delete[] distFrom3Prime;
		delete[] distFrom5Prime;
		initialized = false;
	}
};

inline int TAlignmentParser::toQual(const char & q){
	//TODO: switch to already parsing quality here. Will need to adjust filters in TGenome!!!!
	//return (int) q - 33;
	return (int) q;
};

void TAlignmentParser::parseBasesQualities(BamTools::BamAlignment & bamAlignment){
	// iterate over CigarOps
	d = 0; //index regarding data structures
	k = 0; //index inside read
	p = 0; //index regarding reference position (!= k for indels)

	std::vector<BamTools::CigarOp>::const_iterator cigarIter = bamAlignment.CigarData.begin();
	std::vector<BamTools::CigarOp>::const_iterator cigarEnd  = bamAlignment.CigarData.end();

	for ( ; cigarIter != cigarEnd; ++cigarIter ) {
		const BamTools::CigarOp& op = (*cigarIter);
		switch ( op.Type ) {

			// for 'M', '=' or 'X': just copy
			case (BamTools::Constants::BAM_CIGAR_MATCH_CHAR)    :
			case (BamTools::Constants::BAM_CIGAR_SEQMATCH_CHAR) :
			case (BamTools::Constants::BAM_CIGAR_MISMATCH_CHAR) :
				for(i=0; i<op.Length; ++i, ++d, ++k, ++p){
					base[d] = genoMap.getBaseOnlyCapitals(bamAlignment.QueryBases[k]);
					baseAsChar[d] = bamAlignment.QueryBases[k];
					quality[d] = toQual(bamAlignment.Qualities[k]);
					aligned[d] = true;
					alignedPos[d] = p;
				}
				break;

			//for 'S' - soft clip: ignore bases, but increase k
			case (BamTools::Constants::BAM_CIGAR_SOFTCLIP_CHAR) :
				k += op.Length;
				break;

				for(i=0; i<op.Length; ++i, ++d, ++k, ++p){
					base[d] = genoMap.getBaseOnlyCapitals(bamAlignment.QueryBases[k]);
					baseAsChar[d] = bamAlignment.QueryBases[k];
					quality[d] = toQual(bamAlignment.Qualities[k]);
					aligned[d] = false;
					alignedPos[d] = p;
				}
				break;


			//for 'I' - insertion: copy bases, but put aligned pos to
			case (BamTools::Constants::BAM_CIGAR_INS_CHAR)      :
				for(i=0; i<op.Length; ++i, ++d, ++k){
					base[d] = genoMap.getBaseOnlyCapitals(bamAlignment.QueryBases[k]);
					baseAsChar[d] = bamAlignment.QueryBases[k];
					quality[d] = toQual(bamAlignment.Qualities[k]);
					aligned[d] = false;
					alignedPos[d] = -1;
				}
				break;

			// for 'D' - deletion: just add to postion
			case (BamTools::Constants::BAM_CIGAR_DEL_CHAR) :
				p += op.Length;
				break;

			// for 'N' - skipped region: copy but say that bases were not aligned
			case (BamTools::Constants::BAM_CIGAR_REFSKIP_CHAR) :
				for(i=0; i<op.Length; ++i, ++d, ++k, ++p){
					base[d] = genoMap.getBaseOnlyCapitals(bamAlignment.QueryBases[k]);
					baseAsChar[d] = bamAlignment.QueryBases[k];
					quality[d] = toQual(bamAlignment.Qualities[k]);
					aligned[d] = false;
					alignedPos[d] = p;
				}
				break;

			// for 'H' - hard clip: do nothing as these bases are not present in SEQ
			case (BamTools::Constants::BAM_CIGAR_HARDCLIP_CHAR) :
				break;

			// invalid CIGAR op-code
			default:
				throw (std::string) "CIGAR operation type '" + op.Type + "' not supported!";
		}
	}

	//update length
	length = k;

	if(length != bamAlignment.Length)
		throw "Length mismatch!";
};

void TAlignmentParser::setDistancesFromEnds(BamTools::BamAlignment & bamAlignment){
	//first parse bases and qualities
	parseBasesQualities(bamAlignment);

	//now calculate distances from 3 and 5 prime ends

	//is it paired-end?
	if(bamAlignment.IsProperPair()){
		if(bamAlignment.IsReverseStrand()){
			//reverse (can be either first or second mate, but it's the one that comes second in bam file)
			//hence distance from 3' is given by f(dist since beginning of fragment) = f(insert - len + pos)
			//and distance from 5' is given as f(end of fragment) = f(len - pos - 1)
			k = abs(bamAlignment.InsertSize) - length;
			p = length - 1;
			for(d=0; d<length; ++d){
				distFrom3Prime[d] = k + d;
				distFrom5Prime[d] = p - d;
			}
		} else {
			//forward (can be either first or second mate, but it's the one that comes first in bam file)
			//Hence distance from 3' is given as a function of pos
			//And distance from 5' is given by (length of fragment) - pos -1
			//NOTE! we ignore indels when calculating distance from 5' since we can not know this info.
			//Luckily, this has only minimal effect since these distances are far from fragment ends
			p = abs(bamAlignment.InsertSize) - 1;
			for(d=0; d<length; ++d){
				distFrom3Prime[d] = d;
				distFrom5Prime[d] = p - d;
			}
		}
	} else {
		//treat as single end
		p = length - 1;
		if(bamAlignment.IsReverseStrand()){
			//not in pair & reverse
			//Hence distance from 5' is just pos
			//And distance from 3' is just len - pos - 1
			for(d=0; d<length; ++d){
				distFrom3Prime[d] = p - d;
				distFrom5Prime[d] = d;
			}

		} else {
			//not in pair & forward
			//Hence distance from 3' is given as a function of pos
			//And distance from 5' is given by len - pos - 1
			for(d=0; d<length; ++d){
				distFrom3Prime[d] = d;
				distFrom5Prime[d] = p - d;
			}
		}
	}
};

void TAlignmentParser::fillContext(){
	if(isReverseStrand){
		//reverse
		for(d=0; d<(length-1); ++d)
			context[d] = genoMap.getContextReverseRead(base[d+1], base[d]);
		context[d] = genoMap.getContext(N, base[d]);
	} else {
		//forward
		context[0] = genoMap.getContext(N, base[0]);
		for(d=1; d<length; ++d)
			context[d] = genoMap.getContext(base[d-1], base[d]);
	}
};

void TAlignmentParser::parse(BamTools::BamAlignment & bamAlignment){
	//add basic info
	position = bamAlignment.Position;
	name = bamAlignment.Name; //to be removed, if possible
	isReverseStrand = bamAlignment.IsReverseStrand();

	//Extract Read Group Info
	bamAlignment.GetTag("RG", readGroup);
	readGroupId = readGroupTable->find(readGroup);

	//first parse bases and qualities
	parseBasesQualities(bamAlignment);

	//then update distances from ends
	setDistancesFromEnds(bamAlignment);
};

void TAlignmentParser::print(){
	std::cout << "NAME:\t" << name << std::endl;
	std::cout << "LEN:\t" << length << std::endl;

	//print bases
	std::cout << "SEQ:\t";
	for(d=0; d<length; ++d)
		std::cout << genoMap.getBaseAsChar(base[d]);
	std::cout << std::endl;

	//print qualities
	std::cout << "QUAL:\t";
	for(d=0; d<length; ++d)
		std::cout << (char) (quality[d]+33);
	std::cout << std::endl;

	//print aligned pos
	std::cout << "POS:\t";
	for(d=0; d<length; ++d){
		if(d>0) std::cout << ",";
		if(aligned[d])
			std::cout << alignedPos[d];
		else
			std::cout << "-";
	}
	std::cout << std::endl;

	//print dist from 3'
	std::cout << "dist 3':\t";
	for(d=0; d<length; ++d){
		if(d>0) std::cout << ",";
		std::cout << distFrom3Prime[d];
	}
	std::cout << std::endl;

	//print dist from 5'
	std::cout << "dist 5':\t";
	for(d=0; d<length; ++d){
		if(d>0) std::cout << ",";
		std::cout << distFrom5Prime[d];
	}
	std::cout << std::endl;
}




