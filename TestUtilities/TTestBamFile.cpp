/*
 * TTestBamFile.cpp
 *
 *  Created on: Aug 10, 2020
 *      Author: wegmannd
 */

#include "TTestBamFile.h"

namespace TestUtilities{

//--------------------------------------
// TTestBamFile
//--------------------------------------
TTestBamFile::TTestBamFile(const std::vector<uint32_t> ChrLength, const uint32_t & NumReadGroups){
	_initialize(ChrLength, NumReadGroups);
};

TTestBamFile::TTestBamFile(const std::string & Filename, const std::vector<uint32_t> ChrLength, const uint32_t & NumReadGroups){
	_initialize(ChrLength, NumReadGroups);

	//open BAM file for writing
	_filename = Filename;
	_bamFile.open(_filename, _header, _chromosomes, _readGroups, &_genoMap, &_qualMap);
};

void TTestBamFile::_initialize(const std::vector<uint32_t> ChrLength, const uint32_t & NumReadGroups){
	BAM::TSamHeader header("1.6", "coordinate", "none");
	_initializeChromosomes(ChrLength);
	_initializeReadGroups(NumReadGroups);

	//tmp vars
	_dummySequence =  "AAACCCGGGTTTACGTTGCAAACGTGGCCGTGACACCGTCGACAGGTGCCACACAGTGGCAAATTGGCCGGTGCAAACCAAACCAAGGTTGCCCG";
	_dummyQualities = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
	_dummyMaxLength = _dummySequence.length();
	_dummyBasePos = 0;
	_dummyQualPos = 0;
	_dummyMinLength = 10;
	_dummyMaxLength = 50;
	_dummyLength = _dummyMinLength;
	_dummyReadGroup = 0;
	_dummyIsReverseStrand = false;
};

void TTestBamFile::_initializeChromosomes(const std::vector<uint32_t> ChrLength){
	uint32_t refID = 0;
	for(auto& len : ChrLength){
		++refID;
		_chromosomes.appendChromosome("Chr" + toString(refID), len);
	}
};

void TTestBamFile::_initializeReadGroups(const uint32_t & NumReadGroups){
	for(uint32_t i=0; i<NumReadGroups; ++i){
		const BAM::TReadGroup& rg = _readGroups.add("ReadGroup" + toString(i));
		rg.sequencingCenter_CN = __GLOBAL_APPLICATION_NAME__ + " " + __GLOBAL_APPLICATION_VERSION__;
		rg.description_DS = "Simulated with commit " + __GLOBAL_APPLICATION_COMMIT__;
		rg.sample_SM = "TestSample";
		rg.sequencingTechnology_PL = "ILLUMINA";
	}
};

void TTestBamFile::openOutput(const std::string & Filename){
	//open BAM file for writing
	_filename = Filename;
	_bamFile.open(_filename, _header, _chromosomes, _readGroups, &_genoMap, &_qualMap);
};

void TTestBamFile::closeOutput(){
	_bamFile.close();
};

void TTestBamFile::writeAlignment(const BAM::TAlignment & alignment){
	//store for later comparisons
	_writtenAlignments.emplace_back(alignment);

	//write to BAM
	_bamFile.writeAlignment(alignment);
};

void TTestBamFile::writeDummyAlignment(const BAM::TGenomePosition & position, const uint32_t & length, const uint32_t & readGroup, const bool & isReverseStrand){
	//extract sequence / qualities
	std::string s = _dummySequence.substr(_dummyBasePos, length);
	while(s.length() < length){
		s += _dummySequence.substr(0, length - s.length());
	}
	std::string q = _dummySequence.substr(_dummyQualPos, length);
	while(q.length() < length){
		q += _dummySequence.substr(0, length - q.length());
	}

	//iterate positions
	_dummyBasePos = (_dummyBasePos + 1) % _dummySequence.length();
	_dummyQualPos = (_dummyQualPos  + 3) % _dummyQualities.length();

	//write alignment
	BAM::TAlignment alignment(position);
    alignment.setSequenceQualitiesOnlyMatches(s, q);
    alignment.setReadGroup(readGroup);
    alignment.setIsReverseStrand(isReverseStrand);
    writeAlignment(alignment);

    //update dummy position
    _dummyPosition = alignment;
};

void TTestBamFile::writeDummyAlignment(const BAM::TGenomePosition & position, const uint32_t & length){
	writeDummyAlignment(position, length, _dummyReadGroup, _dummyIsReverseStrand);

	//iterate
	_dummyReadGroup = (_dummyReadGroup + 1) % _readGroups.size();
	_dummyIsReverseStrand = !_dummyIsReverseStrand;
};

void TTestBamFile::writeDummyAlignment(const BAM::TGenomePosition & position){
	writeDummyAlignment(position, _dummyLength);

	_dummyLength += 3;
	if(_dummyLength > _dummyMaxLength){
		_dummyLength = _dummyMinLength + _dummyLength % _dummyMaxLength;
	}
};

void TTestBamFile::writeDummyAlignments(const uint32_t & numAlignments){
	//get distance between alignments
	uint32_t usableLength = _chromosomes.referenceLength() - _chromosomes.size() * _dummyMaxLength;
	uint32_t dist = usableLength / ((double) numAlignments + 1);
	if(dist > _chromosomes.minLength() - _dummyMaxLength){
		dist = _chromosomes.minLength() - _dummyMaxLength;
	}

	BAM::TGenomePosition position;
	auto chr = _chromosomes.begin();

	for(uint32_t i=0; i<numAlignments; ++i){
		//iterate position
		position += dist;
		if(position + _dummyLength > chr->chrEnd){
			++chr;
			if(chr == _chromosomes.end()){
				throw std::runtime_error("void TTestBamFile::writeDummyAlignments(const uint32_t & numAlignments): chromosome reached end!");
			}

			position = chr->chrStart;
		}
		writeDummyAlignment(position);
	}
};


}; //end namespace
