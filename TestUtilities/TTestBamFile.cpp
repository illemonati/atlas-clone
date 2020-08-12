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
	_header.set("1.6", "coordinate", "none", "none");
	_initializeChromosomes(ChrLength);
	_initializeReadGroups(NumReadGroups);

	//tmp vars
	_dummySequence =  "AAACCCGGGTTTACGTTGCAAACGTGGCCGTGACACCGTCGACAGGTGCCACACAGTGGCAAATTGGCCGGTGCAAACCAAACCAAGGTTGCCCG";
	_dummyQualities = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
	_dummyMaxLength = _dummySequence.length();
	_dummyBasePos = 0;
	_dummyQualPos = 0;
	_dummyMinLength = 10;
	_dummyMaxLength = 50;
	_dummyLength = _dummyMinLength;
	_dummyReadGroup = 0;
	_dummyIsReverseStrand = false;
	_dummyCigarChars = "MID";
	_dummyCigarPos = 0;
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

void TTestBamFile::_iterateReadGroupAndReverseStrand(){
    _dummyReadGroup = (_dummyReadGroup + 1) % _readGroups.size();
    _dummyIsReverseStrand = !_dummyIsReverseStrand;
}

void TTestBamFile::_iterateCigar(BAM::TCigar & cigar, uint32_t length){
    std::string s = _dummyCigarChars.substr(_dummyCigarPos, 3);
    while(s.length() < 3){
        s += _dummyCigarChars.substr(0, 3 - s.length());
    }
    // construct cigar: choose number of S, M, I and D
    cigar.add('S', std::floor(length / 8));
    for (char const &c: s) {
        cigar.add(c, std::floor(length / 3));
    }
    if (length > cigar.lengthRead())
        cigar.add('S',  length - cigar.lengthRead()); // fill rest with S (right)

    // update
    _dummyCigarPos = (_dummyCigarPos + 1) % _dummyCigarChars.size();
}

void TTestBamFile::_iterateLength(){
    _dummyLength += 3;
    if(_dummyLength > _dummyMaxLength){
        _dummyLength = _dummyMinLength + _dummyLength % _dummyMaxLength;
    }
}

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

void TTestBamFile::writeDummyAlignment(const BAM::TGenomePosition & position, const BAM::TCigar & cigar, const uint32_t & readGroup, const bool & isReverseStrand){
	//extract sequence / qualities
	std::string s = _dummySequence.substr(_dummyBasePos, cigar.lengthRead());
	while(s.length() < cigar.lengthRead()){
		s += _dummySequence.substr(0, cigar.lengthRead() - s.length());
	}
	std::string q = _dummyQualities.substr(_dummyQualPos, cigar.lengthRead());
	while(q.length() < cigar.lengthRead()){
		q += _dummyQualities.substr(0, cigar.lengthRead() - q.length());
	}

	//iterate positions
	_dummyBasePos = (_dummyBasePos + 1) % _dummySequence.length();
	_dummyQualPos = (_dummyQualPos  + 3) % _dummyQualities.length();

	//write alignment
	BAM::TAlignment alignment(position);
    alignment.setSequenceQualities(cigar, s, q);
    alignment.setReadGroup(readGroup);
    alignment.setIsReverseStrand(isReverseStrand);
    writeAlignment(alignment);
};

void TTestBamFile::writeDummyAlignment(const BAM::TGenomePosition & position, const BAM::TCigar & cigar){
    writeDummyAlignment(position, cigar, _dummyReadGroup, _dummyIsReverseStrand);

    //iterate
    _iterateReadGroupAndReverseStrand();
};

void TTestBamFile::writeDummyAlignment(const BAM::TGenomePosition & position, const uint32_t & length){
    // iterate over order of M, I and D
    BAM::TCigar cigar;
    _iterateCigar(cigar, length);

    // now write
	writeDummyAlignment(position, cigar);
};

void TTestBamFile::writeDummyAlignment(const BAM::TGenomePosition & position){
	writeDummyAlignment(position, _dummyLength);

	_iterateLength();
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
}

void TTestBamFile::writeDummyAlignment(const char &oneBase, const char &oneQual, const BAM::TGenomePosition &position,
                                  const BAM::TCigar &cigar, const uint32_t &readGroup,
                                  const bool &isReverseStrand) {
    // expand oneBase to string
    std::string seq(cigar.lengthRead(), oneBase);
    std::string qual(cigar.lengthRead(), oneQual);

    //write alignment
    BAM::TAlignment alignment(position);
    alignment.setSequenceQualities(cigar, seq, qual);
    alignment.setReadGroup(readGroup);
    alignment.setIsReverseStrand(isReverseStrand);
    writeAlignment(alignment);


    /*
    std::cout << "----------------------------------" << std::endl;
    std::cout << "refID = " << position.refID() << std::endl;
    std::cout << "position = " << position.position() << std::endl;
    std::cout << "length = " << cigar.lengthRead() << std::endl;
    std::cout << "seq = " << seq << std::endl;
    std::cout << "qual = " << qual << std::endl;
    std::string cigg;
    for (auto & it : cigar){
        std::string cigger(it.length, it.type);
        cigg += cigger;
    }
    std::cout << "cigar = " << cigg << std::endl;
    std::cout << "readGroup = " << (int) readGroup << std::endl;
    std::cout << "isReverseStrand = " << isReverseStrand << std::endl;
    */
};


void TTestBamFile::writeDummyAlignment(const char& oneBase, const char& oneQual, const BAM::TGenomePosition & position, const BAM::TCigar & cigar){
    writeDummyAlignment(oneBase, oneQual, position, cigar, _dummyReadGroup, _dummyIsReverseStrand);

    //iterate
    _iterateReadGroupAndReverseStrand();
}

void TTestBamFile::writeDummyAlignment(const char &oneBase, const char &oneQual, const BAM::TGenomePosition &position,
                                  const uint32_t &length) {
    // construct cigar
    BAM::TCigar cigar;
    _iterateCigar(cigar, length);

    writeDummyAlignment(oneBase, oneQual, position, cigar);
};

void TTestBamFile::writeDummyAlignment(const char &oneBase, const char &oneQual, const BAM::TGenomePosition &position) {
    writeDummyAlignment(oneBase, oneQual, position, _dummyLength);

    // iterate length
    _iterateLength();
};


}; //end namespace
